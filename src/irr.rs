//! In-repeat-read (IRR) detection heuristic: decides whether a read's
//! sequence is dominated by repetitions of some short motif, and if so,
//! returns the motif's canonical repeat unit.
//!
//! For a candidate period, the repeat unit is the plain per-phase
//! majority-vote base (see [`extract_consensus_base`]), optionally reduced
//! to a smaller period if the extracted unit is itself an exact
//! smaller-period repeat, then canonicalized (lexicographically-smallest
//! rotation, forward or reverse-complement) and scored against the whole
//! read allowing both orientations, with a per-base match/mismatch/
//! low-quality scoring scheme (see [`score_chunk`]). A period only reaches
//! that scoring step once its raw literal match frequency
//! ([`match_frequency_at_offset`]) clears [`MIN_UNIT_FREQUENCY`], and only
//! survives once its final score clears [`MIN_IRR_SCORE`].
//!
//! Every period in `[motif_min_len, motif_max_len]` is evaluated
//! independently against these two thresholds, and every canonical motif
//! that passes is returned -- a read can plausibly satisfy more than one
//! motif (e.g. a longer compound period that is itself a repetition of a
//! shorter one also scores well), so this doesn't collapse to a single
//! "best" period the way a greedy shortest-period search would.
//!
//! `bases` are expected to be uppercase decoded read sequence bytes (as
//! returned by `rust_htslib::bam::record::Seq::as_bytes`), and `quals` are
//! raw (non-ASCII-offset) PHRED scores (as returned by `Record::qual`).

/// The raw (literal, position-by-position) match frequency a candidate
/// period's offset must clear before its consensus repeat unit is even
/// extracted and scored -- see [`match_frequency_at_offset`].
const MIN_UNIT_FREQUENCY: f64 = 0.8;
/// The final quality-aware match score (as a fraction of read length) a
/// candidate motif must clear to be accepted -- see [`score_chunk`].
const MIN_IRR_SCORE: f64 = 0.90;
/// Below this PHRED quality, a mismatch is treated as an uncertain
/// low-quality mismatch rather than a confident one -- see
/// [`score_chunk`].
const MIN_BASE_QUALITY: u8 = 20;

/// Default shortest repeat-unit (motif) length to consider.
pub const DEFAULT_MOTIF_MIN_LEN: u32 = 2;
/// Default longest repeat-unit (motif) length to consider.
pub const DEFAULT_MOTIF_MAX_LEN: u32 = 20;

/// Returns every canonical repeat unit (motif) that independently clears
/// both the unit-frequency and IRR-score thresholds for this read. A read
/// can plausibly satisfy more than one motif (e.g. a longer compound period
/// that is itself a repetition of a shorter one also scores well), so this
/// doesn't collapse to a single "best" motif (see the module docs). Each
/// returned motif is distinct; order is not significant.
pub fn identify_repeat_motifs(
    bases: &[u8],
    quals: &[u8],
    motif_min_len: u32,
    motif_max_len: u32,
) -> Vec<Vec<u8>> {
    let smallest_period = motif_min_len.max(1) as usize;
    let largest_period = (motif_max_len as usize).min(bases.len() / 2 + 1);
    let bases_len = bases.len() as f64;
    // An internal short-circuit threshold only -- every per-motif accept/
    // reject decision below still uses the exact original `score /
    // bases_len >= MIN_IRR_SCORE` comparison, so this never changes which
    // motifs come back (including the bases_len == 0 edge case, where this
    // threshold is 0.0 but the real comparison's `0.0 / 0.0` is NaN and
    // therefore never passes).
    let early_exit_score = MIN_IRR_SCORE * bases_len;

    // The read's reverse complement doesn't depend on the candidate period,
    // so it's the same for every iteration below; computed at most once,
    // lazily, since a period whose forward orientation alone already clears
    // the score threshold never needs it at all.
    let mut bases_rc_cache: Option<(Vec<u8>, Vec<u8>)> = None;

    // At most one motif per candidate period, so this upper bound avoids
    // any reallocation as `motifs` grows.
    let mut motifs: Vec<Vec<u8>> =
        Vec::with_capacity(largest_period.saturating_sub(smallest_period) + 1);
    for period in smallest_period..=largest_period {
        if match_frequency_at_offset(period, bases) < MIN_UNIT_FREQUENCY {
            continue;
        }

        let mut unit = extract_consensus_repeat_unit(period, bases);

        // Attempt to reduce the motif to a smaller period if one exists,
        // exactly matching a perfect (frequency == 1.0) repeat -- this
        // reduction step always searches a fixed [1, 20] range, regardless
        // of the caller's own motif length bounds.
        const PERFECT_MATCH_FREQUENCY: f64 = 1.0;
        const REDUCTION_MIN_LEN: u32 = 1;
        const REDUCTION_MAX_LEN: u32 = 20;
        if let Some(reduced_period) = smallest_frequent_period(
            PERFECT_MATCH_FREQUENCY,
            &unit,
            REDUCTION_MIN_LEN,
            REDUCTION_MAX_LEN,
        ) && reduced_period != period
        {
            unit = extract_consensus_repeat_unit(reduced_period, &unit);
        }

        if unit.len() < motif_min_len as usize || unit.len() > motif_max_len as usize {
            continue;
        }

        let canonical = compute_canonical_repeat_unit(&unit);
        // An empty or literal "N" unit (an all-no-call homopolymer) is
        // never a meaningful repeat motif.
        if canonical.is_empty() || canonical == b"N" || motifs.contains(&canonical) {
            continue;
        }

        let forward_score =
            best_score_across_shifts(&canonical, bases, quals, MIN_BASE_QUALITY, early_exit_score);
        // The forward score alone already clears the threshold, so the
        // (comparatively expensive) reverse-complement orientation can't
        // change the accept/reject outcome -- skip computing it.
        let score = if forward_score >= early_exit_score {
            forward_score
        } else {
            let (bases_rc, quals_rc) = bases_rc_cache.get_or_insert_with(|| {
                let bases_rc = reverse_complement(bases);
                let mut quals_rc = quals.to_vec();
                quals_rc.reverse();
                (bases_rc, quals_rc)
            });
            let reverse_score = best_score_across_shifts(
                &canonical,
                bases_rc,
                quals_rc,
                MIN_BASE_QUALITY,
                early_exit_score,
            );
            forward_score.max(reverse_score)
        };
        if score / bases_len >= MIN_IRR_SCORE {
            motifs.push(canonical);
        }
    }

    motifs
}

// --- motif detection ---------------------------------------------------

#[inline]
fn max_matches_at_offset(offset: usize, bases: &[u8]) -> usize {
    bases.len().saturating_sub(offset)
}

fn match_frequency_at_offset(offset: usize, bases: &[u8]) -> f64 {
    // A period can be at most half the read length, since we need at
    // least two repetitions to observe a periodic match.
    #[allow(clippy::int_plus_one)]
    if bases.len() / 2 + 1 <= offset {
        return 0.0;
    }

    let max_matches = max_matches_at_offset(offset, bases);
    let num_matches = bases[..max_matches]
        .iter()
        .zip(&bases[offset..])
        .filter(|(a, b)| a == b)
        .count();
    num_matches as f64 / max_matches as f64
}

/// Finds the shortest motif period whose match frequency is at least as
/// good as any longer period's, or `None` if none clears `min_frequency`.
/// Used here only for the perfect-match period-reduction step (see
/// [`identify_repeat_motifs`]).
fn smallest_frequent_period(
    min_frequency: f64,
    bases: &[u8],
    motif_min_len: u32,
    motif_max_len: u32,
) -> Option<usize> {
    let smallest_period = motif_min_len.max(1) as usize;
    let largest_period = (motif_max_len as usize).min(bases.len() / 2 + 1);

    let mut max_match_frequency = min_frequency;
    let mut best_period = None;

    for period in (smallest_period..=largest_period).rev() {
        let frequency = match_frequency_at_offset(period, bases);
        if frequency >= max_match_frequency {
            max_match_frequency = frequency;
            best_period = Some(period);
        }
    }

    best_period
}

/// Number of distinct letters [`ALPHABET_INDEX`] covers (`A`..=`Z`). Every
/// byte [`extract_consensus_base`] ever tallies -- raw read bases (always
/// `A`/`C`/`G`/`T`/`N`) and, during period reduction, a previously-computed
/// consensus unit's bytes -- is an uppercase ASCII letter, so a 26-entry
/// array covers every real input while costing far less to zero and scan
/// per call than a full 256-entry one.
const ALPHABET_SIZE: usize = 26;

/// Vote-tally array index for an uppercase ASCII letter (see
/// [`ALPHABET_SIZE`]).
#[inline]
fn alphabet_index(base: u8) -> usize {
    debug_assert!(
        base.is_ascii_uppercase(),
        "expected an uppercase base byte, got {base}"
    );
    (base - b'A') as usize
}

/// Plain per-phase majority vote: the most frequently observed literal byte
/// at this phase, quality-blind. Ties are broken deterministically by
/// preferring the larger byte value.
fn extract_consensus_base(offset: usize, period: usize, bases: &[u8]) -> u8 {
    let mut counts = [0u32; ALPHABET_SIZE];
    let mut index = offset;
    while index < bases.len() {
        counts[alphabet_index(bases[index])] += 1;
        index += period;
    }

    counts
        .iter()
        .enumerate()
        .filter(|&(_, &count)| count > 0)
        .max_by_key(|&(index, &count)| (count, index))
        .map(|(index, _)| index as u8 + b'A')
        .unwrap_or(b'?')
}

fn extract_consensus_repeat_unit(period: usize, bases: &[u8]) -> Vec<u8> {
    (0..period)
        .map(|offset| extract_consensus_base(offset, period, bases))
        .collect()
}

/// `unit` followed by a second copy of itself, so every cyclic rotation of
/// `unit` is available as a `unit.len()`-wide window
/// `doubled(unit)[offset..offset + unit.len()]` for `offset` in
/// `0..unit.len()` -- one allocation shared across every rotation, instead
/// of materializing each rotation as its own owned buffer.
fn doubled(unit: &[u8]) -> Vec<u8> {
    let mut doubled = unit.to_vec();
    doubled.extend_from_slice(unit);
    doubled
}

fn minimal_unit_under_shift(unit: &[u8]) -> Vec<u8> {
    let len = unit.len();
    let doubled = doubled(unit);
    let best_offset = (0..len)
        .min_by_key(|&offset| &doubled[offset..offset + len])
        .unwrap_or(0);
    doubled[best_offset..best_offset + len].to_vec()
}

/// `A<->T`, `C<->G`; anything else (including a sequencer no-call `N`)
/// complements to `N`.
#[inline]
fn complement_base(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'N',
    }
}

fn reverse_complement(bases: &[u8]) -> Vec<u8> {
    bases
        .iter()
        .rev()
        .map(|&base| complement_base(base))
        .collect()
}

fn compute_canonical_repeat_unit(unit: &[u8]) -> Vec<u8> {
    let minimal = minimal_unit_under_shift(unit);
    let unit_rc = reverse_complement(unit);
    let minimal_rc = minimal_unit_under_shift(&unit_rc);
    if minimal_rc < minimal {
        minimal_rc
    } else {
        minimal
    }
}

// --- quality-aware matching ----------------------------------------------

/// Scores one `unit`-length chunk against `unit`: +1 per matching position,
/// -1 per confident mismatch, or +0.5 for a mismatch whose quality is below
/// `min_baseq` (too uncertain to penalize as confidently wrong).
#[inline]
fn score_chunk(unit: &[u8], bases: &[u8], quals: &[u8], min_baseq: u8) -> f64 {
    const MATCH_SCORE: f64 = 1.0;
    const LOWQUAL_MISMATCH_SCORE: f64 = 0.5;
    const MISMATCH_PENALTY: f64 = -1.0;

    bases
        .iter()
        .zip(quals)
        .zip(unit)
        .map(|((&base, &qual), &unit_base)| {
            if base == unit_base {
                MATCH_SCORE
            } else if qual < min_baseq {
                LOWQUAL_MISMATCH_SCORE
            } else {
                MISMATCH_PENALTY
            }
        })
        .sum()
}

/// Scores the whole read against `unit` repeated end to end: `bases`
/// chunked into `unit.len()`-sized pieces, each scored independently and
/// summed.
fn score_repeat(unit: &[u8], bases: &[u8], quals: &[u8], min_baseq: u8) -> f64 {
    bases
        .chunks(unit.len())
        .zip(quals.chunks(unit.len()))
        .map(|(base_chunk, qual_chunk)| score_chunk(unit, base_chunk, qual_chunk, min_baseq))
        .sum()
}

/// The best [`score_repeat`] across every cyclic rotation of `unit` (e.g.
/// `"ATG"` -> `"ATG"`, `"TGA"`, `"GAT"`), since the read's phase relative to
/// the motif is unknown -- the read's first base isn't necessarily the
/// motif's first base. Rotations are scored as `unit.len()`-wide windows
/// into a single [`doubled`] copy of `unit` rather than each being
/// materialized as its own owned buffer. Stops scanning further rotations
/// as soon as the running best clears `early_exit_at` -- callers only ever
/// compare this result against that same threshold, and once it's cleared,
/// no later rotation (which can only raise, never lower, the running max)
/// could change that outcome.
fn best_score_across_shifts(
    unit: &[u8],
    bases: &[u8],
    quals: &[u8],
    min_baseq: u8,
    early_exit_at: f64,
) -> f64 {
    let len = unit.len();
    let doubled = doubled(unit);

    let mut best = f64::NEG_INFINITY;
    for offset in 0..len {
        let shift = &doubled[offset..offset + len];
        let score = score_repeat(shift, bases, quals, min_baseq);
        if score > best {
            best = score;
        }
        if best >= early_exit_at {
            break;
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn max_matches_at_offset_various() {
        let bases = b"ATCGATCG";
        assert_eq!(max_matches_at_offset(0, bases), 8);
        assert_eq!(max_matches_at_offset(1, bases), 7);
        assert_eq!(max_matches_at_offset(2, bases), 6);
        assert_eq!(max_matches_at_offset(8, bases), 0);
        assert_eq!(max_matches_at_offset(9, bases), 0);
    }

    #[test]
    fn match_frequency_at_offset_various() {
        let bases = b"GGCCCCGGCCCC";
        let expected = [0.73, 0.40, 0.33, 0.25, 0.57, 1.00];
        for offset in 1..=6 {
            let freq = match_frequency_at_offset(offset, bases);
            assert!(
                (freq - expected[offset - 1]).abs() < 0.01,
                "offset {offset}: got {freq}, expected {}",
                expected[offset - 1]
            );
        }
    }

    #[test]
    fn match_frequency_at_offset_imperfect_repeat() {
        let bases = b"ATGATCATGTTGATG";
        let freq = match_frequency_at_offset(3, bases);
        assert!((freq - 8.0 / 12.0).abs() < 1e-9);
    }

    #[test]
    fn smallest_frequent_period_typical() {
        assert_eq!(
            smallest_frequent_period(0.85, b"GGCCCCGGCCCC", 1, 20),
            Some(6)
        );
        assert_eq!(
            smallest_frequent_period(0.85, b"ATGATCATGATGATGATGATG", 1, 20),
            Some(6)
        );
    }

    #[test]
    fn smallest_frequent_period_none_when_no_period_found() {
        assert_eq!(smallest_frequent_period(0.85, b"ATCGGCTA", 1, 20), None);
    }

    #[test]
    fn extract_consensus_base_basic() {
        let bases = b"CGATGACTG";
        assert_eq!(extract_consensus_base(0, 3, bases), b'C');
        assert_eq!(extract_consensus_base(1, 3, bases), b'G');
        assert_eq!(extract_consensus_base(2, 3, bases), b'A');
    }

    #[test]
    fn extract_consensus_base_majority_vote_ignores_minority() {
        // 8 A's and 2 G's at the same phase: plain majority vote, no
        // ambiguity calling, so this resolves to the outright majority (A)
        // regardless of how confident/frequent the minority is.
        let bases: Vec<u8> = (0..10).map(|i| if i < 8 { b'A' } else { b'G' }).collect();
        assert_eq!(extract_consensus_base(0, 1, &bases), b'A');
    }

    #[test]
    fn extract_consensus_repeat_unit_basic() {
        assert_eq!(extract_consensus_repeat_unit(3, b"CGGCGGCGG"), b"CGG");
        assert_eq!(extract_consensus_repeat_unit(3, b"CGGATTATTATTCGG"), b"ATT");
    }

    #[test]
    fn minimal_unit_under_shift_basic() {
        assert_eq!(minimal_unit_under_shift(b"GGC"), b"CGG");
    }

    #[test]
    fn compute_canonical_repeat_unit_basic() {
        assert_eq!(compute_canonical_repeat_unit(b"CGG"), b"CCG");
        assert_eq!(compute_canonical_repeat_unit(b"GCC"), b"CCG");
    }

    /// A read whose composition is 20 A's followed by a single G, repeated:
    /// mostly-A with only rare G interruptions passes the homopolymer "A"
    /// motif's thresholds despite not being a pure homopolymer, while the
    /// exact 21bp repeat unit also passes on its own -- so this read
    /// legitimately qualifies under two distinct, non-harmonic motifs
    /// (unlike a pure homopolymer, where every candidate period reduces
    /// back to the same single-base canonical unit). Every other period
    /// either reduces to the same "A" homopolymer (plain majority vote
    /// always favors the dominant A, with no ambiguity calling to instead
    /// flag it as a mixed position) or fails the raw frequency/score
    /// thresholds outright, so this fixture's two motifs are the only ones
    /// expected back.
    fn mostly_a_with_rare_g() -> Vec<u8> {
        let unit: Vec<u8> = (0..20).map(|_| b'A').chain(std::iter::once(b'G')).collect();
        unit.iter().cloned().cycle().take(21 * 16).collect()
    }

    #[test]
    fn identify_repeat_motifs_returns_single_motif_for_pure_repeat() {
        let bases = "CAG".repeat(20).into_bytes();
        let quals = vec![40u8; bases.len()];
        assert_eq!(
            identify_repeat_motifs(&bases, &quals, 1, 20),
            vec![b"AGC".to_vec()]
        );
    }

    #[test]
    fn identify_repeat_motifs_returns_empty_for_non_repetitive() {
        let bases = b"ACGTTGCAACGGTTCAGTAGCTAGCATCGATCGTAGCTAGGCTAGCATCGTAGCTAGCA";
        let quals = vec![40u8; bases.len()];
        assert!(identify_repeat_motifs(bases, &quals, 1, 20).is_empty());
    }

    #[test]
    fn identify_repeat_motifs_returns_multiple_distinct_motifs() {
        let bases = mostly_a_with_rare_g();
        let quals = vec![40u8; bases.len()];

        let motifs = identify_repeat_motifs(&bases, &quals, 1, 30);

        assert_eq!(
            motifs,
            vec![
                b"A".to_vec(),
                (0..20u8)
                    .map(|_| b'A')
                    .chain(std::iter::once(b'G'))
                    .collect()
            ],
            "expected exactly the mostly-A homopolymer motif and the exact 21bp repeat unit: \
             {motifs:?}"
        );
    }

    #[test]
    fn identify_repeat_motifs_rejects_period_below_min_unit_frequency() {
        // A read with no periodic structure at all: no candidate period's
        // raw match frequency clears MIN_UNIT_FREQUENCY, so nothing is even
        // scored.
        let bases = b"ACGTTGCAACGGTTCAGTAGCTAGCATCGATCGTAGCTAGGCTAGCATCGTAGCTAGCA";
        let quals = vec![40u8; bases.len()];
        assert!(identify_repeat_motifs(bases, &quals, 1, 20).is_empty());
    }

    #[test]
    fn identify_repeat_motifs_rejects_literal_n_unit() {
        // An all-no-call "read": every base is a sequencer no-call, so the
        // consensus unit is literally "N" and must be rejected outright.
        let bases = vec![b'N'; 60];
        let quals = vec![40u8; bases.len()];
        assert!(identify_repeat_motifs(&bases, &quals, 1, 20).is_empty());
    }
}
