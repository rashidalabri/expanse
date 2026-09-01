//! In-repeat-read (IRR) detection heuristic: decides whether a read's
//! sequence is dominated by repetitions of some short motif, and if so,
//! returns the motif's canonical repeat unit.
//!
//! A returned motif may contain IUPAC ambiguity codes (`R`, `Y`, `S`, `W`,
//! `K`, `M`, `B`, `D`, `H`, `V`, `N`) at positions that are consistently
//! mixed across repeat copies rather than a single fixed base -- e.g. `GCN`
//! (an unconstrained 3rd position) or `AARRG` (two purine-only positions).
//! See [`extract_consensus_base_iupac`] for how a position is called
//! degenerate, and [`degenerate_count`] for the per-motif cap callers place
//! on how many positions may be degenerate before a motif is rejected
//! outright -- separately configurable for mononucleotide, dinucleotide,
//! trinucleotide, and longer motifs (see `max_degenerate_mononucleotide`
//! /`max_degenerate_dinucleotide`/`max_degenerate_trinucleotide`
//! /`max_degenerate_other` on
//! [`classify_in_repeat_read`]/[`classify_in_repeat_read_all`]).
//!
//! `bases` are expected to be uppercase decoded read sequence bytes (as
//! returned by `rust_htslib::bam::record::Seq::as_bytes`), and `quals` are
//! raw (non-ASCII-offset) PHRED scores (as returned by `Record::qual`).

const MIN_UNIT_FREQUENCY: f64 = 0.8;
const MIN_IRR_SCORE: f64 = 0.90;
const MIN_BASE_QUALITY: u8 = 20;
const DEFAULT_REDUCTION_MOTIF_RANGE: (u32, u32) = (1, 20);

/// The coverage fraction an IUPAC tier's best code must clear (see
/// [`extract_consensus_base_iupac`]) to be called instead of escalating to
/// a more degenerate tier. Ordinary sequencing error is normally well
/// under 10%, so this comfortably tolerates noise while still catching a
/// genuine, substantial minority allele (e.g. a true ~50/50 split).
const MIN_BASE_PURITY: f64 = 0.9;

/// Minimum number of real (A/C/G/T) observations a phase position needs
/// before degenerate (IUPAC) calling is even attempted; below this, a
/// handful of votes can't distinguish genuine ambiguity from small-sample
/// noise, so [`extract_consensus_base_iupac`] falls back to plain majority
/// vote.
const MIN_DEGENERACY_SAMPLES: u32 = 8;

/// IUPAC ambiguity codes grouped by degeneracy tier (fewest represented
/// bases first), each paired with the literal bases it represents. Used by
/// [`extract_consensus_base_iupac`] to find the smallest code that covers
/// enough of the observed bases at a position, and by [`iupac_matches`] to
/// test whether an observed base satisfies one of these codes.
const IUPAC_TIERS: [&[(u8, &[u8])]; 4] = [
    &[(b'A', b"A"), (b'C', b"C"), (b'G', b"G"), (b'T', b"T")],
    &[
        (b'R', b"AG"),
        (b'Y', b"CT"),
        (b'S', b"GC"),
        (b'W', b"AT"),
        (b'K', b"GT"),
        (b'M', b"AC"),
    ],
    &[(b'B', b"CGT"), (b'D', b"AGT"), (b'H', b"ACT"), (b'V', b"ACG")],
    &[(b'N', b"ACGT")],
];

/// Does an observed literal base (`A`/`C`/`G`/`T`, or occasionally a
/// sequencer no-call `N`) satisfy an IUPAC code from a repeat motif? A
/// plain-letter code only matches itself; an ambiguity code matches any
/// base in its represented set (see [`IUPAC_TIERS`]).
fn iupac_matches(observed: u8, code: u8) -> bool {
    if observed == code {
        return true;
    }
    IUPAC_TIERS
        .iter()
        .flat_map(|tier| tier.iter())
        .find(|&&(tier_code, _)| tier_code == code)
        .is_some_and(|&(_, members)| members.contains(&observed))
}

/// Default shortest repeat-unit (motif) length to consider.
pub const DEFAULT_MOTIF_MIN_LEN: u32 = 2;
/// Default longest repeat-unit (motif) length to consider.
pub const DEFAULT_MOTIF_MAX_LEN: u32 = 20;
/// Default degenerate-position cap for a mononucleotide (1bp) motif.
pub const DEFAULT_MAX_DEGENERATE_MONONUCLEOTIDE: u32 = 0;
/// Default degenerate-position cap for a dinucleotide (2bp) motif.
pub const DEFAULT_MAX_DEGENERATE_DINUCLEOTIDE: u32 = 0;
/// Default degenerate-position cap for a trinucleotide (3bp) motif.
pub const DEFAULT_MAX_DEGENERATE_TRINUCLEOTIDE: u32 = 1;
/// Default degenerate-position cap for any motif of 4bp or longer.
pub const DEFAULT_MAX_DEGENERATE_OTHER: u32 = 2;

/// The number of `unit`'s positions that are an IUPAC ambiguity code
/// rather than a plain A/C/G/T base.
fn degenerate_count(unit: &[u8]) -> u32 {
    unit.iter().filter(|&&b| !matches!(b, b'A' | b'C' | b'G' | b'T')).count() as u32
}

/// Does `unit` have more degenerate positions than allowed for its length?
/// Mononucleotide (1bp), dinucleotide (2bp), and trinucleotide (3bp)
/// motifs are short enough that even one or two degenerate positions
/// dominate the whole motif, so each gets its own, stricter cap; every
/// motif of 4bp or longer shares `max_other`.
fn exceeds_degenerate_limit(
    unit: &[u8],
    max_mononucleotide: u32,
    max_dinucleotide: u32,
    max_trinucleotide: u32,
    max_other: u32,
) -> bool {
    let limit = match unit.len() {
        1 => max_mononucleotide,
        2 => max_dinucleotide,
        3 => max_trinucleotide,
        _ => max_other,
    };
    degenerate_count(unit) > limit
}

/// Returns the canonical repeat unit if `bases` (paired with `quals`) look
/// like an in-repeat read, `None` otherwise. The unit may contain IUPAC
/// ambiguity codes at consistently-mixed positions (see the module docs),
/// but is rejected if it has more degenerate positions than
/// `max_degenerate_mononucleotide`/`max_degenerate_dinucleotide`/
/// `max_degenerate_trinucleotide`/`max_degenerate_other` allow for its
/// length (see [`exceeds_degenerate_limit`]).
#[allow(clippy::too_many_arguments)]
pub fn classify_in_repeat_read(
    bases: &[u8],
    quals: &[u8],
    motif_min_len: u32,
    motif_max_len: u32,
    max_degenerate_mononucleotide: u32,
    max_degenerate_dinucleotide: u32,
    max_degenerate_trinucleotide: u32,
    max_degenerate_other: u32,
) -> Option<Vec<u8>> {
    let unit =
        compute_canonical_repeat_unit_with_frequency(MIN_UNIT_FREQUENCY, bases, quals, motif_min_len, motif_max_len)?;
    if unit.is_empty()
        || exceeds_degenerate_limit(
            &unit,
            max_degenerate_mononucleotide,
            max_degenerate_dinucleotide,
            max_degenerate_trinucleotide,
            max_degenerate_other,
        )
    {
        return None;
    }

    let units_shifts = shift_units(std::slice::from_ref(&unit));
    let score = match_repeat_rc(&units_shifts, bases, quals) / bases.len() as f64;

    if score >= MIN_IRR_SCORE {
        Some(unit)
    } else {
        None
    }
}

/// Returns every canonical repeat unit (motif) that independently clears
/// both the unit-frequency and IRR-score thresholds for this read. A read
/// can plausibly satisfy more than one motif (e.g. a longer compound period
/// that is itself a repetition of a shorter one also scores well), so
/// unlike [`classify_in_repeat_read`] this doesn't collapse to a single
/// "best" motif. Each returned motif is distinct; order is not significant.
/// A motif may contain IUPAC ambiguity codes at consistently-mixed
/// positions (see the module docs), but is excluded if it has more
/// degenerate positions than `max_degenerate_mononucleotide`/
/// `max_degenerate_dinucleotide`/`max_degenerate_trinucleotide`/
/// `max_degenerate_other` allow for its length (see
/// [`exceeds_degenerate_limit`]).
#[allow(clippy::too_many_arguments)]
pub fn classify_in_repeat_read_all(
    bases: &[u8],
    quals: &[u8],
    motif_min_len: u32,
    motif_max_len: u32,
    max_degenerate_mononucleotide: u32,
    max_degenerate_dinucleotide: u32,
    max_degenerate_trinucleotide: u32,
    max_degenerate_other: u32,
) -> Vec<Vec<u8>> {
    let smallest_period = motif_min_len.max(1) as usize;
    let largest_period = (motif_max_len as usize).min(bases.len() / 2 + 1);

    let mut motifs: Vec<Vec<u8>> = Vec::new();
    for period in smallest_period..=largest_period {
        if match_frequency_at_offset(period, bases) < MIN_UNIT_FREQUENCY {
            continue;
        }

        let mut unit = extract_consensus_repeat_unit_iupac(period, bases, quals);

        const PERFECT_MATCH_FREQUENCY: f64 = 1.0;
        let (reduction_min, reduction_max) = DEFAULT_REDUCTION_MOTIF_RANGE;
        if let Some(reduced_period) =
            smallest_frequent_period(PERFECT_MATCH_FREQUENCY, &unit, reduction_min, reduction_max)
            && reduced_period != period
        {
            unit = extract_consensus_repeat_unit(reduced_period, &unit);
        }

        let canonical = compute_canonical_repeat_unit(&unit);
        if canonical.is_empty()
            || exceeds_degenerate_limit(
                &canonical,
                max_degenerate_mononucleotide,
                max_degenerate_dinucleotide,
                max_degenerate_trinucleotide,
                max_degenerate_other,
            )
            || motifs.contains(&canonical)
        {
            continue;
        }

        let units_shifts = shift_units(std::slice::from_ref(&canonical));
        let score = match_repeat_rc(&units_shifts, bases, quals) / bases.len() as f64;
        if score >= MIN_IRR_SCORE {
            motifs.push(canonical);
        }
    }

    motifs
}

// --- motif detection ---------------------------------------------------

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
    let num_matches = bases[..max_matches].iter().zip(&bases[offset..]).filter(|(a, b)| a == b).count();
    num_matches as f64 / max_matches as f64
}

/// Finds the shortest motif period whose match frequency is at least as
/// good as any longer period's, or `None` if none clears `min_frequency`.
fn smallest_frequent_period(min_frequency: f64, bases: &[u8], motif_min_len: u32, motif_max_len: u32) -> Option<usize> {
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

fn extract_consensus_base(offset: usize, period: usize, bases: &[u8]) -> u8 {
    let mut counts = [0u32; 256];
    let mut index = offset;
    while index < bases.len() {
        counts[bases[index] as usize] += 1;
        index += period;
    }

    // Ties (equal counts) are broken deterministically by preferring the
    // larger byte value.
    counts
        .iter()
        .enumerate()
        .filter(|&(_, &count)| count > 0)
        .max_by_key(|&(base, &count)| (count, base))
        .map(|(base, _)| base as u8)
        .unwrap_or(b'?')
}

fn extract_consensus_repeat_unit(period: usize, bases: &[u8]) -> Vec<u8> {
    (0..period).map(|offset| extract_consensus_base(offset, period, bases)).collect()
}

/// Like [`extract_consensus_base`], but may call an IUPAC ambiguity code
/// (see [`IUPAC_TIERS`]) instead of a single base, when a phase position is
/// a genuine mix across repeat copies rather than one dominant base.
///
/// Two things guard against false ambiguity calls:
/// - Only *confidently-called* bases (`quals[i] >= MIN_BASE_QUALITY`) count
///   toward a tier's coverage, so an otherwise-clean position doesn't get
///   marked ambiguous just because its noisy/low-quality observations
///   happen to disagree -- that kind of noise is what the quality-aware
///   scoring in [`match_units`] already exists to tolerate, not something
///   the *motif* itself should absorb.
/// - `period == 1` never escalates: a length-1 repeat unit is a homopolymer
///   by definition, and an "ambiguous homopolymer" covering 2+ bases is
///   nearly a wildcard -- it would happily paper over two unrelated
///   same-length runs (e.g. `AAAAACCCCC`) as one fake "IRR" instead of
///   correctly rejecting them as non-repetitive.
///
/// Below [`MIN_DEGENERACY_SAMPLES`] confident observations, or when no
/// tier below `N` clears [`MIN_BASE_PURITY`], falls back to
/// [`extract_consensus_base`] (i.e. today's plain majority vote, quality
///-blind, over every observed base).
fn extract_consensus_base_iupac(offset: usize, period: usize, bases: &[u8], quals: &[u8]) -> u8 {
    if period >= 2 {
        let mut hq_counts = [0u32; 256];
        let mut index = offset;
        while index < bases.len() {
            if quals[index] >= MIN_BASE_QUALITY {
                hq_counts[bases[index] as usize] += 1;
            }
            index += period;
        }

        let acgt_total: u32 = [b'A', b'C', b'G', b'T'].iter().map(|&b| hq_counts[b as usize]).sum();
        if acgt_total >= MIN_DEGENERACY_SAMPLES {
            for tier in IUPAC_TIERS {
                let (code, covered) = tier
                    .iter()
                    .map(|&(code, members)| (code, members.iter().map(|&b| hq_counts[b as usize]).sum::<u32>()))
                    .max_by_key(|&(code, covered)| (covered, code))
                    .expect("tiers are non-empty");
                if covered as f64 / acgt_total as f64 >= MIN_BASE_PURITY {
                    return code;
                }
            }
        }
    }

    extract_consensus_base(offset, period, bases)
}

/// Like [`extract_consensus_repeat_unit`], but positions may come back as
/// an IUPAC ambiguity code -- see [`extract_consensus_base_iupac`].
fn extract_consensus_repeat_unit_iupac(period: usize, bases: &[u8], quals: &[u8]) -> Vec<u8> {
    (0..period).map(|offset| extract_consensus_base_iupac(offset, period, bases, quals)).collect()
}

fn minimal_unit_under_shift(unit: &[u8]) -> Vec<u8> {
    let len = unit.len();
    let mut doubled = unit.to_vec();
    doubled.extend_from_slice(unit);
    let best_offset = (0..len).min_by_key(|&offset| &doubled[offset..offset + len]).unwrap_or(0);
    doubled[best_offset..best_offset + len].to_vec()
}

/// The IUPAC complement of a single base or ambiguity code: `A<->T`,
/// `C<->G`, `R<->Y`, `K<->M`, `B<->V`, `D<->H`, and `S`/`W`/`N` (each
/// self-complementary, since complementing every base in their represented
/// set yields the same set back).
fn complement_base(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        b'S' | b'W' | b'N' => base,
        _ => b'N',
    }
}

fn reverse_complement(bases: &[u8]) -> Vec<u8> {
    bases.iter().rev().map(|&base| complement_base(base)).collect()
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

fn compute_canonical_repeat_unit_with_frequency(
    min_frequency: f64,
    bases: &[u8],
    quals: &[u8],
    motif_min_len: u32,
    motif_max_len: u32,
) -> Option<Vec<u8>> {
    let period = smallest_frequent_period(min_frequency, bases, motif_min_len, motif_max_len)?;
    let mut motif = extract_consensus_repeat_unit_iupac(period, bases, quals);

    const PERFECT_MATCH_FREQUENCY: f64 = 1.0;
    let (reduction_min, reduction_max) = DEFAULT_REDUCTION_MOTIF_RANGE;
    if let Some(reduced_period) = smallest_frequent_period(PERFECT_MATCH_FREQUENCY, &motif, reduction_min, reduction_max)
        && reduced_period != period
    {
        motif = extract_consensus_repeat_unit(reduced_period, &motif);
    }

    Some(compute_canonical_repeat_unit(&motif))
}

// --- quality-aware matching ----------------------------------------------

fn shift_units(units: &[Vec<u8>]) -> Vec<Vec<Vec<u8>>> {
    let unit_length = units[0].len();
    let extended: Vec<Vec<u8>> = units
        .iter()
        .map(|unit| {
            let mut doubled = unit.clone();
            doubled.extend_from_slice(unit);
            doubled
        })
        .collect();

    (0..unit_length)
        .map(|offset| extended.iter().map(|doubled| doubled[offset..offset + unit_length].to_vec()).collect())
        .collect()
}

fn match_units(units: &[Vec<u8>], bases: &[u8], quals: &[u8], min_baseq: u8) -> f64 {
    const MATCH_SCORE: f64 = 1.0;
    const LOWQUAL_MISMATCH_SCORE: f64 = 0.5;
    const MISMATCH_PENALTY: f64 = -1.0;

    units
        .iter()
        .map(|unit| {
            bases
                .iter()
                .zip(quals)
                .zip(unit)
                .map(|((&base, &qual), &unit_base)| {
                    if iupac_matches(base, unit_base) {
                        MATCH_SCORE
                    } else if qual < min_baseq {
                        LOWQUAL_MISMATCH_SCORE
                    } else {
                        MISMATCH_PENALTY
                    }
                })
                .sum::<f64>()
        })
        .fold(f64::NEG_INFINITY, f64::max)
}

fn match_repeat(units: &[Vec<u8>], bases: &[u8], quals: &[u8], min_baseq: u8) -> f64 {
    let unit_len = units[0].len();
    bases
        .chunks(unit_len)
        .zip(quals.chunks(unit_len))
        .map(|(base_chunk, qual_chunk)| match_units(units, base_chunk, qual_chunk, min_baseq))
        .sum()
}

fn match_repeat_with_shifts(units_shifts: &[Vec<Vec<u8>>], bases: &[u8], quals: &[u8], min_baseq: u8) -> f64 {
    units_shifts
        .iter()
        .map(|units_shift| match_repeat(units_shift, bases, quals, min_baseq))
        .fold(f64::NEG_INFINITY, f64::max)
}

fn match_repeat_rc(units_shifts: &[Vec<Vec<u8>>], bases: &[u8], quals: &[u8]) -> f64 {
    let forward_score = match_repeat_with_shifts(units_shifts, bases, quals, MIN_BASE_QUALITY);

    let bases_rc = reverse_complement(bases);
    let mut quals_rc = quals.to_vec();
    quals_rc.reverse();
    let reverse_score = match_repeat_with_shifts(units_shifts, &bases_rc, &quals_rc, MIN_BASE_QUALITY);

    forward_score.max(reverse_score)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Test qualities are given as ASCII PHRED+33 (as if lifted straight
    /// from a FASTQ/SAM text field); convert to the raw PHRED bytes our
    /// functions expect.
    fn raw_quals(ascii: &str) -> Vec<u8> {
        ascii.bytes().map(|c| c - 33).collect()
    }

    /// `classify_in_repeat_read` at the default degenerate-position caps,
    /// for tests that don't care about that parameter specifically.
    fn classify_default(bases: &[u8], quals: &[u8], motif_min_len: u32, motif_max_len: u32) -> Option<Vec<u8>> {
        classify_in_repeat_read(
            bases,
            quals,
            motif_min_len,
            motif_max_len,
            DEFAULT_MAX_DEGENERATE_MONONUCLEOTIDE,
            DEFAULT_MAX_DEGENERATE_DINUCLEOTIDE,
            DEFAULT_MAX_DEGENERATE_TRINUCLEOTIDE,
            DEFAULT_MAX_DEGENERATE_OTHER,
        )
    }

    /// `classify_in_repeat_read_all` at the default degenerate-position
    /// caps, for tests that don't care about that parameter specifically.
    fn classify_all_default(bases: &[u8], quals: &[u8], motif_min_len: u32, motif_max_len: u32) -> Vec<Vec<u8>> {
        classify_in_repeat_read_all(
            bases,
            quals,
            motif_min_len,
            motif_max_len,
            DEFAULT_MAX_DEGENERATE_MONONUCLEOTIDE,
            DEFAULT_MAX_DEGENERATE_DINUCLEOTIDE,
            DEFAULT_MAX_DEGENERATE_TRINUCLEOTIDE,
            DEFAULT_MAX_DEGENERATE_OTHER,
        )
    }

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
        assert_eq!(smallest_frequent_period(0.85, b"GGCCCCGGCCCC", 1, 20), Some(6));
        assert_eq!(smallest_frequent_period(0.85, b"ATGATCATGATGATGATGATG", 1, 20), Some(6));
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

    /// High, uniform per-base quality, for tests that don't care about
    /// quality-gated degeneracy and just want the plain periodicity/
    /// canonicalization behavior exercised.
    fn hq(len: usize) -> Vec<u8> {
        vec![40u8; len]
    }

    #[test]
    fn compute_canonical_repeat_unit_with_frequency_typical() {
        let bases = b"CGGCGCCGGCGG";
        assert_eq!(
            compute_canonical_repeat_unit_with_frequency(0.8, bases, &hq(bases.len()), 1, 20),
            Some(b"CCG".to_vec())
        );
        assert_eq!(
            compute_canonical_repeat_unit_with_frequency(0.85, bases, &hq(bases.len()), 1, 20),
            None
        );
        let bases = b"ACCCCAACCCCAACCCCAACCCCAACCCCAACCCCA";
        assert_eq!(
            compute_canonical_repeat_unit_with_frequency(0.8, bases, &hq(bases.len()), 1, 20),
            Some(b"AACCCC".to_vec())
        );
    }

    #[test]
    fn compute_canonical_repeat_unit_with_frequency_homopolymer() {
        let bases = b"CCCCCCC";
        assert_eq!(
            compute_canonical_repeat_unit_with_frequency(1.0, bases, &hq(bases.len()), 1, 20),
            Some(b"C".to_vec())
        );
    }

    #[test]
    fn classify_in_repeat_read_typical_cases() {
        assert_eq!(
            classify_default(b"CCCCC", &raw_quals("$$$$$"), 1, 20),
            Some(b"C".to_vec())
        );

        assert_eq!(classify_default(b"AAAAACCCCC", &raw_quals("$$$$$$$$$$"), 1, 20), None);

        let bases: &[u8] = concat!(
            "TCCACCCACCTCACCCCCCCCCCCCCCCGCCCCCCCCCCACCCCCCCCGCCCCCCCCCCCGGCCCCCCACTCCCCCCCCCCGGTCCTCCCC",
            "CCCCCCCACCCTCCCCCCCCGCCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCCC"
        )
        .as_bytes();
        let quals = raw_quals(concat!(
            "------7----7-----7-777-7-F<--777F777F<J-7--7-7-A7-AFJA<<A-<<-7--7A77---7A-77A77A7---7-7-",
            "7--77-7-77-777---7<7A<A-7A)-)-<)7))77A<JJF))--A<F-)-<-)<---7<J"
        ));
        assert_eq!(classify_default(bases, &quals, 1, 70), None);

        let bases: &[u8] = concat!(
            "TCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTT",
            "CATTTCATTTCATTTCATTTCATTTCTTTTTTTTTATTTTTTTTTATTTTATATCGGAT"
        )
        .as_bytes();
        let quals = raw_quals(concat!(
            "((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((",
            "(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((("
        ));
        assert_eq!(classify_default(bases, &quals, 1, 20), Some(b"AAATG".to_vec()));

        let bases: &[u8] = concat!(
            "CCCGCGCCCCGCCCCGCGCCCCGCCCCGCGCCCCGCCCCGCGCCCCGCCCCGCGCCCCGCCCCGCGCCCCGCCCCGCGCCCCGCCCCCCGCCCCGCC",
            "CCGCGCCCCGCCCCGCGCCCCGCCCCGCGCCCCGCCCCGCGCCCCGCCCCGCG"
        )
        .as_bytes();
        let quals = raw_quals(concat!(
            "((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((",
            "(((((((((((((((((((((((((((((((((((((((((((((((((((((((("
        ));
        assert_eq!(classify_default(bases, &quals, 1, 15), Some(b"CCCCGCCCCGCG".to_vec()));

        let bases: &[u8] = concat!(
            "GGGGCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCG",
            "GGGCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCGGGGCG"
        )
        .as_bytes();
        let quals = raw_quals(concat!(
            "((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((",
            "(((((((((((((((((((((((((((((((((((((((((((((((((((((((("
        ));
        assert_eq!(classify_default(bases, &quals, 1, 20), Some(b"CCCCGCCCCGCG".to_vec()));
    }

    #[test]
    fn classify_in_repeat_read_rejects_n_bases() {
        assert_eq!(classify_default(b"NNNNN", &raw_quals("$$$$$"), 1, 20), None);
    }

    /// A read whose composition is 20 A's followed by a single G, repeated:
    /// mostly-A with only rare G interruptions passes the homopolymer "A"
    /// motif's thresholds despite not being a pure homopolymer, while the
    /// exact 21bp repeat unit also passes on its own -- so this read
    /// legitimately qualifies under two distinct, non-harmonic motifs
    /// (unlike a pure homopolymer, where every candidate period reduces
    /// back to the same single-base canonical unit).
    fn mostly_a_with_rare_g() -> Vec<u8> {
        let unit: Vec<u8> = (0..20).map(|_| b'A').chain(std::iter::once(b'G')).collect();
        unit.iter().cloned().cycle().take(21 * 16).collect()
    }

    #[test]
    fn classify_in_repeat_read_all_returns_single_motif_for_pure_repeat() {
        let bases = "CAG".repeat(20).into_bytes();
        let quals = vec![40u8; bases.len()];
        assert_eq!(classify_all_default(&bases, &quals, 1, 20), vec![b"AGC".to_vec()]);
    }

    #[test]
    fn classify_in_repeat_read_all_returns_empty_for_non_repetitive() {
        let bases = b"ACGTTGCAACGGTTCAGTAGCTAGCATCGATCGTAGCTAGGCTAGCATCGTAGCTAGCA";
        let quals = vec![40u8; bases.len()];
        assert!(classify_all_default(bases, &quals, 1, 20).is_empty());
    }

    #[test]
    fn classify_in_repeat_read_all_returns_multiple_distinct_motifs() {
        let bases = mostly_a_with_rare_g();
        let quals = vec![40u8; bases.len()];

        let motifs = classify_all_default(&bases, &quals, 1, 30);

        assert!(
            motifs.contains(&b"A".to_vec()),
            "expected the mostly-A homopolymer motif to qualify: {motifs:?}"
        );
        assert!(
            motifs.iter().any(|m| m.len() == 21),
            "expected the exact 21bp repeat unit to also qualify: {motifs:?}"
        );
        // Not asserting an exact count: this fixture's rare-but-real G
        // interruptions are confident (uniform high quality), so several
        // *other* periods can also legitimately turn up a partially
        // degenerate (R-containing) motif that independently clears the
        // score threshold. That's expected now that degenerate calling
        // exists -- this test only cares that the two motifs it names
        // above are among whatever comes back.

        // The single-motif classifier only ever returns one of them.
        let single = classify_default(&bases, &quals, 1, 30);
        assert!(single.is_some());
        assert!(motifs.contains(single.as_ref().unwrap()));
    }

    // --- IUPAC degenerate-motif calling ---------------------------------

    #[test]
    fn iupac_matches_basic() {
        assert!(iupac_matches(b'A', b'A'));
        assert!(!iupac_matches(b'A', b'C'));

        assert!(iupac_matches(b'A', b'R'));
        assert!(iupac_matches(b'G', b'R'));
        assert!(!iupac_matches(b'C', b'R'));
        assert!(!iupac_matches(b'T', b'R'));

        for base in [b'A', b'C', b'G', b'T'] {
            assert!(iupac_matches(base, b'N'), "N should match {}", base as char);
        }
    }

    #[test]
    fn complement_base_iupac_pairs() {
        assert_eq!(complement_base(b'A'), b'T');
        assert_eq!(complement_base(b'T'), b'A');
        assert_eq!(complement_base(b'C'), b'G');
        assert_eq!(complement_base(b'G'), b'C');
        assert_eq!(complement_base(b'R'), b'Y');
        assert_eq!(complement_base(b'Y'), b'R');
        assert_eq!(complement_base(b'K'), b'M');
        assert_eq!(complement_base(b'M'), b'K');
        assert_eq!(complement_base(b'B'), b'V');
        assert_eq!(complement_base(b'V'), b'B');
        assert_eq!(complement_base(b'D'), b'H');
        assert_eq!(complement_base(b'H'), b'D');
        // Self-complementary: complementing every base in the represented
        // set yields the same set back.
        assert_eq!(complement_base(b'S'), b'S');
        assert_eq!(complement_base(b'W'), b'W');
        assert_eq!(complement_base(b'N'), b'N');
    }

    #[test]
    fn reverse_complement_handles_iupac_codes() {
        assert_eq!(reverse_complement(b"GCN"), b"NGC");
        assert_eq!(reverse_complement(b"AARRG"), b"CYYTT");
    }

    #[test]
    fn extract_consensus_base_iupac_calls_ambiguity_code_for_genuine_mixture() {
        // A period-2 unit where phase 0 evenly alternates A/G copy-to-copy
        // (10 copies, well above the sample-size gate) while phase 1 stays
        // pure C: phase 0 should resolve to R (A/G), not force one or the
        // other.
        let bases: Vec<u8> =
            (0..10).flat_map(|i| [if i % 2 == 0 { b'A' } else { b'G' }, b'C']).collect();
        let quals = vec![40u8; bases.len()];

        assert_eq!(extract_consensus_base_iupac(0, 2, &bases, &quals), b'R');
        assert_eq!(extract_consensus_base_iupac(1, 2, &bases, &quals), b'C');
    }

    #[test]
    fn extract_consensus_base_iupac_ignores_low_quality_votes() {
        // Same 50/50 A/G mixture at phase 0, but every vote is low-quality:
        // not enough *confident* observations to call ambiguity, so this
        // falls back to the plain majority vote (tie broken toward the
        // larger byte value, i.e. G), not R.
        let bases: Vec<u8> =
            (0..10).flat_map(|i| [if i % 2 == 0 { b'A' } else { b'G' }, b'C']).collect();
        let quals = vec![5u8; bases.len()];

        assert_eq!(extract_consensus_base_iupac(0, 2, &bases, &quals), b'G');
    }

    #[test]
    fn extract_consensus_base_iupac_never_escalates_at_period_one() {
        // A clean, confident, well-sampled 50/50 A/G split -- but at
        // period 1 (a homopolymer position) ambiguity calling must never
        // trigger, regardless of sample size or quality.
        let bases: Vec<u8> = (0..10).flat_map(|_| [b'A', b'G']).collect();
        let quals = vec![40u8; bases.len()];

        assert_eq!(extract_consensus_base_iupac(0, 1, &bases, &quals), b'G');
    }

    /// `extract_consensus_repeat_unit_iupac` is tested directly at a fixed
    /// period below (like the plain `extract_consensus_repeat_unit_basic`
    /// test above it), rather than through the full
    /// `classify_in_repeat_read_all` pipeline: a position that's genuinely
    /// unconstrained (as in a literal `GCN`) or split not-quite-evenly
    /// (as in `AARRG`) measurably dilutes the *raw*, IUPAC-unaware
    /// literal-byte periodicity signal that the outer period search relies
    /// on -- correctly so, since from that check's point of view alone,
    /// weak-to-nonexistent periodicity at one whole position out of a
    /// short motif is genuinely weak evidence of any period at all. Real
    /// reads carrying a true degenerate position are usually much longer
    /// relative to one ambiguous slot than these minimal fixtures, so the
    /// aggregate signal clears the period-detection bar fine in practice
    /// (see `classify_in_repeat_read_all_returns_multiple_distinct_motifs`
    /// above for an end-to-end example). These two tests instead isolate
    /// exactly the new piece of logic: given a period, does the consensus
    /// step correctly call each position?

    #[test]
    fn extract_consensus_repeat_unit_iupac_calls_gcn_style_motif() {
        // A `GC` repeat whose 3rd position cycles evenly through all 4
        // bases across repeat copies: no single base, pair, or triple
        // covers enough of it, so it must resolve to the fully degenerate
        // `N` code.
        let third = [b'A', b'C', b'G', b'T'];
        let bases: Vec<u8> = (0..40).flat_map(|i| [b'G', b'C', third[i % 4]]).collect();
        let quals = vec![40u8; bases.len()];

        assert_eq!(extract_consensus_repeat_unit_iupac(3, &bases, &quals), b"GCN");
    }

    #[test]
    fn extract_consensus_repeat_unit_iupac_calls_aarrg_style_motif() {
        // Positions 0, 1, 4 are always A, A, G; positions 2, 3 are an even
        // (10/10), irregularly-ordered mix of A and G across repeat
        // copies -- each should resolve to the purine ambiguity code `R`.
        // (Irregular, not a clean alternation: alternating every copy is
        // itself an exact period twice as long, which is a different,
        // non-degenerate case already covered by the `ATAT`-collapses-to-
        // `AT` reasoning elsewhere in this module.)
        let purines: [u8; 20] = [
            b'A', b'G', b'G', b'A', b'A', b'G', b'A', b'G', b'G', b'A', b'G', b'A', b'A', b'G',
            b'A', b'G', b'G', b'A', b'G', b'A',
        ];
        let bases: Vec<u8> = purines.iter().flat_map(|&p| [b'A', b'A', p, p, b'G']).collect();
        let quals = vec![40u8; bases.len()];

        assert_eq!(extract_consensus_repeat_unit_iupac(5, &bases, &quals), b"AARRG");
    }

    #[test]
    fn degenerate_count_basic() {
        assert_eq!(degenerate_count(b"AAATG"), 0);
        assert_eq!(degenerate_count(b"GCN"), 1);
        assert_eq!(degenerate_count(b"AARRG"), 2);
        assert_eq!(degenerate_count(b"NNNN"), 4);
    }

    #[test]
    fn exceeds_degenerate_limit_uses_length_specific_defaults() {
        let defaults = |unit: &[u8]| {
            exceeds_degenerate_limit(
                unit,
                DEFAULT_MAX_DEGENERATE_MONONUCLEOTIDE,
                DEFAULT_MAX_DEGENERATE_DINUCLEOTIDE,
                DEFAULT_MAX_DEGENERATE_TRINUCLEOTIDE,
                DEFAULT_MAX_DEGENERATE_OTHER,
            )
        };

        // Mononucleotide (1bp): default cap 0.
        assert!(!defaults(b"A"), "a clean mononucleotide motif should never be rejected");
        assert!(defaults(b"N"), "any degenerate position should exceed the mononucleotide default of 0");

        // Dinucleotide (2bp): default cap 0.
        assert!(!defaults(b"AC"), "a clean dinucleotide motif should never be rejected");
        assert!(defaults(b"AR"), "any degenerate position should exceed the dinucleotide default of 0");

        // Trinucleotide (3bp): default cap 1.
        assert!(!defaults(b"GCN"), "one degenerate position should be within the trinucleotide default of 1");
        assert!(defaults(b"RCN"), "two degenerate positions should exceed the trinucleotide default of 1");

        // Everything else (4bp+): default cap 2.
        assert!(!defaults(b"AARRG"), "two degenerate positions should be within the \"other\" default of 2");
        assert!(defaults(b"ARRRG"), "three degenerate positions should exceed the \"other\" default of 2");
    }

    /// An 8bp unit: 5 fixed A's, then 3 positions that mix evenly (10/10)
    /// between A and G, cycling copy-to-copy through an `AAGG` pattern
    /// (period 4 in copy-index, i.e. not aligned with the target period-8
    /// grouping) -- keeping the raw literal periodicity signal strong
    /// enough to be found (5/8 positions always match exactly, and the
    /// 3 mixed positions still match at the adjacent-copy lag about half
    /// the time) while pushing the resulting unit's degenerate count (3)
    /// just over the default "other" cap of 2.
    fn mostly_pure_with_three_mixed_positions() -> Vec<u8> {
        let purine_cycle = [b'A', b'A', b'G', b'G'];
        (0..20usize)
            .flat_map(|i| {
                let p = purine_cycle[i % 4];
                [b'A', b'A', b'A', b'A', b'A', p, p, p]
            })
            .collect()
    }

    #[test]
    fn classify_in_repeat_read_all_rejects_motif_over_default_degenerate_limit() {
        let bases = mostly_pure_with_three_mixed_positions();
        let quals = vec![40u8; bases.len()];

        let motifs = classify_all_default(&bases, &quals, 2, 10);
        assert!(
            !motifs.iter().any(|m| m.len() == 8),
            "an 8bp motif with 3 degenerate positions should be rejected at the default \
             \"other\" cap of 2: {motifs:?}"
        );
    }

    #[test]
    fn classify_in_repeat_read_all_allows_motif_under_relaxed_degenerate_limit() {
        let bases = mostly_pure_with_three_mixed_positions();
        let quals = vec![40u8; bases.len()];

        let motifs = classify_in_repeat_read_all(
            &bases,
            &quals,
            2,
            10,
            DEFAULT_MAX_DEGENERATE_MONONUCLEOTIDE,
            DEFAULT_MAX_DEGENERATE_DINUCLEOTIDE,
            DEFAULT_MAX_DEGENERATE_TRINUCLEOTIDE,
            3,
        );
        let eight_mer = motifs.iter().find(|m| m.len() == 8);
        assert!(
            eight_mer.is_some_and(|m| degenerate_count(m) == 3),
            "expected the 3-degenerate 8bp motif to survive with a relaxed \"other\" cap of 3: {motifs:?}"
        );
    }
}
