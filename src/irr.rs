//! In-repeat-read (IRR) detection heuristic: decides whether a read's
//! sequence is dominated by repetitions of some short motif, and if so,
//! returns the motif's canonical repeat unit.
//!
//! `bases` are expected to be uppercase decoded read sequence bytes (as
//! returned by `rust_htslib::bam::record::Seq::as_bytes`), and `quals` are
//! raw (non-ASCII-offset) PHRED scores (as returned by `Record::qual`).

const MIN_UNIT_FREQUENCY: f64 = 0.8;
const MIN_IRR_SCORE: f64 = 0.90;
const MIN_BASE_QUALITY: u8 = 20;
const DEFAULT_REDUCTION_MOTIF_RANGE: (u32, u32) = (1, 20);

/// Default shortest repeat-unit (motif) length to consider.
pub const DEFAULT_MOTIF_MIN_LEN: u32 = 2;
/// Default longest repeat-unit (motif) length to consider.
pub const DEFAULT_MOTIF_MAX_LEN: u32 = 20;

/// Returns the canonical repeat unit if `bases` (paired with `quals`) look
/// like an in-repeat read, `None` otherwise.
pub fn classify_in_repeat_read(bases: &[u8], quals: &[u8], motif_min_len: u32, motif_max_len: u32) -> Option<Vec<u8>> {
    let unit =
        compute_canonical_repeat_unit_with_frequency(MIN_UNIT_FREQUENCY, bases, motif_min_len, motif_max_len)?;
    if unit.is_empty() || unit == b"N" {
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

fn minimal_unit_under_shift(unit: &[u8]) -> Vec<u8> {
    let len = unit.len();
    let mut doubled = unit.to_vec();
    doubled.extend_from_slice(unit);
    let best_offset = (0..len).min_by_key(|&offset| &doubled[offset..offset + len]).unwrap_or(0);
    doubled[best_offset..best_offset + len].to_vec()
}

fn reverse_complement(bases: &[u8]) -> Vec<u8> {
    bases
        .iter()
        .rev()
        .map(|&base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
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

fn compute_canonical_repeat_unit_with_frequency(
    min_frequency: f64,
    bases: &[u8],
    motif_min_len: u32,
    motif_max_len: u32,
) -> Option<Vec<u8>> {
    let period = smallest_frequent_period(min_frequency, bases, motif_min_len, motif_max_len)?;
    let mut motif = extract_consensus_repeat_unit(period, bases);

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
                    if base == unit_base {
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

    #[test]
    fn compute_canonical_repeat_unit_with_frequency_typical() {
        assert_eq!(
            compute_canonical_repeat_unit_with_frequency(0.8, b"CGGCGCCGGCGG", 1, 20),
            Some(b"CCG".to_vec())
        );
        assert_eq!(compute_canonical_repeat_unit_with_frequency(0.85, b"CGGCGCCGGCGG", 1, 20), None);
        assert_eq!(
            compute_canonical_repeat_unit_with_frequency(0.8, b"ACCCCAACCCCAACCCCAACCCCAACCCCAACCCCA", 1, 20),
            Some(b"AACCCC".to_vec())
        );
    }

    #[test]
    fn compute_canonical_repeat_unit_with_frequency_homopolymer() {
        assert_eq!(
            compute_canonical_repeat_unit_with_frequency(1.0, b"CCCCCCC", 1, 20),
            Some(b"C".to_vec())
        );
    }

    #[test]
    fn classify_in_repeat_read_typical_cases() {
        assert_eq!(
            classify_in_repeat_read(b"CCCCC", &raw_quals("$$$$$"), 1, 20),
            Some(b"C".to_vec())
        );

        assert_eq!(
            classify_in_repeat_read(b"AAAAACCCCC", &raw_quals("$$$$$$$$$$"), 1, 20),
            None
        );

        let bases: &[u8] = concat!(
            "TCCACCCACCTCACCCCCCCCCCCCCCCGCCCCCCCCCCACCCCCCCCGCCCCCCCCCCCGGCCCCCCACTCCCCCCCCCCGGTCCTCCCC",
            "CCCCCCCACCCTCCCCCCCCGCCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCCC"
        )
        .as_bytes();
        let quals = raw_quals(concat!(
            "------7----7-----7-777-7-F<--777F777F<J-7--7-7-A7-AFJA<<A-<<-7--7A77---7A-77A77A7---7-7-",
            "7--77-7-77-777---7<7A<A-7A)-)-<)7))77A<JJF))--A<F-)-<-)<---7<J"
        ));
        assert_eq!(classify_in_repeat_read(bases, &quals, 1, 70), None);

        let bases: &[u8] = concat!(
            "TCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTT",
            "CATTTCATTTCATTTCATTTCATTTCTTTTTTTTTATTTTTTTTTATTTTATATCGGAT"
        )
        .as_bytes();
        let quals = raw_quals(concat!(
            "((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((",
            "(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((("
        ));
        assert_eq!(classify_in_repeat_read(bases, &quals, 1, 20), Some(b"AAATG".to_vec()));

        let bases: &[u8] = concat!(
            "CCCGCGCCCCGCCCCGCGCCCCGCCCCGCGCCCCGCCCCGCGCCCCGCCCCGCGCCCCGCCCCGCGCCCCGCCCCGCGCCCCGCCCCCCGCCCCGCC",
            "CCGCGCCCCGCCCCGCGCCCCGCCCCGCGCCCCGCCCCGCGCCCCGCCCCGCG"
        )
        .as_bytes();
        let quals = raw_quals(concat!(
            "((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((",
            "(((((((((((((((((((((((((((((((((((((((((((((((((((((((("
        ));
        assert_eq!(classify_in_repeat_read(bases, &quals, 1, 15), Some(b"CCCCGCCCCGCG".to_vec()));

        let bases: &[u8] = concat!(
            "GGGGCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCG",
            "GGGCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCGGGGCG"
        )
        .as_bytes();
        let quals = raw_quals(concat!(
            "((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((",
            "(((((((((((((((((((((((((((((((((((((((((((((((((((((((("
        ));
        assert_eq!(classify_in_repeat_read(bases, &quals, 1, 20), Some(b"CCCCGCCCCGCG".to_vec()));
    }

    #[test]
    fn classify_in_repeat_read_rejects_n_bases() {
        assert_eq!(classify_in_repeat_read(b"NNNNN", &raw_quals("$$$$$"), 1, 20), None);
    }
}
