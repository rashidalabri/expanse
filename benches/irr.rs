//! Benchmarks for `expanse::irr::classify_in_repeat_read`, covering both
//! ends of its early-exit behavior (a clean repeat vs. a read with no
//! periodic structure), a noisy-but-still-IRR read, and a read whose
//! confident, real minority allele (not noise) exercises the IUPAC
//! degenerate-motif consensus path, at a couple of representative read
//! lengths.

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use std::hint::black_box;

use expanse::irr::{
    DEFAULT_MAX_DEGENERATE_DINUCLEOTIDE, DEFAULT_MAX_DEGENERATE_MONONUCLEOTIDE, DEFAULT_MAX_DEGENERATE_OTHER,
    DEFAULT_MAX_DEGENERATE_TRINUCLEOTIDE, DEFAULT_MOTIF_MAX_LEN, DEFAULT_MOTIF_MIN_LEN, classify_in_repeat_read,
};

/// Test qualities are given as ASCII PHRED+33 (as if lifted straight from a
/// FASTQ/SAM text field); convert to the raw PHRED bytes `classify_in_repeat_read`
/// expects.
fn raw_quals(ascii: &str) -> Vec<u8> {
    ascii.bytes().map(|c| c - 33).collect()
}

/// Builds a `read_len`-byte read by repeating `unit`, plus a matching
/// high-quality (Q40) qualities vector.
fn repeat_read(unit: &[u8], read_len: usize) -> (Vec<u8>, Vec<u8>) {
    let bases: Vec<u8> = unit.iter().copied().cycle().take(read_len).collect();
    let quals = vec![b'I' - 33; read_len]; // Q40
    (bases, quals)
}

/// Builds a `read_len`-byte read with no periodic structure (a de Bruijn-ish
/// shuffle of the four bases), plus matching Q40 qualities.
fn random_read(read_len: usize) -> (Vec<u8>, Vec<u8>) {
    const ALPHABET: [u8; 4] = [b'A', b'C', b'G', b'T'];
    // A small LCG is enough to avoid accidental periodicity without pulling
    // in a `rand` dependency just for benchmark fixtures.
    let mut state: u64 = 0x2545F4914F6CDD1D;
    let bases: Vec<u8> = (0..read_len)
        .map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            ALPHABET[(state >> 62) as usize]
        })
        .collect();
    let quals = vec![b'I' - 33; read_len]; // Q40
    (bases, quals)
}

/// Builds a `read_len`-byte read that is 20 A's followed by a single G,
/// repeated, plus matching Q40 qualities. The rare G is confident (not
/// noise), so it exercises `extract_consensus_base_iupac`'s tiered
/// coverage check on every period tried, the same fixture shape used in
/// `src/irr.rs`'s and `tests/profile_integration.rs`'s degenerate-motif
/// tests.
fn degenerate_read(read_len: usize) -> (Vec<u8>, Vec<u8>) {
    let unit: Vec<u8> = (0..20).map(|_| b'A').chain(std::iter::once(b'G')).collect();
    repeat_read(&unit, read_len)
}

fn bench_classify_in_repeat_read(c: &mut Criterion) {
    let mut group = c.benchmark_group("classify_in_repeat_read");

    for &read_len in &[150usize, 1000] {
        let (clean_bases, clean_quals) = repeat_read(b"CAG", read_len);
        group.bench_with_input(BenchmarkId::new("clean_repeat", read_len), &read_len, |b, _| {
            b.iter(|| {
                classify_in_repeat_read(
                    black_box(&clean_bases),
                    black_box(&clean_quals),
                    DEFAULT_MOTIF_MIN_LEN,
                    DEFAULT_MOTIF_MAX_LEN,
                    DEFAULT_MAX_DEGENERATE_MONONUCLEOTIDE,
                    DEFAULT_MAX_DEGENERATE_DINUCLEOTIDE,
                    DEFAULT_MAX_DEGENERATE_TRINUCLEOTIDE,
                    DEFAULT_MAX_DEGENERATE_OTHER,
                )
            })
        });

        let (nonrepeat_bases, nonrepeat_quals) = random_read(read_len);
        group.bench_with_input(BenchmarkId::new("non_repeat", read_len), &read_len, |b, _| {
            b.iter(|| {
                classify_in_repeat_read(
                    black_box(&nonrepeat_bases),
                    black_box(&nonrepeat_quals),
                    DEFAULT_MOTIF_MIN_LEN,
                    DEFAULT_MOTIF_MAX_LEN,
                    DEFAULT_MAX_DEGENERATE_MONONUCLEOTIDE,
                    DEFAULT_MAX_DEGENERATE_DINUCLEOTIDE,
                    DEFAULT_MAX_DEGENERATE_TRINUCLEOTIDE,
                    DEFAULT_MAX_DEGENERATE_OTHER,
                )
            })
        });

        let (degenerate_bases, degenerate_quals) = degenerate_read(read_len);
        group.bench_with_input(BenchmarkId::new("degenerate_repeat", read_len), &read_len, |b, _| {
            b.iter(|| {
                classify_in_repeat_read(
                    black_box(&degenerate_bases),
                    black_box(&degenerate_quals),
                    DEFAULT_MOTIF_MIN_LEN,
                    DEFAULT_MOTIF_MAX_LEN,
                    DEFAULT_MAX_DEGENERATE_MONONUCLEOTIDE,
                    DEFAULT_MAX_DEGENERATE_DINUCLEOTIDE,
                    DEFAULT_MAX_DEGENERATE_TRINUCLEOTIDE,
                    DEFAULT_MAX_DEGENERATE_OTHER,
                )
            })
        });
    }

    // A real-world-shaped noisy IRR: mostly a AAATG repeat with a low-quality
    // tail run of A/T homopolymer, matching the kind of read the pass-1 scan
    // in `commands::profile` actually has to classify.
    let noisy_bases: &[u8] = concat!(
        "TCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTT",
        "CATTTCATTTCATTTCATTTCATTTCTTTTTTTTTATTTTTTTTTATTTTATATCGGAT"
    )
    .as_bytes();
    let noisy_quals = raw_quals(concat!(
        "((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((",
        "(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((("
    ));
    group.bench_function("noisy_repeat", |b| {
        b.iter(|| {
            classify_in_repeat_read(
                black_box(noisy_bases),
                black_box(&noisy_quals),
                DEFAULT_MOTIF_MIN_LEN,
                DEFAULT_MOTIF_MAX_LEN,
                DEFAULT_MAX_DEGENERATE_MONONUCLEOTIDE,
                DEFAULT_MAX_DEGENERATE_DINUCLEOTIDE,
                DEFAULT_MAX_DEGENERATE_TRINUCLEOTIDE,
                DEFAULT_MAX_DEGENERATE_OTHER,
            )
        })
    });

    group.finish();
}

criterion_group!(benches, bench_classify_in_repeat_read);
criterion_main!(benches);
