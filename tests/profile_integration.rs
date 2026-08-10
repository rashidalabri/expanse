use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::index::{self, Type};
use rust_htslib::bam::record::{Cigar, CigarString};
use rust_htslib::bam::{Format, Header, Read, Reader, Record, Writer};

use expanse::commands::profile::{run, OutputFormat, ProfileArgs};

const PAIRED: u16 = 1;
const MUNMAP: u16 = 8;
const READ1: u16 = 64;
const READ2: u16 = 128;

/// A clean CAG-repeat read: passes the IRR purity filter.
fn irr_seq() -> Vec<u8> {
    "CAG".repeat(20).into_bytes()
}

/// An aperiodic read: fails the IRR purity filter despite low MAPQ.
fn non_repetitive_seq() -> Vec<u8> {
    b"ACGTTGCAACGGTTCAGTAGCTAGCATCGATCGTAGCTAGGCTAGCATCGTAGCTAGCA".to_vec()
}

#[allow(clippy::too_many_arguments)]
fn make_record(qname: &str, tid: i32, pos: i64, mapq: u8, flags: u16, mtid: i32, mpos: i64, seq: &[u8]) -> Record {
    let mut record = Record::new();
    let qual = vec![40u8; seq.len()];
    let cigar = CigarString(vec![Cigar::Match(seq.len() as u32)]);
    record.set(qname.as_bytes(), Some(&cigar), seq, &qual);
    record.set_tid(tid);
    record.set_pos(pos);
    record.set_mapq(mapq);
    record.set_mtid(mtid);
    record.set_mpos(mpos);
    record.set_flags(flags);
    record
}

fn scratch_path(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("expanse-profile-test-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir.join(name)
}

/// Builds a small coordinate-sorted, indexed BAM on chr1 with several read
/// pairs designed to exercise the two-pass extraction logic end to end:
///  - "pairA": read1 is low-mapq, IRR-classified (CAG-repeat sequence)
///    inside the BED region; mate is high-mapq far away -> both should end
///    up in the output (pass 1 + pass 2).
///  - "pairB": both ends high-mapq inside the BED region -> excluded.
///  - "pairC": read1 is low-mapq, IRR-classified, inside the BED region,
///    but its mate is unmapped -> excluded entirely in pass 1 (an IRR read
///    with no mate to anchor it is never a candidate).
///  - "pairD": low-mapq, IRR-classified read1 entirely outside the BED
///    region -> never seen by pass 1 (region-restricted fetch), excluded.
///  - "pairE": read1 is low-mapq and inside the BED region, but its
///    sequence is NOT repetitive -> must be excluded by the IRR purity
///    filter despite passing the MAPQ threshold, and its mate must not be
///    fetched either.
///  - "pairF": read1 is low-mapq, IRR-classified, inside the BED region,
///    but its mate's own MAPQ is below `min_anchor_mapq` -> the mate is
///    excluded from pass 2, and read1 itself must then be dropped in the
///    final filtering step since it has no passing anchor.
///  - "decoyX": high-mapq read that merely overlaps pairA's merged mate
///    region -> must NOT be pulled in by pass 2 (qname filter).
fn fixture_header() -> Header {
    let mut header = Header::new();
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", "chr1");
    sq.push_tag(b"LN", 10_000);
    header.push_record(&sq);
    header
}

fn fixture_records() -> Vec<Record> {
    let irr = irr_seq();
    let non_repetitive = non_repetitive_seq();
    let filler = vec![b'A'; 50];

    vec![
        make_record("pairA", 0, 100, 10, PAIRED | READ1, 0, 5000, &irr),
        make_record("pairF", 0, 105, 12, PAIRED | READ1, 0, 8000, &irr),
        make_record("pairB", 0, 120, 60, PAIRED | READ1, 0, 6000, &irr),
        make_record("pairC", 0, 130, 5, PAIRED | READ1 | MUNMAP, -1, -1, &irr),
        make_record("pairE", 0, 140, 8, PAIRED | READ1, 0, 7000, &non_repetitive),
        make_record("decoyX", 0, 4950, 60, PAIRED | READ1, -1, -1, &[b'A'; 100]),
        make_record("pairA", 0, 5000, 60, PAIRED | READ2, 0, 100, &filler),
        make_record("pairB", 0, 6000, 60, PAIRED | READ2, 0, 120, &filler),
        make_record("pairE", 0, 7000, 60, PAIRED | READ2, 0, 140, &filler),
        make_record("pairF", 0, 8000, 30, PAIRED | READ2, 0, 105, &filler),
        make_record("pairD", 0, 9000, 5, PAIRED | READ1, 0, 9500, &irr),
        make_record("pairD", 0, 9500, 60, PAIRED | READ2, 0, 9000, &filler),
    ]
}

fn build_fixture_bed() -> PathBuf {
    let bed_path = scratch_path("fixture.bed");
    let mut bed_file = File::create(&bed_path).unwrap();
    writeln!(bed_file, "chr1\t50\t150").unwrap();
    bed_path
}

fn build_fixture_bam() -> (PathBuf, PathBuf) {
    let bam_path = scratch_path("fixture.bam");
    let header = fixture_header();

    {
        let mut writer = Writer::from_path(&bam_path, &header, Format::Bam).unwrap();
        for record in &fixture_records() {
            writer.write(record).unwrap();
        }
    }

    index::build(&bam_path, None, Type::Bai, 1).unwrap();

    (bam_path, build_fixture_bed())
}

/// Writes a single-contig reference FASTA (with a hand-written `.fai`)
/// covering `chr1:1-10000`, for use as the CRAM reference.
fn build_fixture_reference() -> PathBuf {
    let fasta_path = scratch_path("fixture.fa");
    let sequence = vec![b'A'; 10_000];

    let mut fasta_file = File::create(&fasta_path).unwrap();
    writeln!(fasta_file, ">chr1").unwrap();
    fasta_file.write_all(&sequence).unwrap();
    fasta_file.write_all(b"\n").unwrap();

    let fai_path = scratch_path("fixture.fa.fai");
    let mut fai_file = File::create(&fai_path).unwrap();
    // name, length, byte offset of sequence, bases per line, bytes per line
    writeln!(fai_file, "chr1\t10000\t6\t10000\t10001").unwrap();

    fasta_path
}

/// Same fixture data as `build_fixture_bam`, but written as an indexed CRAM
/// (with a real `.crai`) so the CRAM-slice-merging path in `profile::run`
/// gets exercised against htslib's own index format, not just BAM/BAI.
fn build_fixture_cram() -> (PathBuf, PathBuf, PathBuf) {
    let cram_path = scratch_path("fixture.cram");
    let reference_path = build_fixture_reference();
    let header = fixture_header();

    {
        let mut writer = Writer::from_path(&cram_path, &header, Format::Cram).unwrap();
        writer.set_reference(&reference_path).unwrap();
        for record in &fixture_records() {
            writer.write(record).unwrap();
        }
    }

    // sam_index_build3 (which this calls into) dispatches on the file's
    // actual format, so this produces a `.crai` despite the `Type::Bai` hint.
    index::build(&cram_path, None, Type::Bai, 1).unwrap();

    (cram_path, reference_path, build_fixture_bed())
}

fn read_output_qnames(path: &PathBuf) -> Vec<(String, u16)> {
    let mut reader = Reader::from_path(path).unwrap();
    let mut out = Vec::new();
    for result in reader.records() {
        let record = result.unwrap();
        out.push((String::from_utf8_lossy(record.qname()).to_string(), record.flags()));
    }
    out
}

#[test]
fn profile_extracts_irr_reads_and_their_mates() {
    let (bam_path, bed_path) = build_fixture_bam();
    let output_path = scratch_path("output.bam");

    let args = ProfileArgs {
        bed: bed_path,
        input: bam_path.to_str().unwrap().to_string(),
        output: output_path.clone(),
        max_irr_mapq: 40,
        motif_min_len: 2,
        motif_max_len: 20,
        min_anchor_mapq: 50,
        reference: None,
        output_format: None,
        threads: 1,
    };

    run(args).expect("profile run should succeed");

    let written = read_output_qnames(&output_path);

    let mut counts: HashMap<String, usize> = HashMap::new();
    for (qname, _) in &written {
        *counts.entry(qname.clone()).or_default() += 1;
    }

    assert_eq!(written.len(), 2, "expected exactly 2 records, got {written:?}");
    assert_eq!(counts.get("pairA"), Some(&2), "both pairA reads should be present");
    assert!(!counts.contains_key("pairB"), "high-mapq pairB should be excluded");
    assert!(!counts.contains_key("pairC"), "IRR read with an unmapped mate should be excluded entirely");
    assert!(!counts.contains_key("pairD"), "out-of-region pairD should be excluded");
    assert!(
        !counts.contains_key("pairE"),
        "low-mapq but non-repetitive pairE should be excluded by the IRR purity filter, \
         and its mate should never be looked up"
    );
    assert!(!counts.contains_key("decoyX"), "decoy overlapping the mate region must not be pulled in");
    assert!(
        !counts.contains_key("pairF"),
        "pairF's read1 must be dropped in the final step since its mate never passed the anchor filter"
    );

    let pair_a_flags: Vec<u16> = written.iter().filter(|(q, _)| q == "pairA").map(|(_, f)| *f).collect();
    assert!(pair_a_flags.contains(&(PAIRED | READ1)));
    assert!(pair_a_flags.contains(&(PAIRED | READ2)));
}

/// Same scenario as `profile_extracts_irr_reads_and_their_mates`, but against
/// a CRAM input, exercising the CRAI-backed slice-merging fetch path (real
/// `hopen`/`hread2` FFI reads of an htslib-built `.crai`, region expansion,
/// and the post-fetch region-membership filter) end to end.
#[test]
fn profile_extracts_irr_reads_and_their_mates_cram() {
    let (cram_path, reference_path, bed_path) = build_fixture_cram();
    let output_path = scratch_path("output_cram.bam");

    // Confirm the CRAI-backed path actually engages (rather than
    // profile::run silently falling back to per-region fetch on a load
    // error) before checking output correctness below.
    let slices = expanse::crai::load(cram_path.to_str().unwrap()).expect("CRAI index should load");
    assert!(!slices.is_empty(), "expected at least one CRAM slice in the index");

    let args = ProfileArgs {
        bed: bed_path,
        input: cram_path.to_str().unwrap().to_string(),
        output: output_path.clone(),
        max_irr_mapq: 40,
        motif_min_len: 2,
        motif_max_len: 20,
        min_anchor_mapq: 50,
        reference: Some(reference_path),
        output_format: Some(OutputFormat::Bam),
        threads: 1,
    };

    run(args).expect("profile run should succeed");

    let written = read_output_qnames(&output_path);

    let mut counts: HashMap<String, usize> = HashMap::new();
    for (qname, _) in &written {
        *counts.entry(qname.clone()).or_default() += 1;
    }

    assert_eq!(written.len(), 2, "expected exactly 2 records, got {written:?}");
    assert_eq!(counts.get("pairA"), Some(&2), "both pairA reads should be present");
    assert!(!counts.contains_key("pairB"), "high-mapq pairB should be excluded");
    assert!(!counts.contains_key("pairC"), "IRR read with an unmapped mate should be excluded entirely");
    assert!(!counts.contains_key("pairD"), "out-of-region pairD should be excluded");
    assert!(
        !counts.contains_key("pairE"),
        "low-mapq but non-repetitive pairE should be excluded by the IRR purity filter, \
         and its mate should never be looked up"
    );
    assert!(!counts.contains_key("decoyX"), "decoy overlapping the mate region must not be pulled in");
    assert!(
        !counts.contains_key("pairF"),
        "pairF's read1 must be dropped in the final step since its mate never passed the anchor filter"
    );

    let pair_a_flags: Vec<u16> = written.iter().filter(|(q, _)| q == "pairA").map(|(_, f)| *f).collect();
    assert!(pair_a_flags.contains(&(PAIRED | READ1)));
    assert!(pair_a_flags.contains(&(PAIRED | READ2)));
}
