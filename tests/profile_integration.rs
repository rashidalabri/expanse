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

/// A clean CAG-repeat read: passes the IRR purity filter.
fn irr_seq() -> Vec<u8> {
    "CAG".repeat(20).into_bytes()
}

/// A GATA-repeat read: passes the IRR purity filter with a different (4bp)
/// motif than `irr_seq`'s CAG (3bp).
fn gata_seq() -> Vec<u8> {
    "GATA".repeat(15).into_bytes()
}

/// An aperiodic read: fails the IRR purity filter despite low MAPQ.
fn non_repetitive_seq() -> Vec<u8> {
    b"ACGTTGCAACGGTTCAGTAGCTAGCATCGATCGTAGCTAGGCTAGCATCGTAGCTAGCA".to_vec()
}

/// The canonical motif `expanse::irr::classify_in_repeat_read` assigns to
/// `seq`, used so tests assert against the real canonicalization instead of
/// a hardcoded (and possibly wrong) guess at its output.
fn canonical_motif(seq: &[u8]) -> String {
    let quals = vec![40u8; seq.len()];
    let motif = expanse::irr::classify_in_repeat_read(
        seq,
        &quals,
        2,
        20,
        expanse::irr::DEFAULT_MAX_DEGENERATE_MONONUCLEOTIDE,
        expanse::irr::DEFAULT_MAX_DEGENERATE_DINUCLEOTIDE,
        expanse::irr::DEFAULT_MAX_DEGENERATE_TRINUCLEOTIDE,
        expanse::irr::DEFAULT_MAX_DEGENERATE_OTHER,
    )
    .expect("fixture sequence should classify as an IRR");
    String::from_utf8(motif).expect("motif should be ASCII bases")
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

/// Each `#[test]` fn runs on its own dedicated thread under the default
/// harness, so keying the scratch directory by thread id (in addition to
/// process id) gives every test its own sandbox: fixture builders that
/// reuse fixed names like "fixture.bam"/"fixture.bed" across several tests
/// (e.g. relying on htslib's own path-based `.bai`/`.fai`/`.crai`
/// conventions) can't race with each other when tests run in parallel.
fn scratch_path(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "expanse-profile-test-{}-{:?}",
        std::process::id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir.join(name)
}

/// Builds a small coordinate-sorted, indexed BAM on chr1 with several reads
/// designed to exercise the pass-1 candidate extraction logic end to end:
///  - "irrIn": low-mapq, IRR-classified (CAG-repeat sequence), inside the
///    BED region, mate mapped at 5000 -> should end up in the --output BAM
///    and anchor its mate's location (not its own) in the summary.
///  - "irrUnmappedMate": low-mapq, IRR-classified, inside the BED region,
///    but its mate is unmapped -> still written to --output, but has no
///    anchor to contribute to the summary.
///  - "highMapq": high-mapq, inside the BED region -> excluded by the MAPQ
///    filter.
///  - "nonRepetitive": low-mapq, inside the BED region, but its sequence is
///    NOT repetitive -> excluded by the IRR purity filter.
///  - "outsideRegion": low-mapq, IRR-classified, entirely outside the BED
///    region -> never seen by the region-restricted fetch, excluded.
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

    vec![
        make_record("irrIn", 0, 100, 10, PAIRED | READ1, 0, 5000, &irr),
        make_record("irrUnmappedMate", 0, 105, 12, PAIRED | READ1 | MUNMAP, -1, -1, &irr),
        make_record("highMapq", 0, 120, 60, PAIRED | READ1, 0, 6000, &irr),
        make_record("nonRepetitive", 0, 140, 8, PAIRED | READ1, 0, 7000, &non_repetitive),
        make_record("outsideRegion", 0, 9000, 5, PAIRED | READ1, 0, 9500, &irr),
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
/// so the CRAM fetch path in `profile::run` gets exercised, not just
/// BAM/BAI.
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

fn read_output_qnames(path: &PathBuf) -> Vec<String> {
    let mut reader = Reader::from_path(path).unwrap();
    let mut out = Vec::new();
    for result in reader.records() {
        let record = result.unwrap();
        out.push(String::from_utf8_lossy(record.qname()).to_string());
    }
    out
}

#[test]
fn profile_extracts_irr_candidates() {
    let (bam_path, bed_path) = build_fixture_bam();
    let output_path = scratch_path("output.bam");
    let summary_path = scratch_path("summary.json");

    let args = ProfileArgs {
        bed: bed_path,
        input: bam_path.to_str().unwrap().to_string(),
        summary: summary_path.clone(),
        output: Some(output_path.clone()),
        max_irr_mapq: 40,
        motif_min_len: 2,
        motif_max_len: 20,
        max_degenerate_mononucleotide: expanse::irr::DEFAULT_MAX_DEGENERATE_MONONUCLEOTIDE,
        max_degenerate_dinucleotide: expanse::irr::DEFAULT_MAX_DEGENERATE_DINUCLEOTIDE,
        max_degenerate_trinucleotide: expanse::irr::DEFAULT_MAX_DEGENERATE_TRINUCLEOTIDE,
        max_degenerate_other: expanse::irr::DEFAULT_MAX_DEGENERATE_OTHER,
        anchor_merge_distance: 500,
        read_length: 150,
        reference: None,
        output_format: None,
        threads: 1,
    };

    run(args).expect("profile run should succeed");

    let written = read_output_qnames(&output_path);
    let mut counts: HashMap<String, usize> = HashMap::new();
    for qname in &written {
        *counts.entry(qname.clone()).or_default() += 1;
    }

    assert_eq!(written.len(), 2, "expected exactly 2 records, got {written:?}");
    assert_eq!(counts.get("irrIn"), Some(&1), "low-mapq IRR read inside the region should be written");
    assert_eq!(
        counts.get("irrUnmappedMate"),
        Some(&1),
        "an unmapped mate should no longer exclude an IRR candidate"
    );
    assert!(!counts.contains_key("highMapq"), "high-mapq read should be excluded");
    assert!(
        !counts.contains_key("nonRepetitive"),
        "low-mapq but non-repetitive read should be excluded by the IRR purity filter"
    );
    assert!(!counts.contains_key("outsideRegion"), "out-of-region read should be excluded");

    // The summary reports the *anchor* (mate) location, not the IRR read's
    // own location: only "irrIn" has a mapped mate (at 5000), so it's the
    // only entry, spanning [mate_pos, mate_pos + read_length).
    let summary_text = std::fs::read_to_string(&summary_path).expect("summary JSON should be written");
    let summary: serde_json::Value = serde_json::from_str(&summary_text).expect("summary should be valid JSON");
    let regions = summary.as_array().expect("summary should be a JSON array");
    assert_eq!(
        regions.len(),
        1,
        "only irrIn has a mapped mate to anchor on: {summary:#}"
    );
    assert_eq!(regions[0]["chrom"], "chr1");
    assert_eq!(regions[0]["start"], 5000);
    assert_eq!(regions[0]["end"], 5150);
    assert_eq!(regions[0]["irr_count"], 1);
}

/// With `--output` omitted, no alignment file should be written at all,
/// even though the run otherwise succeeds and still produces the JSON
/// summary.
#[test]
fn profile_writes_no_bam_output_by_default() {
    let (bam_path, bed_path) = build_fixture_bam();
    let summary_path = scratch_path("summary_no_output.json");

    let args = ProfileArgs {
        bed: bed_path,
        input: bam_path.to_str().unwrap().to_string(),
        summary: summary_path.clone(),
        output: None,
        max_irr_mapq: 40,
        motif_min_len: 2,
        motif_max_len: 20,
        max_degenerate_mononucleotide: expanse::irr::DEFAULT_MAX_DEGENERATE_MONONUCLEOTIDE,
        max_degenerate_dinucleotide: expanse::irr::DEFAULT_MAX_DEGENERATE_DINUCLEOTIDE,
        max_degenerate_trinucleotide: expanse::irr::DEFAULT_MAX_DEGENERATE_TRINUCLEOTIDE,
        max_degenerate_other: expanse::irr::DEFAULT_MAX_DEGENERATE_OTHER,
        anchor_merge_distance: 500,
        read_length: 150,
        reference: None,
        output_format: None,
        threads: 1,
    };

    run(args).expect("profile run should succeed");

    assert!(summary_path.exists(), "summary JSON should always be written");
}

/// Same scenario as `profile_extracts_irr_candidates`, but against a CRAM
/// input, exercising the indexed per-region CRAM fetch path end to end.
#[test]
fn profile_extracts_irr_candidates_cram() {
    let (cram_path, reference_path, bed_path) = build_fixture_cram();
    let output_path = scratch_path("output_cram.bam");
    let summary_path = scratch_path("summary_cram.json");

    let args = ProfileArgs {
        bed: bed_path,
        input: cram_path.to_str().unwrap().to_string(),
        summary: summary_path.clone(),
        output: Some(output_path.clone()),
        max_irr_mapq: 40,
        motif_min_len: 2,
        motif_max_len: 20,
        max_degenerate_mononucleotide: expanse::irr::DEFAULT_MAX_DEGENERATE_MONONUCLEOTIDE,
        max_degenerate_dinucleotide: expanse::irr::DEFAULT_MAX_DEGENERATE_DINUCLEOTIDE,
        max_degenerate_trinucleotide: expanse::irr::DEFAULT_MAX_DEGENERATE_TRINUCLEOTIDE,
        max_degenerate_other: expanse::irr::DEFAULT_MAX_DEGENERATE_OTHER,
        anchor_merge_distance: 500,
        read_length: 150,
        reference: Some(reference_path),
        output_format: Some(OutputFormat::Bam),
        threads: 1,
    };

    run(args).expect("profile run should succeed");

    let written = read_output_qnames(&output_path);
    let mut counts: HashMap<String, usize> = HashMap::new();
    for qname in &written {
        *counts.entry(qname.clone()).or_default() += 1;
    }

    assert_eq!(written.len(), 2, "expected exactly 2 records, got {written:?}");
    assert_eq!(counts.get("irrIn"), Some(&1), "low-mapq IRR read inside the region should be written");
    assert_eq!(
        counts.get("irrUnmappedMate"),
        Some(&1),
        "an unmapped mate should no longer exclude an IRR candidate"
    );
    assert!(!counts.contains_key("highMapq"), "high-mapq read should be excluded");
    assert!(
        !counts.contains_key("nonRepetitive"),
        "low-mapq but non-repetitive read should be excluded by the IRR purity filter"
    );
    assert!(!counts.contains_key("outsideRegion"), "out-of-region read should be excluded");

    let summary_text = std::fs::read_to_string(&summary_path).expect("summary JSON should be written");
    let summary: serde_json::Value = serde_json::from_str(&summary_text).expect("summary should be valid JSON");
    let regions = summary.as_array().expect("summary should be a JSON array");
    assert_eq!(
        regions.len(),
        1,
        "only irrIn has a mapped mate to anchor on: {summary:#}"
    );
    assert_eq!(regions[0]["chrom"], "chr1");
    assert_eq!(regions[0]["start"], 5000);
    assert_eq!(regions[0]["end"], 5150);
    assert_eq!(regions[0]["irr_count"], 1);
}

/// Four candidate reads, all with their own position inside one BED region
/// (50-150, so a single fetch picks all of them up), but whose *mate*
/// (anchor) locations define two separate clusters:
///  - "motifA_1"/"motifA_2" (CAG) and "motifB_1" (GATA) anchor at 3000,
///    3005, 3010 -- close enough that their (default 150bp `--read-length`)
///    spans overlap and merge into a single near anchor region spanning
///    [3000, 3160), with counts broken down 2 (CAG) / 1 (GATA).
///  - "farCluster_1" (CAG) anchors at 4160, 1000bp past the near cluster's
///    merged end (3160) -- far enough that the default 500bp
///    `--anchor-merge-distance` keeps it as its own anchor region ([4160,
///    4310)), but a large enough distance (3000) merges it into the near
///    one.
fn summary_fixture_records() -> Vec<Record> {
    let cag = irr_seq();
    let gata = gata_seq();

    vec![
        make_record("motifA_1", 0, 60, 10, PAIRED | READ1, 0, 3000, &cag),
        make_record("motifA_2", 0, 65, 10, PAIRED | READ1, 0, 3005, &cag),
        make_record("motifB_1", 0, 70, 10, PAIRED | READ1, 0, 3010, &gata),
        make_record("farCluster_1", 0, 75, 10, PAIRED | READ1, 0, 4160, &cag),
    ]
}

fn build_summary_fixture_bam() -> (PathBuf, PathBuf) {
    let bam_path = scratch_path("summary_fixture.bam");
    let header = fixture_header();

    {
        let mut writer = Writer::from_path(&bam_path, &header, Format::Bam).unwrap();
        for record in &summary_fixture_records() {
            writer.write(record).unwrap();
        }
    }

    index::build(&bam_path, None, Type::Bai, 1).unwrap();

    (bam_path, build_fixture_bed())
}

#[test]
fn profile_summary_keeps_distant_anchors_separate_by_default() {
    let (bam_path, bed_path) = build_summary_fixture_bam();
    let summary_path = scratch_path("summary_distant.json");

    let args = ProfileArgs {
        bed: bed_path,
        input: bam_path.to_str().unwrap().to_string(),
        summary: summary_path.clone(),
        output: None,
        max_irr_mapq: 40,
        motif_min_len: 2,
        motif_max_len: 20,
        max_degenerate_mononucleotide: expanse::irr::DEFAULT_MAX_DEGENERATE_MONONUCLEOTIDE,
        max_degenerate_dinucleotide: expanse::irr::DEFAULT_MAX_DEGENERATE_DINUCLEOTIDE,
        max_degenerate_trinucleotide: expanse::irr::DEFAULT_MAX_DEGENERATE_TRINUCLEOTIDE,
        max_degenerate_other: expanse::irr::DEFAULT_MAX_DEGENERATE_OTHER,
        anchor_merge_distance: 500,
        read_length: 150,
        reference: None,
        output_format: None,
        threads: 1,
    };

    run(args).expect("profile run should succeed");

    let cag_motif = canonical_motif(&irr_seq());
    let gata_motif = canonical_motif(&gata_seq());
    assert_ne!(cag_motif, gata_motif, "fixture motifs should be distinct");

    let summary_text = std::fs::read_to_string(&summary_path).expect("summary JSON should be written");
    let summary: serde_json::Value = serde_json::from_str(&summary_text).expect("summary should be valid JSON");
    let mut regions = summary.as_array().expect("summary should be a JSON array").clone();
    regions.sort_by_key(|region| region["start"].as_i64().unwrap());

    assert_eq!(regions.len(), 2, "expected two separate anchor regions, got {summary:#}");

    let near = &regions[0];
    assert_eq!(near["chrom"], "chr1");
    assert_eq!(near["start"], 3000);
    assert_eq!(near["end"], 3160);
    assert_eq!(near["irr_count"], 3);
    assert_eq!(near["motifs"][&cag_motif], 2, "expected 2 CAG-motif IRRs: {summary:#}");
    assert_eq!(near["motifs"][&gata_motif], 1, "expected 1 GATA-motif IRR: {summary:#}");

    let far = &regions[1];
    assert_eq!(far["chrom"], "chr1");
    assert_eq!(far["start"], 4160);
    assert_eq!(far["end"], 4310);
    assert_eq!(far["irr_count"], 1);
    assert_eq!(far["motifs"][&cag_motif], 1);
}

#[test]
fn profile_summary_merges_anchors_within_custom_distance() {
    let (bam_path, bed_path) = build_summary_fixture_bam();
    let summary_path = scratch_path("summary_merged.json");

    let args = ProfileArgs {
        bed: bed_path,
        input: bam_path.to_str().unwrap().to_string(),
        summary: summary_path.clone(),
        output: None,
        max_irr_mapq: 40,
        motif_min_len: 2,
        motif_max_len: 20,
        max_degenerate_mononucleotide: expanse::irr::DEFAULT_MAX_DEGENERATE_MONONUCLEOTIDE,
        max_degenerate_dinucleotide: expanse::irr::DEFAULT_MAX_DEGENERATE_DINUCLEOTIDE,
        max_degenerate_trinucleotide: expanse::irr::DEFAULT_MAX_DEGENERATE_TRINUCLEOTIDE,
        max_degenerate_other: expanse::irr::DEFAULT_MAX_DEGENERATE_OTHER,
        anchor_merge_distance: 3000,
        read_length: 150,
        reference: None,
        output_format: None,
        threads: 1,
    };

    run(args).expect("profile run should succeed");

    let cag_motif = canonical_motif(&irr_seq());
    let gata_motif = canonical_motif(&gata_seq());

    let summary_text = std::fs::read_to_string(&summary_path).expect("summary JSON should be written");
    let summary: serde_json::Value = serde_json::from_str(&summary_text).expect("summary should be valid JSON");
    let regions = summary.as_array().expect("summary should be a JSON array");

    assert_eq!(
        regions.len(),
        1,
        "expected the two clusters to merge into one anchor region, got {summary:#}"
    );
    let region = &regions[0];
    assert_eq!(region["chrom"], "chr1");
    assert_eq!(region["start"], 3000);
    assert_eq!(region["end"], 4310);
    assert_eq!(region["irr_count"], 4);
    assert_eq!(region["motifs"][&cag_motif], 3, "expected 3 CAG-motif IRRs: {summary:#}");
    assert_eq!(region["motifs"][&gata_motif], 1, "expected 1 GATA-motif IRR: {summary:#}");
}

/// A read whose composition is 20 A's followed by a single G, repeated:
/// mostly-A with only rare G interruptions passes the homopolymer "A"
/// motif's thresholds despite not being a pure homopolymer, while the exact
/// 21bp repeat unit also passes on its own -- so
/// `irr::classify_in_repeat_read_all` legitimately returns two distinct
/// motifs for it (see the equivalent unit test in `src/irr.rs` for why a
/// pure homopolymer doesn't trigger this).
fn mostly_a_with_rare_g_seq() -> Vec<u8> {
    let unit: Vec<u8> = (0..20).map(|_| b'A').chain(std::iter::once(b'G')).collect();
    unit.iter().cloned().cycle().take(21 * 16).collect()
}

/// A single IRR read that qualifies under two distinct motifs should be
/// counted once in its anchor region's `irr_count`, but once under each of
/// the two motifs in the `motifs` breakdown.
#[test]
fn profile_summary_counts_multi_motif_read_once_per_motif_not_per_read() {
    let bam_path = scratch_path("multi_motif_fixture.bam");
    let header = fixture_header();
    let seq = mostly_a_with_rare_g_seq();

    {
        let mut writer = Writer::from_path(&bam_path, &header, Format::Bam).unwrap();
        writer.write(&make_record("multiMotif", 0, 60, 10, PAIRED | READ1, 0, 5000, &seq)).unwrap();
    }
    index::build(&bam_path, None, Type::Bai, 1).unwrap();

    let bed_path = scratch_path("multi_motif_fixture.bed");
    let mut bed_file = File::create(&bed_path).unwrap();
    writeln!(bed_file, "chr1\t50\t{}", 60 + seq.len() as i64 + 50).unwrap();

    let summary_path = scratch_path("summary_multi_motif.json");

    let args = ProfileArgs {
        bed: bed_path,
        input: bam_path.to_str().unwrap().to_string(),
        summary: summary_path.clone(),
        output: None,
        max_irr_mapq: 40,
        motif_min_len: 1,
        motif_max_len: 30,
        max_degenerate_mononucleotide: expanse::irr::DEFAULT_MAX_DEGENERATE_MONONUCLEOTIDE,
        max_degenerate_dinucleotide: expanse::irr::DEFAULT_MAX_DEGENERATE_DINUCLEOTIDE,
        max_degenerate_trinucleotide: expanse::irr::DEFAULT_MAX_DEGENERATE_TRINUCLEOTIDE,
        max_degenerate_other: expanse::irr::DEFAULT_MAX_DEGENERATE_OTHER,
        anchor_merge_distance: 500,
        read_length: 150,
        reference: None,
        output_format: None,
        threads: 1,
    };

    run(args).expect("profile run should succeed");

    let summary_text = std::fs::read_to_string(&summary_path).expect("summary JSON should be written");
    let summary: serde_json::Value = serde_json::from_str(&summary_text).expect("summary should be valid JSON");
    let regions = summary.as_array().expect("summary should be a JSON array");

    assert_eq!(regions.len(), 1, "expected a single anchor region, got {summary:#}");
    let region = &regions[0];

    assert_eq!(
        region["irr_count"], 1,
        "a single read must only be counted once toward the region total, even though it \
         qualifies under multiple motifs: {summary:#}"
    );

    // Not asserting an exact motif count: with degenerate (IUPAC) motif
    // calling, this fixture's confident, real G interruptions can also
    // independently satisfy other periods via a partially degenerate
    // (R-containing) motif. That doesn't matter for what this test checks
    // -- every motif entry a single read contributes to must show count 1,
    // never more, no matter how many entries there are.
    let motifs = region["motifs"].as_object().expect("motifs should be a JSON object");
    assert!(
        motifs.len() >= 2,
        "expected at least the read's two originally-intended qualifying motifs: {summary:#}"
    );
    for (motif, count) in motifs {
        assert_eq!(count, 1, "motif {motif:?} should have exactly 1 IRR: {summary:#}");
    }
}

/// A motif containing an IUPAC ambiguity code (e.g. an `R` for a
/// consistently purine-mixed position) should flow untouched through the
/// whole pipeline -- dedup, clustering, and JSON serialization -- and show
/// up as a plain ASCII string key in the `--summary` output, with no
/// `profile.rs` code needing to know anything about the wider alphabet.
#[test]
fn profile_summary_reports_iupac_ambiguity_code_in_motif() {
    let bam_path = scratch_path("iupac_motif_fixture.bam");
    let header = fixture_header();
    let seq = mostly_a_with_rare_g_seq();

    {
        let mut writer = Writer::from_path(&bam_path, &header, Format::Bam).unwrap();
        writer.write(&make_record("iupacMotif", 0, 60, 10, PAIRED | READ1, 0, 5000, &seq)).unwrap();
    }
    index::build(&bam_path, None, Type::Bai, 1).unwrap();

    let bed_path = scratch_path("iupac_motif_fixture.bed");
    let mut bed_file = File::create(&bed_path).unwrap();
    writeln!(bed_file, "chr1\t50\t{}", 60 + seq.len() as i64 + 50).unwrap();

    let summary_path = scratch_path("summary_iupac_motif.json");

    let args = ProfileArgs {
        bed: bed_path,
        input: bam_path.to_str().unwrap().to_string(),
        summary: summary_path.clone(),
        output: None,
        max_irr_mapq: 40,
        motif_min_len: 1,
        motif_max_len: 30,
        max_degenerate_mononucleotide: expanse::irr::DEFAULT_MAX_DEGENERATE_MONONUCLEOTIDE,
        max_degenerate_dinucleotide: expanse::irr::DEFAULT_MAX_DEGENERATE_DINUCLEOTIDE,
        max_degenerate_trinucleotide: expanse::irr::DEFAULT_MAX_DEGENERATE_TRINUCLEOTIDE,
        max_degenerate_other: expanse::irr::DEFAULT_MAX_DEGENERATE_OTHER,
        anchor_merge_distance: 500,
        read_length: 150,
        reference: None,
        output_format: None,
        threads: 1,
    };

    run(args).expect("profile run should succeed");

    let summary_text = std::fs::read_to_string(&summary_path).expect("summary JSON should be written");
    let summary: serde_json::Value = serde_json::from_str(&summary_text).expect("summary should be valid JSON");
    let regions = summary.as_array().expect("summary should be a JSON array");
    let motifs = regions[0]["motifs"].as_object().expect("motifs should be a JSON object");

    let is_ambiguous = |motif: &str| motif.bytes().any(|b| !matches!(b, b'A' | b'C' | b'G' | b'T'));
    assert!(
        motifs.keys().any(|motif| is_ambiguous(motif)),
        "expected at least one motif with an IUPAC ambiguity code in the summary: {summary:#}"
    );
}
