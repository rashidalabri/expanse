use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::{Cigar, CigarString};
use rust_htslib::bam::{Format, Header, Record, Writer};

use expanse::commands::find_sinks::{FindSinksArgs, run};

const PAIRED: u16 = 1;
const UNMAP: u16 = 4;
const MUNMAP: u16 = 8;
const READ1: u16 = 64;
const READ2: u16 = 128;

/// A clean CAG-repeat read: passes the IRR purity filter.
fn irr_seq() -> Vec<u8> {
    "CAG".repeat(20).into_bytes()
}

/// A GATA-repeat read: passes the IRR purity filter with a different motif.
fn gata_seq() -> Vec<u8> {
    "GATA".repeat(15).into_bytes()
}

/// An aperiodic read: fails the IRR purity filter despite low MAPQ.
fn non_repetitive_seq() -> Vec<u8> {
    b"ACGTTGCAACGGTTCAGTAGCTAGCATCGATCGTAGCTAGGCTAGCATCGTAGCTAGCA".to_vec()
}

/// The canonical motif `expanse::irr::classify_in_repeat_read` assigns to
/// `seq`, used so assertions check against the real canonicalization
/// instead of a hardcoded (and possibly wrong) guess at its output.
fn canonical_motif(seq: &[u8]) -> String {
    let quals = vec![40u8; seq.len()];
    let motif =
        expanse::irr::classify_in_repeat_read(seq, &quals, 2, 20).expect("fixture sequence should classify as an IRR");
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
/// process id) gives every test its own sandbox.
fn scratch_path(name: &str) -> PathBuf {
    let dir =
        std::env::temp_dir().join(format!("expanse-find-sinks-test-{}-{:?}", std::process::id(), std::thread::current().id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir.join(name)
}

fn fixture_header() -> Header {
    let mut header = Header::new();
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", "chr1");
    sq.push_tag(b"LN", 10_000);
    header.push_record(&sq);
    header
}

/// A battery of read pairs designed to exercise the single-pass,
/// cache-based mate-resolution logic in `find_sinks::run` end to end:
///  - "pairA": IRR read1 (low mapq, CAG) arrives before its high-mapq
///    (anchor) mate -> one anchored sink at read1's own span.
///  - "pairB": the anchor mate arrives *before* the IRR read -> same
///    outcome, proving pairing is order-independent.
///  - "pairC": IRR read1 whose mate is flagged unmapped (plus the literal
///    unmapped mate record, which must be skipped entirely) -> unanchored,
///    resolved via the mate-unmapped fast path without ever touching the
///    pairing cache.
///  - "pairD": IRR read1 whose mate is mapped but only mapq 45 (between
///    max_irr_mapq and min_anchor_mapq) -> unanchored.
///  - "pairE": low-mapq but non-repetitive read1 with a high-mapq mate ->
///    excluded entirely (fails the IRR purity filter).
///  - "pairF": both mates are independently IRR (same CAG motif) with
///    overlapping spans -> both counted, unanchored, merged into one sink.
///  - "pairH"/"pairI": two separately-anchored IRR pairs whose own spans
///    overlap -> merge into a single sink region with irr_count 2,
///    anchored_count 2.
///  - "pairJ": IRR read1 whose mate is claimed mapped but never appears in
///    the file at all -> resolved as unanchored via the EOF cache flush.
///  - "pairK": an unpaired (non-PAIRED) IRR read -> unanchored, resolved
///    immediately without touching the cache.
///  - "pairL": an anchored IRR pair with a different (GATA) motif -> proves
///    the motif column and per-motif grouping.
fn fixture_records() -> Vec<Record> {
    let irr = irr_seq();
    let gata = gata_seq();
    let non_repetitive = non_repetitive_seq();
    let filler = vec![b'A'; 50];

    vec![
        // pairA: IRR-then-anchor.
        make_record("pairA", 0, 100, 10, PAIRED | READ1, 0, 5000, &irr),
        make_record("pairA", 0, 5000, 60, PAIRED | READ2, 0, 100, &filler),
        // pairB: anchor-then-IRR.
        make_record("pairB", 0, 2000, 60, PAIRED | READ1, 0, 2200, &filler),
        make_record("pairB", 0, 2200, 10, PAIRED | READ2, 0, 2000, &irr),
        // pairC: mate-unmapped fast path (plus the literal unmapped mate).
        make_record("pairC", 0, 3000, 10, PAIRED | READ1 | MUNMAP, -1, -1, &irr),
        make_record("pairC", 0, 3000, 0, PAIRED | READ2 | UNMAP, 0, 3000, &filler),
        // pairD: mid-mapq ("other") mate.
        make_record("pairD", 0, 4000, 10, PAIRED | READ1, 0, 4500, &irr),
        make_record("pairD", 0, 4500, 45, PAIRED | READ2, 0, 4000, &filler),
        // pairE: non-repetitive despite low mapq.
        make_record("pairE", 0, 4700, 8, PAIRED | READ1, 0, 4800, &non_repetitive),
        make_record("pairE", 0, 4800, 60, PAIRED | READ2, 0, 4700, &filler),
        // pairF: both mates independently IRR, overlapping spans.
        make_record("pairF", 0, 6000, 10, PAIRED | READ1, 0, 6030, &irr),
        make_record("pairF", 0, 6030, 12, PAIRED | READ2, 0, 6000, &irr),
        // pairH/pairI: two anchored IRR pairs with overlapping IRR spans.
        make_record("pairH", 0, 7000, 10, PAIRED | READ1, 0, 9000, &irr),
        make_record("pairH", 0, 9000, 55, PAIRED | READ2, 0, 7000, &filler),
        make_record("pairI", 0, 7040, 10, PAIRED | READ1, 0, 9100, &irr),
        make_record("pairI", 0, 9100, 55, PAIRED | READ2, 0, 7040, &filler),
        // pairJ: mate claimed mapped but never present in the file.
        make_record("pairJ", 0, 8000, 10, PAIRED | READ1, 0, 8500, &irr),
        // pairK: unpaired IRR read.
        make_record("pairK", 0, 8600, 10, 0, -1, -1, &irr),
        // pairL: anchored IRR pair with a different motif.
        make_record("pairL", 0, 9500, 10, PAIRED | READ1, 0, 9800, &gata),
        make_record("pairL", 0, 9800, 60, PAIRED | READ2, 0, 9500, &filler),
    ]
}

fn build_fixture_bam() -> PathBuf {
    let bam_path = scratch_path("fixture.bam");
    let header = fixture_header();

    let mut writer = Writer::from_path(&bam_path, &header, Format::Bam).unwrap();
    for record in &fixture_records() {
        writer.write(record).unwrap();
    }
    drop(writer);

    bam_path
}

struct BedRow {
    chrom: String,
    start: i64,
    end: i64,
    motif: String,
    irr_count: usize,
    anchored_count: usize,
}

fn read_bed(path: &PathBuf) -> Vec<BedRow> {
    let file = File::open(path).unwrap();
    BufReader::new(file)
        .lines()
        .map(|line| {
            let line = line.unwrap();
            let mut fields = line.split('\t');
            BedRow {
                chrom: fields.next().unwrap().to_string(),
                start: fields.next().unwrap().parse().unwrap(),
                end: fields.next().unwrap().parse().unwrap(),
                motif: fields.next().unwrap().to_string(),
                irr_count: fields.next().unwrap().parse().unwrap(),
                anchored_count: fields.next().unwrap().parse().unwrap(),
            }
        })
        .collect()
}

fn find_row(rows: &[BedRow], start: i64) -> &BedRow {
    rows.iter().find(|r| r.start == start).unwrap_or_else(|| panic!("no BED row starting at {start} in {:?}", rows_debug(rows)))
}

fn rows_debug(rows: &[BedRow]) -> Vec<(String, i64, i64, &str, usize, usize)> {
    rows.iter().map(|r| (r.chrom.clone(), r.start, r.end, r.motif.as_str(), r.irr_count, r.anchored_count)).collect()
}

#[test]
fn find_sinks_end_to_end() {
    let bam_path = build_fixture_bam();
    let output_path = scratch_path("sinks.bed");

    let args = FindSinksArgs {
        input: bam_path.to_str().unwrap().to_string(),
        output: output_path.clone(),
        max_irr_mapq: 40,
        min_anchor_mapq: 50,
        motif_min_len: 2,
        motif_max_len: 20,
        reference: None,
        threads: 1,
    };

    run(args).expect("find-sinks run should succeed");

    let rows = read_bed(&output_path);
    let cag_motif = canonical_motif(&irr_seq());
    let gata_motif = canonical_motif(&gata_seq());

    // pairA: IRR-then-anchor -> anchored.
    let row = find_row(&rows, 100);
    assert_eq!(row.chrom, "chr1");
    assert_eq!(row.end, 160);
    assert_eq!(row.motif, cag_motif);
    assert_eq!(row.irr_count, 1);
    assert_eq!(row.anchored_count, 1);

    // pairB: anchor-then-IRR -> anchored (order-independence).
    let row = find_row(&rows, 2200);
    assert_eq!(row.end, 2260);
    assert_eq!(row.irr_count, 1);
    assert_eq!(row.anchored_count, 1);

    // pairC: mate-unmapped fast path -> unanchored.
    let row = find_row(&rows, 3000);
    assert_eq!(row.irr_count, 1);
    assert_eq!(row.anchored_count, 0);

    // pairD: mid-mapq ("other") mate -> unanchored.
    let row = find_row(&rows, 4000);
    assert_eq!(row.irr_count, 1);
    assert_eq!(row.anchored_count, 0);

    // pairE: non-repetitive -> excluded entirely.
    assert!(rows.iter().all(|r| r.start != 4700), "pairE should not classify as IRR: {:?}", rows_debug(&rows));

    // pairF: both mates independently IRR, overlapping -> merged, unanchored.
    let row = find_row(&rows, 6000);
    assert_eq!(row.end, 6090, "pairF's overlapping IRR spans should merge into one region");
    assert_eq!(row.irr_count, 2);
    assert_eq!(row.anchored_count, 0);

    // pairH/pairI: two anchored IRR pairs with overlapping spans -> merged.
    let row = find_row(&rows, 7000);
    assert_eq!(row.end, 7100);
    assert_eq!(row.irr_count, 2);
    assert_eq!(row.anchored_count, 2);

    // pairJ: mate never appears in the file -> resolved unanchored at EOF.
    let row = find_row(&rows, 8000);
    assert_eq!(row.end, 8060);
    assert_eq!(row.irr_count, 1);
    assert_eq!(row.anchored_count, 0);

    // pairK: unpaired IRR read -> unanchored.
    let row = find_row(&rows, 8600);
    assert_eq!(row.end, 8660);
    assert_eq!(row.irr_count, 1);
    assert_eq!(row.anchored_count, 0);

    // pairL: anchored IRR pair, distinct (GATA) motif.
    let row = find_row(&rows, 9500);
    assert_eq!(row.end, 9560);
    assert_eq!(row.motif, gata_motif);
    assert_eq!(row.irr_count, 1);
    assert_eq!(row.anchored_count, 1);
    assert_ne!(cag_motif, gata_motif);

    // No spurious rows: exactly one row per distinct sink region above.
    assert_eq!(rows.len(), 9, "unexpected BED rows: {:?}", rows_debug(&rows));
}

#[test]
fn find_sinks_requires_reference_for_cram_input() {
    let output_path = scratch_path("sinks_cram.bed");
    let args = FindSinksArgs {
        input: "input.cram".to_string(),
        output: output_path,
        max_irr_mapq: 40,
        min_anchor_mapq: 50,
        motif_min_len: 2,
        motif_max_len: 20,
        reference: None,
        threads: 1,
    };

    let err = run(args).expect_err("CRAM input without --reference should fail");
    assert!(err.to_string().contains("--reference"));
}
