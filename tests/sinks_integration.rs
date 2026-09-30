use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::index::{self, Type};
use rust_htslib::bam::record::{Cigar, CigarString};
use rust_htslib::bam::{Format, Header, Record, Writer};

use expanse::commands::sinks::{SinksArgs, run};

const PAIRED: u16 = 1;
const UNMAP: u16 = 4;
const MUNMAP: u16 = 8;
const READ1: u16 = 64;
const READ2: u16 = 128;
const SECONDARY: u16 = 256;
const SUPPLEMENTARY: u16 = 2048;

/// A clean CAG-repeat read: passes the IRR purity filter.
fn cag_seq() -> Vec<u8> {
    "CAG".repeat(20).into_bytes()
}

/// A GATA-repeat read: passes the IRR purity filter with a different (4bp)
/// motif than `cag_seq`'s CAG (3bp).
fn gata_seq() -> Vec<u8> {
    "GATA".repeat(15).into_bytes()
}

/// An aperiodic read: fails the IRR purity filter despite low MAPQ.
fn non_repetitive_seq() -> Vec<u8> {
    b"ACGTTGCAACGGTTCAGTAGCTAGCATCGATCGTAGCTAGGCTAGCATCGTAGCTAGCA".to_vec()
}

/// The canonical motif `expanse::irr::identify_repeat_motifs` assigns to
/// `seq`, used so tests assert against the real canonicalization instead of
/// a hardcoded (and possibly wrong) guess at its output. `seq` is expected
/// to be a clean, single-motif fixture, so exactly one motif should come
/// back.
fn canonical_motif(seq: &[u8]) -> String {
    let quals = vec![40u8; seq.len()];
    let motifs = expanse::irr::identify_repeat_motifs(seq, &quals, 2, 20);
    assert_eq!(
        motifs.len(),
        1,
        "expected exactly one motif for this clean fixture, got {motifs:?}"
    );
    String::from_utf8(motifs.into_iter().next().unwrap()).expect("motif should be ASCII bases")
}

fn make_record(qname: &str, tid: i32, pos: i64, mapq: u8, flags: u16, seq: &[u8]) -> Record {
    let mut record = Record::new();
    let qual = vec![40u8; seq.len()];
    let cigar = CigarString(vec![Cigar::Match(seq.len() as u32)]);
    record.set(qname.as_bytes(), Some(&cigar), seq, &qual);
    record.set_tid(tid);
    record.set_pos(pos);
    record.set_mapq(mapq);
    record.set_flags(flags);
    record
}

/// A high-MAPQ (60), non-repetitive mate record for `qname`, sharing no
/// content requirements with its IRR partner since anchor classification
/// is MAPQ-only -- used to give a fixture IRR read a mate that clears
/// `--min-anchor-mapq` (default 50), which `sinks` now requires before an
/// IRR read is counted at all.
fn anchor_mate(qname: &str, tid: i32, pos: i64) -> Record {
    make_record(qname, tid, pos, 60, PAIRED | READ2, &non_repetitive_seq())
}

/// Each `#[test]` fn runs on its own dedicated thread under the default
/// harness, so keying the scratch directory by thread id (in addition to
/// process id) gives every test its own sandbox.
fn scratch_path(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "expanse-sinks-test-{}-{:?}",
        std::process::id(),
        std::thread::current().id()
    ));
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

/// A small unsorted-safe (but coordinate-sorted anyway) BAM on chr1 with:
///  - "cagNear1"/"cagNear2": low-mapq CAG-repeat reads, each given an
///    anchor mate (so they're counted at all -- see `anchor_mate`), whose
///    own spans overlap -> should merge into one region with irr_count 2.
///  - "cagFar": low-mapq CAG-repeat read (with an anchor mate) far from the
///    near cluster -> a separate region with irr_count 1.
///  - "gataSameSpot": low-mapq GATA-repeat read (with an anchor mate) at
///    the same position as "cagNear1" -> a distinct motif, so a separate
///    output row at (mostly) the same coordinates.
///  - "highMapq": high-mapq, otherwise IRR-classified -> excluded (not an
///    IRR candidate at all, regardless of pairing).
///  - "nonRepetitive": low-mapq, non-repetitive sequence -> excluded (not
///    an IRR, regardless of pairing).
///  - "unmapped"/"secondary"/"supplementary": low-mapq, IRR-classified, but
///    filtered out by their flags before motif classification even runs.
fn fixture_records() -> Vec<Record> {
    let cag = cag_seq();
    let gata = gata_seq();
    let non_repetitive = non_repetitive_seq();

    vec![
        make_record("cagNear1", 0, 100, 10, PAIRED | READ1, &cag),
        anchor_mate("cagNear1", 0, 8000),
        make_record("cagNear2", 0, 105, 10, PAIRED | READ1, &cag),
        anchor_mate("cagNear2", 0, 8010),
        make_record("cagFar", 0, 5000, 10, PAIRED | READ1, &cag),
        anchor_mate("cagFar", 0, 8020),
        make_record("gataSameSpot", 0, 100, 10, PAIRED | READ1, &gata),
        anchor_mate("gataSameSpot", 0, 8030),
        make_record("highMapq", 0, 200, 60, PAIRED, &cag),
        make_record("nonRepetitive", 0, 300, 8, PAIRED, &non_repetitive),
        make_record("unmapped", 0, 400, 10, PAIRED | UNMAP, &cag),
        make_record("secondary", 0, 500, 10, PAIRED | SECONDARY, &cag),
        make_record("supplementary", 0, 600, 10, PAIRED | SUPPLEMENTARY, &cag),
    ]
}

fn build_fixture_bam() -> PathBuf {
    let bam_path = scratch_path("fixture.bam");
    let header = fixture_header();

    let mut writer = Writer::from_path(&bam_path, &header, Format::Bam).unwrap();
    for record in &fixture_records() {
        writer.write(record).unwrap();
    }

    bam_path
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
    writeln!(fai_file, "chr1\t10000\t6\t10000\t10001").unwrap();

    fasta_path
}

fn build_fixture_cram() -> (PathBuf, PathBuf) {
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
    index::build(&cram_path, None, Type::Bai, 1).unwrap();

    (cram_path, reference_path)
}

struct BedRow {
    chrom: String,
    start: i64,
    end: i64,
    motif: String,
    irr_count: usize,
}

fn read_bed(path: &PathBuf) -> Vec<BedRow> {
    let text = std::fs::read_to_string(path).expect("output BED should be written");
    text.lines()
        .map(|line| {
            let mut fields = line.split('\t');
            BedRow {
                chrom: fields.next().unwrap().to_string(),
                start: fields.next().unwrap().parse().unwrap(),
                end: fields.next().unwrap().parse().unwrap(),
                motif: fields.next().unwrap().to_string(),
                irr_count: fields.next().unwrap().parse().unwrap(),
            }
        })
        .collect()
}

fn default_args(input: String, output: PathBuf, reference: Option<PathBuf>) -> SinksArgs {
    SinksArgs {
        input,
        output,
        max_irr_mapq: 40,
        min_anchor_mapq: 50,
        motif_min_len: 2,
        motif_max_len: 20,
        merge_distance: 0,
        reference,
        threads: 1,
    }
}

#[test]
fn sinks_merges_overlapping_same_motif_reads_and_separates_others() {
    let bam_path = build_fixture_bam();
    let output_path = scratch_path("sinks.bed");

    let args = default_args(
        bam_path.to_str().unwrap().to_string(),
        output_path.clone(),
        None,
    );
    run(args).expect("sinks run should succeed");

    let cag_motif = canonical_motif(&cag_seq());
    let gata_motif = canonical_motif(&gata_seq());
    assert_ne!(cag_motif, gata_motif, "fixture motifs should be distinct");

    let rows = read_bed(&output_path);
    assert_eq!(
        rows.len(),
        3,
        "expected 3 output rows, got {:?}",
        rows_debug(&rows)
    );

    // The merged near-CAG cluster: cagNear1 [100,160) + cagNear2 [105,165).
    let near_cag = rows
        .iter()
        .find(|r| r.motif == cag_motif && r.start == 100)
        .unwrap_or_else(|| panic!("expected a near CAG row, got {:?}", rows_debug(&rows)));
    assert_eq!(near_cag.chrom, "chr1");
    assert_eq!(near_cag.end, 165);
    assert_eq!(near_cag.irr_count, 2);

    // gataSameSpot: same coordinates as cagNear1 alone, different motif.
    let gata_row = rows
        .iter()
        .find(|r| r.motif == gata_motif)
        .unwrap_or_else(|| panic!("expected a GATA row, got {:?}", rows_debug(&rows)));
    assert_eq!(gata_row.chrom, "chr1");
    assert_eq!(gata_row.start, 100);
    assert_eq!(gata_row.end, 160);
    assert_eq!(gata_row.irr_count, 1);

    // The distant CAG read stays its own region.
    let far_cag = rows
        .iter()
        .find(|r| r.motif == cag_motif && r.start == 5000)
        .unwrap_or_else(|| panic!("expected a far CAG row, got {:?}", rows_debug(&rows)));
    assert_eq!(far_cag.chrom, "chr1");
    assert_eq!(far_cag.end, 5060);
    assert_eq!(far_cag.irr_count, 1);
}

fn rows_debug(rows: &[BedRow]) -> Vec<(String, i64, i64, String, usize)> {
    rows.iter()
        .map(|r| {
            (
                r.chrom.clone(),
                r.start,
                r.end,
                r.motif.clone(),
                r.irr_count,
            )
        })
        .collect()
}

#[test]
fn sinks_merge_distance_bridges_distant_same_motif_regions() {
    let bam_path = build_fixture_bam();
    let output_path = scratch_path("sinks_merged.bed");

    let mut args = default_args(
        bam_path.to_str().unwrap().to_string(),
        output_path.clone(),
        None,
    );
    // cagNear cluster ends at 165; cagFar starts at 5000. A merge distance
    // that bridges the ~4835bp gap should fold all three CAG reads into one
    // region.
    args.merge_distance = 5_000;
    run(args).expect("sinks run should succeed");

    let cag_motif = canonical_motif(&cag_seq());
    let rows = read_bed(&output_path);
    let cag_rows: Vec<&BedRow> = rows.iter().filter(|r| r.motif == cag_motif).collect();

    assert_eq!(
        cag_rows.len(),
        1,
        "expected all CAG regions bridged into one, got {:?}",
        rows_debug(&rows)
    );
    assert_eq!(cag_rows[0].start, 100);
    assert_eq!(cag_rows[0].end, 5060);
    assert_eq!(cag_rows[0].irr_count, 3);
}

#[test]
fn sinks_requires_reference_for_cram_input() {
    let (cram_path, _reference_path) = build_fixture_cram();
    let output_path = scratch_path("sinks_no_ref.bed");

    let args = default_args(cram_path.to_str().unwrap().to_string(), output_path, None);
    let result = run(args);

    assert!(
        result.is_err(),
        "expected an error without --reference for CRAM input"
    );
    assert!(result.unwrap_err().to_string().contains("--reference"));
}

/// A read whose composition is 20 A's followed by a single G, repeated:
/// qualifies under both the "A" homopolymer motif and the exact 21bp
/// repeat unit (see the equivalent fixture in `src/irr.rs`), but the 21bp
/// unit's exact match (score 1.0) beats the homopolymer's diluted score
/// (~0.90) -- so `sinks` should attribute the read to the 21bp motif only,
/// not fan it out across both.
fn mostly_a_with_rare_g_seq() -> Vec<u8> {
    let unit: Vec<u8> = (0..20).map(|_| b'A').chain(std::iter::once(b'G')).collect();
    unit.iter().cloned().cycle().take(21 * 16).collect()
}

#[test]
fn sinks_attributes_a_multi_motif_read_to_its_single_best_motif() {
    let bam_path = scratch_path("multi_motif_fixture.bam");
    let header = fixture_header();
    let seq = mostly_a_with_rare_g_seq();

    {
        let mut writer = Writer::from_path(&bam_path, &header, Format::Bam).unwrap();
        writer
            .write(&make_record("multiMotif", 0, 100, 10, PAIRED | READ1, &seq))
            .unwrap();
        writer.write(&anchor_mate("multiMotif", 0, 8000)).unwrap();
    }

    let output_path = scratch_path("multi_motif.bed");
    let mut args = default_args(
        bam_path.to_str().unwrap().to_string(),
        output_path.clone(),
        None,
    );
    args.motif_min_len = 1;
    args.motif_max_len = 30;
    run(args).expect("sinks run should succeed");

    let rows = read_bed(&output_path);
    assert_eq!(
        rows.len(),
        1,
        "expected the read to be attributed to a single motif, not fanned out across every \
         qualifying one: {:?}",
        rows_debug(&rows)
    );
    assert_eq!(rows[0].irr_count, 1);

    let quals = vec![40u8; seq.len()];
    let best_motif = expanse::irr::best_repeat_motif(&seq, &quals, 1, 30)
        .expect("fixture should qualify under at least one motif");
    assert_eq!(
        rows[0].motif.as_bytes(),
        best_motif.as_slice(),
        "expected the read attributed to its single best-scoring motif"
    );
}

/// A BAM exercising every pairing outcome `sinks` now has to distinguish:
///  - "anchoredIrr": low-mapq CAG IRR + a high-mapq (60) mate -> counted
///    (an anchored IRR pair).
///  - "pairedIrr": two low-mapq CAG IRR mates of each other -> both
///    counted (a paired IRR), even though their spans overlap and merge
///    into one output row.
///  - "droppedLowMapqMate": low-mapq CAG IRR whose mate is mapped but
///    neither an anchor (mapq 30, below the default 50) nor itself an IRR
///    (non-repetitive) -> dropped.
///  - "droppedNoMateInFile": low-mapq CAG IRR whose mate never appears in
///    the scan at all (unlike `sinks_counts_paired_irr_with_unmapped_mate`
///    below, no mate record is written here, not even an unmapped one) ->
///    left pending and dropped, since its pairing can never be resolved.
fn build_pairing_fixture_bam() -> PathBuf {
    let bam_path = scratch_path("pairing_fixture.bam");
    let header = fixture_header();
    let cag = cag_seq();
    let non_repetitive = non_repetitive_seq();

    let records = vec![
        make_record("anchoredIrr", 0, 2000, 10, PAIRED | READ1, &cag),
        anchor_mate("anchoredIrr", 0, 9000),
        make_record("pairedIrr", 0, 3000, 10, PAIRED | READ1, &cag),
        make_record("pairedIrr", 0, 3050, 10, PAIRED | READ2, &cag),
        make_record("droppedLowMapqMate", 0, 4000, 10, PAIRED | READ1, &cag),
        make_record(
            "droppedLowMapqMate",
            0,
            9100,
            30,
            PAIRED | READ2,
            &non_repetitive,
        ),
        make_record(
            "droppedNoMateInFile",
            0,
            5000,
            10,
            PAIRED | READ1 | MUNMAP,
            &cag,
        ),
    ];

    let mut writer = Writer::from_path(&bam_path, &header, Format::Bam).unwrap();
    for record in &records {
        writer.write(record).unwrap();
    }

    bam_path
}

#[test]
fn sinks_counts_anchored_and_paired_irrs_but_drops_unpaired_ones() {
    let bam_path = build_pairing_fixture_bam();
    let output_path = scratch_path("pairing.bed");

    let args = default_args(
        bam_path.to_str().unwrap().to_string(),
        output_path.clone(),
        None,
    );
    run(args).expect("sinks run should succeed");

    let cag_motif = canonical_motif(&cag_seq());
    let rows = read_bed(&output_path);

    // anchoredIrr (1) + both of pairedIrr's mates (2), whether or not
    // their overlapping spans end up merged into one row.
    let total_irr_count: usize = rows
        .iter()
        .filter(|r| r.motif == cag_motif)
        .map(|r| r.irr_count)
        .sum();
    assert_eq!(
        total_irr_count,
        3,
        "expected only the anchored IRR and both paired-IRR mates to be counted, got {:?}",
        rows_debug(&rows)
    );

    // droppedLowMapqMate (starting at 4000) and droppedNoMateInFile
    // (starting at 5000) should not appear at all.
    assert!(
        rows.iter().all(|r| r.start != 4000 && r.start != 5000),
        "reads with a disqualifying or missing mate should not appear in the output: {:?}",
        rows_debug(&rows)
    );
}

/// An IRR read paired with a mate that is genuinely unmapped (an actual
/// unmapped record in the file, unlike `droppedNoMateInFile` above where no
/// mate record exists at all) -- but that unmapped mate still
/// independently qualifies as an IRR by its own motif content, so the pair
/// should still be recognized as a paired IRR. Only the mapped side has a
/// position to report, so exactly one region is counted, not two.
#[test]
fn sinks_counts_paired_irr_with_unmapped_mate() {
    let bam_path = scratch_path("unmapped_mate_pairing_fixture.bam");
    let header = fixture_header();
    let cag = cag_seq();

    {
        let mut writer = Writer::from_path(&bam_path, &header, Format::Bam).unwrap();
        writer
            .write(&make_record(
                "unmappedMatePair",
                0,
                6000,
                10,
                PAIRED | READ1,
                &cag,
            ))
            .unwrap();
        writer
            .write(&make_record(
                "unmappedMatePair",
                0,
                6000,
                0,
                PAIRED | READ2 | UNMAP,
                &cag,
            ))
            .unwrap();
    }

    let output_path = scratch_path("unmapped_mate_pairing.bed");
    let args = default_args(
        bam_path.to_str().unwrap().to_string(),
        output_path.clone(),
        None,
    );
    run(args).expect("sinks run should succeed");

    let cag_motif = canonical_motif(&cag_seq());
    let rows = read_bed(&output_path);
    let cag_rows: Vec<&BedRow> = rows.iter().filter(|r| r.motif == cag_motif).collect();

    assert_eq!(
        cag_rows.len(),
        1,
        "expected the mapped side of the pair to be counted, got {:?}",
        rows_debug(&rows)
    );
    assert_eq!(cag_rows[0].start, 6000);
    assert_eq!(cag_rows[0].end, 6060);
    assert_eq!(
        cag_rows[0].irr_count, 1,
        "only the mapped mate has a position to report, even though both reads qualify as IRRs"
    );
}

#[test]
fn sinks_scans_cram_input() {
    let (cram_path, reference_path) = build_fixture_cram();
    let output_path = scratch_path("sinks_cram.bed");

    let args = default_args(
        cram_path.to_str().unwrap().to_string(),
        output_path.clone(),
        Some(reference_path),
    );
    run(args).expect("sinks run should succeed on CRAM input");

    let cag_motif = canonical_motif(&cag_seq());
    let rows = read_bed(&output_path);
    let far_cag = rows
        .iter()
        .find(|r| r.motif == cag_motif && r.start == 5000)
        .expect("expected the far CAG region to be present");
    assert_eq!(far_cag.end, 5060);
    assert_eq!(far_cag.irr_count, 1);
}
