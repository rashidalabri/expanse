use std::fs;
use std::path::PathBuf;

use expanse::commands::merge::{MergeArgs, run};

fn scratch_path(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "expanse-merge-test-{}-{:?}",
        std::process::id(),
        std::thread::current().id()
    ));
    fs::create_dir_all(&dir).unwrap();
    dir.join(name)
}

fn write_summary(name: &str, json: &str) -> PathBuf {
    let path = scratch_path(name);
    fs::write(&path, json).unwrap();
    path
}

fn write_manifest(name: &str, entries: &[(&str, &PathBuf)]) -> PathBuf {
    let path = scratch_path(name);
    let contents: String = entries
        .iter()
        .map(|(sample_id, summary_path)| format!("{sample_id}\t{}\n", summary_path.display()))
        .collect();
    fs::write(&path, contents).unwrap();
    path
}

fn read_merged(path: &PathBuf) -> Vec<serde_json::Value> {
    let text = fs::read_to_string(path).expect("merged summary should be written");
    let mut value: Vec<serde_json::Value> =
        serde_json::from_str(&text).expect("merged summary should be valid JSON");
    value.sort_by_key(|locus| locus["start"].as_i64().unwrap());
    value
}

/// Two samples with overlapping-ish loci on the same locus should merge
/// into one, with `irr_count`/`motifs` broken out per sample.
#[test]
fn merge_combines_nearby_loci_across_samples() {
    let sample_a = write_summary(
        "a.json",
        r#"[{"chrom":"chr1","start":1000,"end":1100,"irr_count":5,"motifs":{"CAG":5}}]"#,
    );
    let sample_b = write_summary(
        "b.json",
        r#"[{"chrom":"chr1","start":1090,"end":1200,"irr_count":3,"motifs":{"CAG":2,"GATA":1}}]"#,
    );
    let manifest = write_manifest("manifest_basic.tsv", &[("sampleA", &sample_a), ("sampleB", &sample_b)]);
    let output = scratch_path("merged_basic.json");

    run(MergeArgs {
        manifest,
        output: output.clone(),
        merge_distance: 500,
        batch_size: 100,
    })
    .expect("merge run should succeed");

    let merged = read_merged(&output);
    assert_eq!(merged.len(), 1, "expected one merged locus: {merged:#?}");
    let locus = &merged[0];
    assert_eq!(locus["chrom"], "chr1");
    assert_eq!(locus["start"], 1000);
    assert_eq!(locus["end"], 1200);
    assert_eq!(locus["irr_count"]["sampleA"], 5);
    assert_eq!(locus["irr_count"]["sampleB"], 3);
    assert_eq!(locus["motifs"]["CAG"]["sampleA"], 5);
    assert_eq!(locus["motifs"]["CAG"]["sampleB"], 2);
    assert_eq!(locus["motifs"]["GATA"]["sampleB"], 1);
}

/// Loci far enough apart should stay as separate merged loci.
#[test]
fn merge_keeps_distant_loci_separate() {
    let sample_a = write_summary(
        "far_a.json",
        r#"[{"chrom":"chr1","start":1000,"end":1100,"irr_count":1,"motifs":{"CAG":1}}]"#,
    );
    let sample_b = write_summary(
        "far_b.json",
        r#"[{"chrom":"chr1","start":5000,"end":5100,"irr_count":1,"motifs":{"CAG":1}}]"#,
    );
    let manifest = write_manifest("manifest_far.tsv", &[("sampleA", &sample_a), ("sampleB", &sample_b)]);
    let output = scratch_path("merged_far.json");

    run(MergeArgs {
        manifest,
        output: output.clone(),
        merge_distance: 500,
        batch_size: 100,
    })
    .expect("merge run should succeed");

    let merged = read_merged(&output);
    assert_eq!(merged.len(), 2, "expected two separate loci: {merged:#?}");
    assert_eq!(merged[0]["start"], 1000);
    assert_eq!(merged[1]["start"], 5000);
}

/// With `--batch-size 1`, each sample is folded into the running merge one
/// at a time. Three samples' loci are individually more than
/// `--merge-distance` apart from their non-adjacent neighbor (A to C is
/// 1300bp), but each is within distance of the next (A-B gap is exactly
/// 500bp, B-C gap is exactly 500bp), so processing them one batch at a time
/// must still chain all three into a single merged locus -- proving the
/// iterative batching doesn't lose transitive merges across batch
/// boundaries.
#[test]
fn merge_chains_transitively_across_batch_boundaries() {
    let sample_a = write_summary(
        "chain_a.json",
        r#"[{"chrom":"chr1","start":1000,"end":1100,"irr_count":1,"motifs":{"CAG":1}}]"#,
    );
    let sample_b = write_summary(
        "chain_b.json",
        r#"[{"chrom":"chr1","start":1600,"end":1700,"irr_count":2,"motifs":{"CAG":2}}]"#,
    );
    let sample_c = write_summary(
        "chain_c.json",
        r#"[{"chrom":"chr1","start":2200,"end":2300,"irr_count":4,"motifs":{"CAG":4}}]"#,
    );
    let manifest = write_manifest(
        "manifest_chain.tsv",
        &[("sampleA", &sample_a), ("sampleB", &sample_b), ("sampleC", &sample_c)],
    );
    let output = scratch_path("merged_chain.json");

    run(MergeArgs {
        manifest,
        output: output.clone(),
        merge_distance: 500,
        batch_size: 1,
    })
    .expect("merge run should succeed");

    let merged = read_merged(&output);
    assert_eq!(
        merged.len(),
        1,
        "expected all three samples to chain into one merged locus: {merged:#?}"
    );
    let locus = &merged[0];
    assert_eq!(locus["start"], 1000);
    assert_eq!(locus["end"], 2300);
    assert_eq!(locus["irr_count"]["sampleA"], 1);
    assert_eq!(locus["irr_count"]["sampleB"], 2);
    assert_eq!(locus["irr_count"]["sampleC"], 4);
    assert_eq!(locus["motifs"]["CAG"]["sampleA"], 1);
    assert_eq!(locus["motifs"]["CAG"]["sampleB"], 2);
    assert_eq!(locus["motifs"]["CAG"]["sampleC"], 4);
}

/// The same three-sample fixture as the chaining test above, but run with
/// a single batch covering all samples at once -- the result should be
/// identical regardless of `--batch-size`.
#[test]
fn merge_result_is_independent_of_batch_size() {
    let sample_a = write_summary(
        "onebatch_a.json",
        r#"[{"chrom":"chr1","start":1000,"end":1100,"irr_count":1,"motifs":{"CAG":1}}]"#,
    );
    let sample_b = write_summary(
        "onebatch_b.json",
        r#"[{"chrom":"chr1","start":1600,"end":1700,"irr_count":2,"motifs":{"CAG":2}}]"#,
    );
    let sample_c = write_summary(
        "onebatch_c.json",
        r#"[{"chrom":"chr1","start":2200,"end":2300,"irr_count":4,"motifs":{"CAG":4}}]"#,
    );
    let manifest = write_manifest(
        "manifest_onebatch.tsv",
        &[("sampleA", &sample_a), ("sampleB", &sample_b), ("sampleC", &sample_c)],
    );
    let output = scratch_path("merged_onebatch.json");

    run(MergeArgs {
        manifest,
        output: output.clone(),
        merge_distance: 500,
        batch_size: 100,
    })
    .expect("merge run should succeed");

    let merged = read_merged(&output);
    assert_eq!(merged.len(), 1, "expected one merged locus: {merged:#?}");
    let locus = &merged[0];
    assert_eq!(locus["start"], 1000);
    assert_eq!(locus["end"], 2300);
    assert_eq!(locus["irr_count"]["sampleA"], 1);
    assert_eq!(locus["irr_count"]["sampleB"], 2);
    assert_eq!(locus["irr_count"]["sampleC"], 4);
}

/// Distinct contigs must never be bridged together, no matter how the
/// coordinate ranges compare numerically.
#[test]
fn merge_does_not_bridge_across_contigs() {
    let sample_a = write_summary(
        "contig_a.json",
        r#"[{"chrom":"chr1","start":1000,"end":1100,"irr_count":1,"motifs":{"CAG":1}}]"#,
    );
    let sample_b = write_summary(
        "contig_b.json",
        r#"[{"chrom":"chr2","start":1050,"end":1150,"irr_count":1,"motifs":{"CAG":1}}]"#,
    );
    let manifest = write_manifest(
        "manifest_contig.tsv",
        &[("sampleA", &sample_a), ("sampleB", &sample_b)],
    );
    let output = scratch_path("merged_contig.json");

    run(MergeArgs {
        manifest,
        output: output.clone(),
        merge_distance: 500,
        batch_size: 100,
    })
    .expect("merge run should succeed");

    let merged = read_merged(&output);
    assert_eq!(merged.len(), 2, "expected loci on different contigs to stay separate: {merged:#?}");
}

/// Same sample, two of its own loci that start out more than
/// `--merge-distance` apart (so they land in separate merged loci after the
/// first batch) but get bridged together in a later batch by another
/// sample's intervening locus. The merge must recognize both pieces as the
/// same sample within the resulting single locus and sum their
/// `irr_count`/`motifs`, rather than only keeping one.
#[test]
fn merge_sums_same_sample_contributions_bridged_across_batches() {
    let sample_a = write_summary(
        "bridge_a.json",
        r#"[
            {"chrom":"chr1","start":1000,"end":1100,"irr_count":1,"motifs":{"CAG":1}},
            {"chrom":"chr1","start":2200,"end":2300,"irr_count":4,"motifs":{"CAG":4}}
        ]"#,
    );
    let sample_b = write_summary(
        "bridge_b.json",
        r#"[{"chrom":"chr1","start":1600,"end":1700,"irr_count":2,"motifs":{"GATA":2}}]"#,
    );
    let manifest = write_manifest(
        "manifest_bridge.tsv",
        &[("sampleA", &sample_a), ("sampleB", &sample_b)],
    );
    let output = scratch_path("merged_bridge.json");

    run(MergeArgs {
        manifest,
        output: output.clone(),
        merge_distance: 500,
        batch_size: 1,
    })
    .expect("merge run should succeed");

    let merged = read_merged(&output);
    assert_eq!(
        merged.len(),
        1,
        "expected sampleA's two loci and sampleB's bridging locus to merge into one: {merged:#?}"
    );
    let locus = &merged[0];
    assert_eq!(locus["start"], 1000);
    assert_eq!(locus["end"], 2300);
    assert_eq!(
        locus["irr_count"]["sampleA"], 5,
        "sampleA's two separate contributions (1 + 4) should be summed: {merged:#?}"
    );
    assert_eq!(locus["motifs"]["CAG"]["sampleA"], 5);
    assert_eq!(locus["irr_count"]["sampleB"], 2);
    assert_eq!(locus["motifs"]["GATA"]["sampleB"], 2);
}
