use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

use expanse::commands::sample::{SampleArgs, run};

/// Each `#[test]` fn runs on its own dedicated thread under the default
/// harness, so keying the scratch directory by thread id (in addition to
/// process id) gives every test its own sandbox.
fn scratch_path(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "expanse-sample-test-{}-{:?}",
        std::process::id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir.join(name)
}

fn write_fixture(path: &PathBuf, header: Option<&str>, row_count: usize) {
    let mut file = File::create(path).unwrap();
    if let Some(header) = header {
        writeln!(file, "{header}").unwrap();
    }
    for i in 0..row_count {
        writeln!(file, "row{i}\t{i}").unwrap();
    }
}

fn default_args(input: PathBuf, output: PathBuf, confidence: f64, probability: f64) -> SampleArgs {
    SampleArgs {
        input,
        output,
        confidence,
        probability,
        no_header: false,
        seed: 42,
    }
}

fn read_lines(path: &PathBuf) -> Vec<String> {
    std::fs::read_to_string(path)
        .expect("output should be written")
        .lines()
        .map(str::to_string)
        .collect()
}

#[test]
fn sample_draws_formula_sized_sample_and_keeps_header_first() {
    let input_path = scratch_path("in.tsv");
    let output_path = scratch_path("out.tsv");
    write_fixture(&input_path, Some("id\tvalue"), 1000);

    let args = default_args(input_path, output_path.clone(), 0.95, 0.01);
    run(args).expect("sample run should succeed");

    let lines = read_lines(&output_path);
    // ceil(ln(1 - 0.95) / ln(1 - 0.01)) = ceil(ln(0.05)/ln(0.99)) = 299.
    assert_eq!(lines[0], "id\tvalue");
    assert_eq!(lines.len(), 300, "expected header + 299 sampled rows, got {}", lines.len());
}

#[test]
fn sample_is_reproducible_given_the_same_seed() {
    let input_path = scratch_path("repro_in.tsv");
    write_fixture(&input_path, Some("id\tvalue"), 500);

    let output_a = scratch_path("repro_a.tsv");
    run(default_args(input_path.clone(), output_a.clone(), 0.9, 0.05)).unwrap();

    let output_b = scratch_path("repro_b.tsv");
    run(default_args(input_path, output_b.clone(), 0.9, 0.05)).unwrap();

    assert_eq!(read_lines(&output_a), read_lines(&output_b));
}

#[test]
fn sample_different_seeds_draw_different_samples() {
    let input_path = scratch_path("seed_in.tsv");
    write_fixture(&input_path, Some("id\tvalue"), 500);

    let output_a = scratch_path("seed_a.tsv");
    run(default_args(input_path.clone(), output_a.clone(), 0.9, 0.05)).unwrap();

    let output_b = scratch_path("seed_b.tsv");
    let mut args_b = default_args(input_path, output_b.clone(), 0.9, 0.05);
    args_b.seed = 7;
    run(args_b).unwrap();

    assert_ne!(read_lines(&output_a), read_lines(&output_b));
}

#[test]
fn sample_preserves_original_file_order() {
    let input_path = scratch_path("order_in.tsv");
    let output_path = scratch_path("order_out.tsv");
    write_fixture(&input_path, Some("id\tvalue"), 500);

    run(default_args(input_path, output_path.clone(), 0.9, 0.05)).unwrap();

    let sampled_indices: Vec<usize> = read_lines(&output_path)[1..]
        .iter()
        .map(|line| line.strip_prefix("row").unwrap().split('\t').next().unwrap().parse().unwrap())
        .collect();

    let mut sorted = sampled_indices.clone();
    sorted.sort_unstable();
    assert_eq!(sampled_indices, sorted, "sampled lines should stay in original file order");
}

#[test]
fn sample_no_header_treats_every_line_as_data() {
    let input_path = scratch_path("no_header_in.tsv");
    let output_path = scratch_path("no_header_out.tsv");
    write_fixture(&input_path, None, 500);

    let mut args = default_args(input_path, output_path.clone(), 0.9, 0.05);
    args.no_header = true;
    run(args).unwrap();

    let lines = read_lines(&output_path);
    assert!(lines.iter().all(|line| line.starts_with("row")));
}

#[test]
fn sample_errors_when_requesting_more_lines_than_available() {
    let input_path = scratch_path("small_in.tsv");
    let output_path = scratch_path("small_out.tsv");
    write_fixture(&input_path, Some("id\tvalue"), 3);

    // n = ceil(ln(0.05)/ln(0.5)) = 5, more than the 3 available data rows.
    let result = run(default_args(input_path, output_path.clone(), 0.95, 0.5));

    assert!(result.is_err());
    assert!(!output_path.exists(), "no output should be written when there aren't enough lines");
}

#[test]
fn sample_rejects_out_of_range_confidence_and_probability() {
    let input_path = scratch_path("bounds_in.tsv");
    write_fixture(&input_path, Some("id\tvalue"), 10);

    let output_path = scratch_path("bounds_out.tsv");
    let bad_confidence = default_args(input_path.clone(), output_path.clone(), 1.0, 0.5);
    assert!(run(bad_confidence).is_err());

    let bad_probability = default_args(input_path, output_path, 0.9, 0.0);
    assert!(run(bad_probability).is_err());
}

#[test]
fn sample_requires_header_unless_no_header_given() {
    let input_path = scratch_path("empty_in.tsv");
    File::create(&input_path).unwrap();
    let output_path = scratch_path("empty_out.tsv");

    let result = run(default_args(input_path, output_path, 0.9, 0.5));
    assert!(result.is_err());
    assert!(result.unwrap_err().to_string().contains("--no-header"));
}
