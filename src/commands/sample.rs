//! Samples `N` lines at random from a line-delimited input file, where `N`
//! is the standard "at least one success" sample size:
//! `N = ceil(ln(1 - c) / ln(1 - p))`, the number of samples needed for at
//! least confidence `c` of observing, at least once, an event that occurs
//! independently with per-sample probability `p`. Here that's sized to
//! give `c` confidence of observing an expansion of motif `M` at least
//! once, given `M` is expanded in a `p` fraction of samples.
//!
//! Sampling is done in a single streaming pass with reservoir sampling
//! (Algorithm R), so memory use is bounded by `N` rather than the input
//! file's size, and the file never needs to be read twice or loaded whole.
//!
//! The header line (unless `--no-header`) is never counted toward `N` and
//! is always copied to the output first; sampled data lines are written in
//! their original file order, not the order they were drawn in.
//!
//! It is an error for the input to have fewer than `N` data lines: no
//! output is written in that case.

use std::fs::File;
use std::io::{BufRead, BufReader, Write as IoWrite};
use std::path::PathBuf;

use anyhow::{Context, Result, bail};
use clap::Args;
use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};

#[derive(Args, Debug)]
pub struct SampleArgs {
    /// Line-delimited input file to sample from.
    #[arg(short = 'i', long)]
    pub input: PathBuf,

    /// Sampled output path: the header line (unless `--no-header`)
    /// followed by the sampled data lines, in their original file order.
    #[arg(short = 'o', long)]
    pub output: PathBuf,

    /// Confidence `c`, in `[0, 1)`, of observing an expansion of motif `M`
    /// at least once across the sample.
    #[arg(short = 'c', long)]
    pub confidence: f64,

    /// Probability `p`, in `(0, 1)`, that motif `M` is expanded in a given
    /// sample.
    #[arg(short = 'p', long)]
    pub probability: f64,

    /// Treat the input as headerless: every line is eligible for
    /// sampling. By default the first line is treated as a header,
    /// excluded from sampling, and always copied to the output first.
    #[arg(long)]
    pub no_header: bool,

    /// Random seed for the sample draw. Fixed by default so runs are
    /// reproducible given the same input and `N`; override to draw a
    /// different sample.
    #[arg(long, default_value_t = 42)]
    pub seed: u64,
}

/// `N = ceil(ln(1 - c) / ln(1 - p))`: the number of independent samples
/// needed for at least confidence `c` of observing, at least once, an
/// event with per-sample probability `p`. Computed via `ln_1p` (`ln(1+x)`)
/// rather than `ln` directly, for better precision when `c` or `p` is close
/// to zero.
fn sample_size(confidence: f64, probability: f64) -> Result<usize> {
    if !(0.0..1.0).contains(&confidence) {
        bail!("--confidence must be in [0, 1), got {confidence}");
    }
    if !(probability > 0.0 && probability < 1.0) {
        bail!("--probability must be in (0, 1), got {probability}");
    }

    let n = ((-confidence).ln_1p() / (-probability).ln_1p()).ceil();
    Ok(n as usize)
}

pub fn run(args: SampleArgs) -> Result<()> {
    let n = sample_size(args.confidence, args.probability)?;

    let file =
        File::open(&args.input).with_context(|| format!("failed to open input {:?}", args.input))?;
    let mut lines = BufReader::new(file).lines();

    let header = if args.no_header {
        None
    } else {
        match lines.next() {
            Some(line) => Some(
                line.with_context(|| format!("failed to read header line of {:?}", args.input))?,
            ),
            None => bail!(
                "input {:?} is empty; expected a header line (pass --no-header if it has none)",
                args.input
            ),
        }
    };

    let mut rng = StdRng::seed_from_u64(args.seed);

    // Reservoir sampling (Algorithm R): a single streaming pass that draws
    // a uniform-random sample of `n` lines without needing to know the
    // total line count up front, keeping memory use bounded by `n` rather
    // than the input's size. `line_no` (0-based, over data lines only) is
    // kept alongside each sampled line so the reservoir can be restored to
    // original file order before writing.
    let mut reservoir: Vec<(usize, String)> = Vec::with_capacity(n);
    let mut seen: usize = 0;
    for (line_no, line) in lines.enumerate() {
        let line = line
            .with_context(|| format!("failed to read {:?} at data line {}", args.input, line_no + 1))?;
        seen += 1;

        if reservoir.len() < n {
            reservoir.push((line_no, line));
        } else {
            let j = rng.random_range(0..seen);
            if j < n {
                reservoir[j] = (line_no, line);
            }
        }
    }

    if seen < n {
        bail!(
            "sample: requested {n} lines (c={}, p={}) but {:?} only has {seen} data line{}",
            args.confidence,
            args.probability,
            args.input,
            if seen == 1 { "" } else { "s" },
        );
    }

    reservoir.sort_by_key(|(line_no, _)| *line_no);

    let mut writer = std::io::BufWriter::new(
        File::create(&args.output)
            .with_context(|| format!("failed to create output {:?}", args.output))?,
    );
    if let Some(header) = &header {
        writeln!(writer, "{header}").context("failed to write header line")?;
    }
    for (_, line) in &reservoir {
        writeln!(writer, "{line}").context("failed to write sampled line")?;
    }
    writer
        .flush()
        .with_context(|| format!("failed to flush output {:?}", args.output))?;

    log::info!(
        "sample: drew {} of {seen} data lines from {:?} (n={n}, c={}, p={}, seed={}), written to \
         {:?}",
        reservoir.len(),
        args.input,
        args.confidence,
        args.probability,
        args.seed,
        args.output,
    );

    Ok(())
}
