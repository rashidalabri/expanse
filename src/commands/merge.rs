//! Merges `profile --summary` outputs from many samples into one combined
//! summary: loci within `--merge-distance` bp of each other (across *all*
//! samples, not just within one) are folded into a single merged locus,
//! whose `irr_count` (`{sample_id: count}`) and `motifs`
//! (`{motif: {sample_id: count}}`) break each contributing sample's own
//! counts out individually.
//!
//! Manifest samples are folded in `--batch-size` at a time rather than all
//! at once, so peak memory is bounded by one batch's parsed summaries plus
//! the running merged-locus set, not by the full manifest.
//!
//! If `--depths` is given, every count is normalized to `--target-depth`
//! via `(count * target_depth) / sample_depth` before being written out;
//! otherwise counts are written raw and a warning is logged.

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use clap::Args;
use serde::{Deserialize, Serialize};

use crate::bed::{self, Region};

#[derive(Args, Debug)]
pub struct MergeArgs {
    /// Headerless TSV manifest: one sample per line, `sample_id<TAB>path`,
    /// where `path` is a `profile --summary` JSON output for that sample.
    #[arg(long)]
    pub manifest: PathBuf,

    /// Merged JSON summary output path.
    #[arg(long)]
    pub output: PathBuf,

    /// Merge loci within this many bp of each other, across all samples,
    /// into one merged locus.
    #[arg(long, default_value_t = 500)]
    pub merge_distance: i64,

    /// Number of samples to fold into the running merge at a time, so
    /// memory usage stays bounded regardless of manifest size.
    #[arg(long, default_value_t = 100)]
    pub batch_size: usize,

    /// Optional headerless TSV of per-sample sequencing depths:
    /// `sample_id<TAB>depth` (float). When given, every sample's IRR
    /// counts are normalized to `--target-depth` before being written out.
    /// When omitted, counts are written raw (unnormalized) and a warning
    /// is logged.
    #[arg(long)]
    pub depths: Option<PathBuf>,

    /// Depth that IRR counts are normalized to when `--depths` is given:
    /// `(count * target_depth) / sample_depth`. Has no effect without
    /// `--depths`.
    #[arg(long, default_value_t = 30.0)]
    pub target_depth: f64,
}

/// One entry of a `profile --summary` input file. Mirrors
/// `commands::profile::AnchorRegionSummary`'s shape.
#[derive(Deserialize, Debug)]
struct AnchorRegionSummaryInput {
    chrom: String,
    start: i64,
    end: i64,
    irr_count: usize,
    motifs: BTreeMap<String, usize>,
}

/// One sample's contribution to a merged locus: its own IRR count and
/// motif breakdown, unmodified from its `profile --summary` input.
#[derive(Debug, Clone)]
struct SampleLocusSummary {
    irr_count: usize,
    motifs: BTreeMap<String, usize>,
}

/// A locus formed by merging one or more samples' anchor regions that fall
/// within `--merge-distance` bp of each other, keyed by sample id. This is
/// the running/internal representation carried across batches; see
/// [`MergedLocusOutput`] for the pivoted shape actually written out.
struct MergedLocus {
    chrom: String,
    start: i64,
    end: i64,
    samples: BTreeMap<String, SampleLocusSummary>,
}

/// One entry in the merge output. `irr_count` and `motifs` are pivoted from
/// [`MergedLocus`]'s per-sample map onto per-sample maps of their own:
/// `irr_count` is `{sample_id: count}`, and `motifs` is
/// `{motif: {sample_id: count}}`. Counts are depth-normalized (see
/// [`normalize_count`]) when a depths file was given.
#[derive(Serialize, Debug)]
struct MergedLocusOutput {
    chrom: String,
    start: i64,
    end: i64,
    irr_count: BTreeMap<String, f64>,
    motifs: BTreeMap<String, BTreeMap<String, f64>>,
}

/// Normalizes a raw IRR count for `sample_id` to `target_depth`, i.e.
/// `(count * target_depth) / sample_depth`. Passing `depths: None` (no
/// `--depths` file given) leaves the count unnormalized, just cast to
/// `f64`.
fn normalize_count(
    count: usize,
    sample_id: &str,
    depths: Option<&HashMap<String, f64>>,
    target_depth: f64,
) -> Result<f64> {
    match depths {
        Some(depths) => {
            let depth = depths
                .get(sample_id)
                .with_context(|| format!("sample {sample_id:?} has no entry in the depths file"))?;
            Ok(count as f64 * target_depth / depth)
        }
        None => Ok(count as f64),
    }
}

fn locus_to_output(
    locus: MergedLocus,
    depths: Option<&HashMap<String, f64>>,
    target_depth: f64,
) -> Result<MergedLocusOutput> {
    let mut irr_count = BTreeMap::new();
    let mut motifs: BTreeMap<String, BTreeMap<String, f64>> = BTreeMap::new();
    for (sample_id, summary) in locus.samples {
        irr_count.insert(
            sample_id.clone(),
            normalize_count(summary.irr_count, &sample_id, depths, target_depth)?,
        );
        for (motif, count) in summary.motifs {
            let normalized = normalize_count(count, &sample_id, depths, target_depth)?;
            motifs.entry(motif).or_default().insert(sample_id.clone(), normalized);
        }
    }
    Ok(MergedLocusOutput {
        chrom: locus.chrom,
        start: locus.start,
        end: locus.end,
        irr_count,
        motifs,
    })
}

/// A region plus everything it contributes to whichever merged locus it
/// ends up clustered into: either a previous merge iteration's already
/// merged locus (folded in as a whole) or one sample's single input locus
/// (folded in as a single-sample contribution).
struct Contribution {
    region: Region,
    samples: BTreeMap<String, SampleLocusSummary>,
}

fn parse_manifest(path: &Path) -> Result<Vec<(String, PathBuf)>> {
    let file = File::open(path).with_context(|| format!("failed to open manifest {path:?}"))?;
    let reader = BufReader::new(file);

    let mut entries = Vec::new();
    for (line_no, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("failed to read {path:?} at line {}", line_no + 1))?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        let mut fields = line.splitn(2, '\t');
        let sample_id = fields
            .next()
            .with_context(|| format!("{path:?}:{}: missing sample_id column", line_no + 1))?;
        let summary_path = fields
            .next()
            .with_context(|| format!("{path:?}:{}: missing summary path column", line_no + 1))?;

        entries.push((sample_id.to_string(), PathBuf::from(summary_path)));
    }

    Ok(entries)
}

fn parse_depths(path: &Path) -> Result<HashMap<String, f64>> {
    let file = File::open(path).with_context(|| format!("failed to open depths file {path:?}"))?;
    let reader = BufReader::new(file);

    let mut depths = HashMap::new();
    for (line_no, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("failed to read {path:?} at line {}", line_no + 1))?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        let mut fields = line.splitn(2, '\t');
        let sample_id = fields
            .next()
            .with_context(|| format!("{path:?}:{}: missing sample_id column", line_no + 1))?;
        let depth: f64 = fields
            .next()
            .with_context(|| format!("{path:?}:{}: missing depth column", line_no + 1))?
            .trim()
            .parse()
            .with_context(|| format!("{path:?}:{}: invalid depth value", line_no + 1))?;

        depths.insert(sample_id.to_string(), depth);
    }

    Ok(depths)
}

pub fn run(args: MergeArgs) -> Result<()> {
    if args.batch_size == 0 {
        bail!("--batch-size must be greater than zero");
    }

    let manifest_entries = parse_manifest(&args.manifest)?;
    if manifest_entries.is_empty() {
        log::warn!("manifest {:?} contained no samples", args.manifest);
    }

    let depths = match &args.depths {
        Some(path) => Some(parse_depths(path)?),
        None => {
            log::warn!(
                "no --depths file given: IRR counts will be written raw, not normalized to a \
                 target depth"
            );
            None
        }
    };

    // Contig name -> synthetic tid, assigned in order of first appearance
    // and kept stable across batches so merged loci from earlier batches
    // still resolve to the same tid as matching contigs in later ones.
    // `chrom_names` is the inverse, indexed by tid.
    let mut chrom_ids: HashMap<String, i32> = HashMap::new();
    let mut chrom_names: Vec<String> = Vec::new();

    let mut merged: Vec<MergedLocus> = Vec::new();

    let batch_count = manifest_entries.len().div_ceil(args.batch_size);
    for (batch_idx, batch) in manifest_entries.chunks(args.batch_size).enumerate() {
        let mut contributions: Vec<Contribution> = Vec::with_capacity(merged.len() + batch.len());

        // Fold in the running merge from prior batches as whole
        // contributions, so it can be re-clustered together with this
        // batch's new loci.
        for locus in merged.drain(..) {
            let tid = *chrom_ids
                .get(&locus.chrom)
                .expect("chrom of a previously merged locus must already be registered");
            contributions.push(Contribution {
                region: Region {
                    tid,
                    start: locus.start,
                    end: locus.end,
                },
                samples: locus.samples,
            });
        }

        for (sample_id, summary_path) in batch {
            let file = File::open(summary_path).with_context(|| {
                format!("failed to open summary {summary_path:?} for sample {sample_id:?}")
            })?;
            let loci: Vec<AnchorRegionSummaryInput> = serde_json::from_reader(BufReader::new(file))
                .with_context(|| {
                    format!("failed to parse summary {summary_path:?} for sample {sample_id:?}")
                })?;

            for locus in loci {
                let next_id = chrom_names.len() as i32;
                let tid = *chrom_ids.entry(locus.chrom.clone()).or_insert(next_id);
                if tid == next_id {
                    chrom_names.push(locus.chrom);
                }

                let mut samples = BTreeMap::new();
                samples.insert(
                    sample_id.clone(),
                    SampleLocusSummary {
                        irr_count: locus.irr_count,
                        motifs: locus.motifs,
                    },
                );

                contributions.push(Contribution {
                    region: Region {
                        tid,
                        start: locus.start,
                        end: locus.end,
                    },
                    samples,
                });
            }
        }

        let regions: Vec<Region> = contributions.iter().map(|c| c.region).collect();
        let clusters = bed::merge_within(&regions, args.merge_distance);

        let mut next_merged: Vec<MergedLocus> = clusters
            .iter()
            .map(|region| MergedLocus {
                chrom: chrom_names[region.tid as usize].clone(),
                start: region.start,
                end: region.end,
                samples: BTreeMap::new(),
            })
            .collect();

        for contribution in contributions {
            let cluster_idx = bed::locate(&clusters, contribution.region.tid, contribution.region.start)
                .expect("every contribution's region must fall within its own merged cluster");
            let target_samples = &mut next_merged[cluster_idx].samples;

            for (sample_id, sample_summary) in contribution.samples {
                // A sample can only land twice in the same cluster if two
                // of its own loci that were still distinct going into this
                // batch get bridged together now, e.g. via another
                // sample's locus landing between them within
                // `--merge-distance`.
                match target_samples.entry(sample_id) {
                    std::collections::btree_map::Entry::Vacant(entry) => {
                        entry.insert(sample_summary);
                    }
                    std::collections::btree_map::Entry::Occupied(mut entry) => {
                        let existing = entry.get_mut();
                        existing.irr_count += sample_summary.irr_count;
                        for (motif, count) in sample_summary.motifs {
                            *existing.motifs.entry(motif).or_insert(0) += count;
                        }
                    }
                }
            }
        }

        merged = next_merged;

        log::info!(
            "merge: batch {}/{batch_count} ({} samples) folded in, {} merged loci so far",
            batch_idx + 1,
            batch.len(),
            merged.len(),
        );
    }

    let merged_count = merged.len();
    let output: Vec<MergedLocusOutput> = merged
        .into_iter()
        .map(|locus| locus_to_output(locus, depths.as_ref(), args.target_depth))
        .collect::<Result<_>>()?;

    let json =
        serde_json::to_string_pretty(&output).context("failed to serialize merged locus summary")?;
    std::fs::write(&args.output, json)
        .with_context(|| format!("failed to write merged summary {:?}", args.output))?;

    log::info!(
        "merge: {} samples from {:?} merged into {merged_count} loci, written to {:?}",
        manifest_entries.len(),
        args.manifest,
        args.output,
    );

    Ok(())
}
