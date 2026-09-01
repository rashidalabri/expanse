use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use clap::{Args, ValueEnum};
use rust_htslib::bam::{self, Format, Header, IndexedReader, Read as BamRead, Record};
use serde::Serialize;
use url::Url;

use crate::bed::{self, Region};
use crate::hts_io::is_cram_path;
use crate::irr;

#[derive(Args, Debug)]
pub struct ProfileArgs {
    /// BED file of IRR-mapping regions to scan for candidate low-mapq reads.
    #[arg(long)]
    pub bed: PathBuf,

    /// Input CRAM/BAM: local path or s3:// / gs:// / https:// URL.
    #[arg(short = 'i', long)]
    pub input: String,

    /// JSON summary output path: for each merged anchor region, its
    /// coordinates (0-based, half-open) and a breakdown of IRR counts by
    /// motif. Always written.
    #[arg(long)]
    pub summary: PathBuf,

    /// Optional output CRAM/BAM path for the candidate IRR reads themselves.
    /// No alignment file is written unless this is given.
    #[arg(short = 'o', long)]
    pub output: Option<PathBuf>,

    /// Reads overlapping an IRR region with MAPQ below this are candidates.
    #[arg(long, default_value_t = 40)]
    pub max_irr_mapq: u8,

    /// Shortest repeat-unit (motif) length considered when checking whether
    /// a candidate read is an in-repeat read (IRR).
    #[arg(long, default_value_t = irr::DEFAULT_MOTIF_MIN_LEN)]
    pub motif_min_len: u32,

    /// Longest repeat-unit (motif) length considered.
    #[arg(long, default_value_t = irr::DEFAULT_MOTIF_MAX_LEN)]
    pub motif_max_len: u32,

    /// Merge IRR candidate locations within this many bp of each other into
    /// one anchor region for `--summary`.
    #[arg(long, default_value_t = 500)]
    pub anchor_merge_distance: i64,

    /// Reference FASTA. Required when the input or output uses CRAM.
    #[arg(short = 'r', long)]
    pub reference: Option<PathBuf>,

    /// Override output format inference from the --output extension.
    #[arg(long, value_enum)]
    pub output_format: Option<OutputFormat>,

    /// Number of htslib I/O threads to use for reading and writing.
    #[arg(long, default_value_t = 1)]
    pub threads: usize,
}

#[derive(ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
pub enum OutputFormat {
    Bam,
    Cram,
}

/// One entry in the `--summary` JSON output: a merged anchor region (IRR
/// candidate locations within `--anchor-merge-distance` bp of each other),
/// with its IRR support broken down by (canonical) motif. A read that
/// qualifies under more than one motif is counted once in `irr_count` but
/// once per motif in `motifs`, so the `motifs` values can sum to more than
/// `irr_count`. A motif key may contain IUPAC ambiguity codes (e.g. `GCN`,
/// `AARRG`) at positions that are consistently mixed across repeat copies;
/// see `irr::classify_in_repeat_read_all`.
#[derive(Serialize, Debug)]
struct AnchorRegionSummary {
    chrom: String,
    /// 0-based, half-open, matching this crate's internal `Region` convention.
    start: i64,
    end: i64,
    irr_count: usize,
    motifs: BTreeMap<String, usize>,
}

pub fn run(args: ProfileArgs) -> Result<()> {
    let input_is_cram = is_cram_path(&args.input);
    let resolved_output_format = args
        .output
        .as_ref()
        .map(|output| resolve_output_format(output, args.output_format));

    if (input_is_cram || resolved_output_format == Some(Format::Cram)) && args.reference.is_none() {
        bail!(
            "--reference is required when CRAM is involved (input={}, output={:?})",
            args.input,
            args.output
        );
    }

    let mut reader = match Url::parse(&args.input) {
        Ok(url) => IndexedReader::from_url(&url)
            .with_context(|| format!("failed to open indexed input {}", args.input))?,
        Err(_) => IndexedReader::from_path(&args.input)
            .with_context(|| format!("failed to open indexed input {}", args.input))?,
    };
    if let Some(reference) = &args.reference {
        reader
            .set_reference(reference)
            .with_context(|| format!("failed to set reference {reference:?} on reader"))?;
    }
    if args.threads > 1 {
        reader
            .set_threads(args.threads)
            .context("failed to set reader thread count")?;
    }

    let header = Header::from_template(reader.header());

    let bed_regions = bed::parse_bed(&args.bed, reader.header())
        .with_context(|| format!("failed to parse BED file {:?}", args.bed))?;
    if bed_regions.is_empty() {
        log::warn!("BED file {:?} contained no usable regions", args.bed);
    }

    let bed_fetch_regions = bed::merge_regions(&bed_regions);

    let mut writer = match (&args.output, resolved_output_format) {
        (Some(output_path), Some(format)) => {
            let mut writer = bam::Writer::from_path(output_path, &header, format)
                .with_context(|| format!("failed to open output {output_path:?}"))?;
            if format == Format::Cram
                && let Some(reference) = &args.reference
            {
                writer
                    .set_reference(reference)
                    .with_context(|| format!("failed to set reference {reference:?} on writer"))?;
            }
            if args.threads > 1 {
                writer
                    .set_threads(args.threads)
                    .context("failed to set writer thread count")?;
            }
            Some(writer)
        }
        _ => None,
    };

    let mut candidate_keys: HashSet<(i32, i64, Vec<u8>, u16)> = HashSet::new();
    // Each surviving candidate's own location and every canonical motif it
    // qualifies under, used to build the merged anchor-region summary below.
    let mut candidate_regions: Vec<Region> = Vec::new();
    let mut candidate_motifs: Vec<Vec<Vec<u8>>> = Vec::new();

    let mut pass_count: u64 = 0;
    let mut written_count: u64 = 0;
    let mut record = Record::new();
    for region in &bed_fetch_regions {
        reader
            .fetch((region.tid, region.start, region.end))
            .with_context(|| {
                format!(
                    "failed to seek to region tid={} {}-{}",
                    region.tid, region.start, region.end
                )
            })?;

        while let Some(result) = reader.read(&mut record) {
            result.context("failed to read record in pass 1")?;

            if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
                continue;
            }
            if record.mapq() > args.max_irr_mapq {
                continue;
            }
            let motifs = irr::classify_in_repeat_read_all(
                &record.seq().as_bytes(),
                record.qual(),
                args.motif_min_len,
                args.motif_max_len,
            );
            if motifs.is_empty() {
                continue;
            }

            pass_count += 1;

            let key = (
                record.tid(),
                record.pos(),
                record.qname().to_vec(),
                record.flags(),
            );
            if candidate_keys.insert(key) {
                candidate_regions.push(Region {
                    tid: record.tid(),
                    start: record.pos(),
                    end: record.pos() + 1,
                });
                candidate_motifs.push(motifs);

                if let Some(writer) = writer.as_mut() {
                    writer
                        .write(&record)
                        .context("failed to write pass-1 record")?;
                    written_count += 1;
                }
            }
        }
    }

    // Only needed for the dedup check above; free it before the clustering
    // step below allocates its own structures.
    drop(candidate_keys);

    let anchor_clusters = bed::merge_within(&candidate_regions, args.anchor_merge_distance);
    let candidate_cluster: Vec<usize> = candidate_regions
        .iter()
        .map(|region| {
            bed::locate(&anchor_clusters, region.tid, region.start)
                .expect("every candidate must fall within its own merged cluster")
        })
        .collect();
    // Save the count for the log message below, then free the rest: no
    // longer needed now that every candidate has been assigned a cluster.
    let candidate_count = candidate_regions.len();
    drop(candidate_regions);

    // `cluster_irr_counts` counts each candidate read once per cluster,
    // regardless of how many motifs it qualified under; `cluster_motif_counts`
    // counts it once per qualifying motif, so a multi-motif read is
    // reflected in more than one motif's tally without inflating the
    // region's overall IRR total.
    let mut cluster_irr_counts = vec![0usize; anchor_clusters.len()];
    let mut cluster_motif_counts: Vec<HashMap<Vec<u8>, usize>> =
        vec![HashMap::new(); anchor_clusters.len()];
    for (&cluster_idx, motifs) in candidate_cluster.iter().zip(candidate_motifs.into_iter()) {
        cluster_irr_counts[cluster_idx] += 1;
        for motif in motifs {
            *cluster_motif_counts[cluster_idx].entry(motif).or_insert(0) += 1;
        }
    }

    let summaries: Vec<AnchorRegionSummary> = anchor_clusters
        .iter()
        .enumerate()
        .map(|(idx, region)| {
            let motifs: BTreeMap<String, usize> = cluster_motif_counts[idx]
                .iter()
                .map(|(motif, &count)| (String::from_utf8_lossy(motif).into_owned(), count))
                .collect();
            AnchorRegionSummary {
                chrom: String::from_utf8_lossy(reader.header().tid2name(region.tid as u32))
                    .into_owned(),
                start: region.start,
                end: region.end,
                irr_count: cluster_irr_counts[idx],
                motifs,
            }
        })
        .collect();

    let json = serde_json::to_string_pretty(&summaries)
        .context("failed to serialize anchor region summary")?;
    std::fs::write(&args.summary, json)
        .with_context(|| format!("failed to write anchor region summary {:?}", args.summary))?;

    log::info!(
        "profile: {} region fetches ({} BED regions), {pass_count} candidate IRR reads \
         ({} distinct), {} anchor regions summarized to {:?}{}",
        bed_fetch_regions.len(),
        bed_regions.len(),
        candidate_count,
        anchor_clusters.len(),
        args.summary,
        if args.output.is_some() {
            format!(", {written_count} IRR reads written to {:?}", args.output)
        } else {
            String::new()
        },
    );

    Ok(())
}

fn resolve_output_format(output: &Path, override_format: Option<OutputFormat>) -> Format {
    match override_format {
        Some(OutputFormat::Bam) => Format::Bam,
        Some(OutputFormat::Cram) => Format::Cram,
        None => {
            let ext = output
                .extension()
                .and_then(|e| e.to_str())
                .unwrap_or("")
                .to_ascii_lowercase();
            if ext == "cram" {
                Format::Cram
            } else {
                Format::Bam
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resolve_output_format_infers_from_extension() {
        assert_eq!(
            resolve_output_format(Path::new("out.bam"), None),
            Format::Bam
        );
        assert_eq!(
            resolve_output_format(Path::new("out.cram"), None),
            Format::Cram
        );
        assert_eq!(
            resolve_output_format(Path::new("out.CRAM"), None),
            Format::Cram
        );
        assert_eq!(resolve_output_format(Path::new("out"), None), Format::Bam);
    }

    #[test]
    fn resolve_output_format_override_wins() {
        assert_eq!(
            resolve_output_format(Path::new("out.bam"), Some(OutputFormat::Cram)),
            Format::Cram
        );
        assert_eq!(
            resolve_output_format(Path::new("out.cram"), Some(OutputFormat::Bam)),
            Format::Bam
        );
    }
}
