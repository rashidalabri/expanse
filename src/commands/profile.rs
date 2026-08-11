use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use clap::{Args, ValueEnum};
use rust_htslib::bam::{self, Format, Header, IndexedReader, Read as BamRead, Record};
use serde::Serialize;
use url::Url;

use crate::bed::{self, Region};
use crate::crai;
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

    /// Output CRAM/BAM path.
    #[arg(short = 'o', long)]
    pub output: PathBuf,

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

    /// Minimum MAPQ required of a mate (anchor read) for it to be retrieved.
    #[arg(long, default_value_t = 50)]
    pub min_anchor_mapq: u8,

    /// Merge IRR anchor (mate) locations within this many bp of each other
    /// into one anchor region before applying --min-irrs-per-anchor.
    #[arg(long, default_value_t = 1000)]
    pub anchor_merge_distance: i64,

    /// Minimum number of distinct IRR candidates a merged anchor region
    /// must have to be retrieved in pass 2. Anchor regions below this
    /// threshold are discarded, along with their associated IRR candidates.
    #[arg(long, default_value_t = 2)]
    pub min_irrs_per_anchor: usize,

    /// Write every pass-1 IRR candidate as-is, skipping mate retrieval
    /// (pass 2) and the passing-anchor requirement entirely.
    #[arg(long)]
    pub irr_only: bool,

    /// Write a JSON summary to this path: for each surviving merged anchor
    /// region, its coordinates (0-based, half-open) and a breakdown of IRR
    /// counts by motif. Has no effect with --irr-only, since anchor
    /// clustering doesn't run in that mode.
    #[arg(long)]
    pub anchor_motif_summary: Option<PathBuf>,

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

/// One entry in the `--anchor-motif-summary` JSON output: a merged anchor
/// region that survived `--min-irrs-per-anchor`, with its IRR support
/// broken down by (canonical) motif.
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
    let resolved_output_format = resolve_output_format(&args.output, args.output_format);

    if (input_is_cram || resolved_output_format == Format::Cram) && args.reference.is_none() {
        bail!(
            "--reference is required when CRAM is involved (input={}, output={})",
            args.input,
            args.output.display()
        );
    }

    if args.irr_only && args.anchor_motif_summary.is_some() {
        log::warn!(
            "--anchor-motif-summary has no effect with --irr-only (anchor clustering is skipped); no summary will be written"
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

    // For CRAM inputs, group regions by the underlying CRAM slice they land
    // in and merge each group down to the bounding box of its own member
    // regions (not the slice's own, generally wider, span), so htslib only
    // seeks/decompresses each slice once rather than once per (possibly much
    // narrower) BED/mate region. A group can still cover a gap between its
    // members that no original region touches, so callers must re-check
    // membership against the pre-merge region set (via `bed::contains`)
    // while iterating.
    let crai_slices = if input_is_cram {
        match crai::load(&args.input) {
            Ok(slices) => Some(slices),
            Err(err) => {
                log::warn!(
                    "failed to load CRAI slice index for {}, falling back to per-region fetch: {err:#}",
                    args.input
                );
                None
            }
        }
    } else {
        None
    };

    let bed_regions_merged = bed::merge_regions(&bed_regions);
    let bed_fetch_regions = match &crai_slices {
        Some(slices) => bed::merge_regions(&crai::merge_by_slice(&bed_regions, slices)),
        None => bed_regions.clone(),
    };

    let mut writer = bam::Writer::from_path(&args.output, &header, resolved_output_format)
        .with_context(|| format!("failed to open output {:?}", args.output))?;
    if resolved_output_format == Format::Cram
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

    // Pass-1 IRR candidates are buffered rather than written immediately:
    // whether each one is kept depends on whether its mate turns out to
    // pass the anchor filter in pass 2 (see the final filtering step below).
    let mut candidate_records: Vec<Record> = Vec::new();
    // Each candidate's canonical motif, parallel to `candidate_records`
    // (same length, same order) -- used to break IRR counts down by motif
    // for `--anchor-motif-summary`.
    let mut candidate_motifs: Vec<Vec<u8>> = Vec::new();
    let mut candidate_keys: HashSet<(i32, i64, Vec<u8>, u16)> = HashSet::new();

    let mut pass1_count: u64 = 0;
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

            if crai_slices.is_some()
                && !bed::contains(&bed_regions_merged, record.tid(), record.pos())
            {
                continue;
            }

            if record.is_unmapped()
                || record.is_secondary()
                || record.is_supplementary()
                || record.is_mate_unmapped()
            {
                continue;
            }
            if record.mapq() >= args.max_irr_mapq {
                continue;
            }
            let Some(motif) = irr::classify_in_repeat_read(
                &record.seq().as_bytes(),
                record.qual(),
                args.motif_min_len,
                args.motif_max_len,
            ) else {
                continue;
            };

            pass1_count += 1;

            let key = (
                record.tid(),
                record.pos(),
                record.qname().to_vec(),
                record.flags(),
            );
            if candidate_keys.insert(key) {
                candidate_records.push(record.clone());
                candidate_motifs.push(motif);
            }
        }
    }

    let mut written_keys: HashSet<(i32, i64, Vec<u8>, u16)> = HashSet::new();
    // Mate identities (qname, is_first_in_template) that were actually
    // found and passed --min-anchor-mapq, i.e. anchors the candidates above
    // are allowed to keep. Left empty (and never consulted, since
    // `args.irr_only` short-circuits the check below) when pass 2 is
    // skipped entirely.
    let mut satisfied_mates: HashSet<(Vec<u8>, bool)> = HashSet::new();

    let mut anchor_cluster_count = 0;
    let mut anchor_cluster_dropped_count = 0;
    let mut anchor_dropped_irr_count: u64 = 0;
    let mut merged_mate_region_count = 0;
    let mut mate_fetch_region_count = 0;
    let mut pass2_count: u64 = 0;

    if !args.irr_only {
        // Cluster each candidate's mate location by proximity
        // (--anchor-merge-distance) rather than just overlap/touching, then
        // require --min-irrs-per-anchor distinct IRR candidates per cluster.
        // Clusters (and their candidates) below that threshold are dropped
        // before pass 2 ever looks for them.
        let mate_targets: Vec<Region> = candidate_records
            .iter()
            .map(|candidate| Region {
                tid: candidate.mtid(),
                start: candidate.mpos(),
                end: candidate.mpos() + 1,
            })
            .collect();

        let anchor_clusters = bed::merge_within(&mate_targets, args.anchor_merge_distance);
        anchor_cluster_count = anchor_clusters.len();

        let candidate_cluster: Vec<usize> = mate_targets
            .iter()
            .map(|target| {
                bed::locate(&anchor_clusters, target.tid, target.start)
                    .expect("every mate target must fall within its own merged cluster")
            })
            .collect();
        let mut cluster_irr_counts = vec![0usize; anchor_clusters.len()];
        let mut cluster_motif_counts: Vec<HashMap<Vec<u8>, usize>> =
            vec![HashMap::new(); anchor_clusters.len()];
        for (&cluster_idx, motif) in candidate_cluster.iter().zip(candidate_motifs.iter()) {
            cluster_irr_counts[cluster_idx] += 1;
            *cluster_motif_counts[cluster_idx]
                .entry(motif.clone())
                .or_insert(0) += 1;
        }

        let surviving_cluster_indices: Vec<usize> = (0..anchor_clusters.len())
            .filter(|&idx| cluster_irr_counts[idx] >= args.min_irrs_per_anchor)
            .collect();
        let merged_mate_regions: Vec<Region> = surviving_cluster_indices
            .iter()
            .map(|&idx| anchor_clusters[idx])
            .collect();
        anchor_cluster_dropped_count = anchor_cluster_count - merged_mate_regions.len();

        if let Some(summary_path) = &args.anchor_motif_summary {
            let summaries: Vec<AnchorRegionSummary> = surviving_cluster_indices
                .iter()
                .map(|&idx| {
                    let region = anchor_clusters[idx];
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
                .context("failed to serialize anchor motif summary")?;
            std::fs::write(summary_path, json).with_context(|| {
                format!("failed to write anchor motif summary {summary_path:?}")
            })?;
        }

        // (qname, is_first_in_template) identifying the specific mate we
        // still need; only populated for candidates whose anchor cluster
        // survived the --min-irrs-per-anchor threshold above.
        let mut wanted_mates: HashSet<(Vec<u8>, bool)> = HashSet::new();
        let mut surviving_candidates = Vec::with_capacity(candidate_records.len());
        for (candidate, cluster_idx) in candidate_records.into_iter().zip(candidate_cluster) {
            if cluster_irr_counts[cluster_idx] >= args.min_irrs_per_anchor {
                wanted_mates.insert((
                    candidate.qname().to_vec(),
                    !candidate.is_first_in_template(),
                ));
                surviving_candidates.push(candidate);
            } else {
                anchor_dropped_irr_count += 1;
            }
        }
        candidate_records = surviving_candidates;

        let mate_fetch_regions = match &crai_slices {
            Some(slices) => bed::merge_regions(&crai::merge_by_slice(&merged_mate_regions, slices)),
            None => merged_mate_regions.clone(),
        };
        merged_mate_region_count = merged_mate_regions.len();
        mate_fetch_region_count = mate_fetch_regions.len();

        for region in &mate_fetch_regions {
            reader
                .fetch((region.tid, region.start, region.end))
                .with_context(|| {
                    format!(
                        "failed to seek to mate region tid={} {}-{}",
                        region.tid, region.start, region.end
                    )
                })?;

            while let Some(result) = reader.read(&mut record) {
                result.context("failed to read record in pass 2")?;

                if crai_slices.is_some()
                    && !bed::contains(&merged_mate_regions, record.tid(), record.pos())
                {
                    continue;
                }

                if record.is_secondary() || record.is_supplementary() {
                    continue;
                }
                if record.mapq() < args.min_anchor_mapq {
                    continue;
                }

                let mate_key = (record.qname().to_vec(), record.is_first_in_template());
                if !wanted_mates.contains(&mate_key) {
                    continue;
                }

                satisfied_mates.insert(mate_key);

                let key = (
                    record.tid(),
                    record.pos(),
                    record.qname().to_vec(),
                    record.flags(),
                );
                if written_keys.insert(key) {
                    writer
                        .write(&record)
                        .context("failed to write pass-2 mate record")?;
                    pass2_count += 1;
                }
            }
        }
    }

    // Final step: in the default two-pass mode, drop IRR candidates whose
    // mate was mapped but never made it into satisfied_mates (not found in
    // its merged region, or found but below --min-anchor-mapq); candidates
    // with an unmapped mate never had an anchor to look for, so the filter
    // doesn't apply to them. In --irr-only mode every candidate is written
    // as-is, since pass 2 never ran and satisfied_mates is empty.
    let mut pass1_written_count: u64 = 0;
    let mut pass1_dropped_count: u64 = 0;
    for candidate in &candidate_records {
        if !args.irr_only && !candidate.is_mate_unmapped() {
            let expected_mate_key = (
                candidate.qname().to_vec(),
                !candidate.is_first_in_template(),
            );
            if !satisfied_mates.contains(&expected_mate_key) {
                pass1_dropped_count += 1;
                continue;
            }
        }

        let key = (
            candidate.tid(),
            candidate.pos(),
            candidate.qname().to_vec(),
            candidate.flags(),
        );
        if written_keys.insert(key) {
            writer
                .write(candidate)
                .context("failed to write pass-1 record")?;
            pass1_written_count += 1;
        }
    }

    if args.irr_only {
        log::info!(
            "profile: {} region fetches ({} BED regions), {pass1_count} candidate IRR reads, \
             {pass1_written_count} IRR reads written (--irr-only, pass 2 skipped)",
            bed_fetch_regions.len(),
            bed_regions.len(),
        );
    } else {
        log::info!(
            "profile: {} region fetches ({} BED regions), {pass1_count} candidate IRR reads, \
             {anchor_cluster_count} anchor clusters ({anchor_cluster_dropped_count} dropped for having \
             fewer than {} IRR reads, taking {anchor_dropped_irr_count} IRR reads with them), \
             {merged_mate_region_count} surviving anchor regions ({mate_fetch_region_count} mate fetches), \
             {pass2_count} mate reads written, \
             {pass1_dropped_count} IRR reads dropped for lacking a passing anchor, {pass1_written_count} IRR reads written, \
             {} total records written",
            bed_fetch_regions.len(),
            bed_regions.len(),
            args.min_irrs_per_anchor,
            written_keys.len(),
        );
    }

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
