use std::collections::HashSet;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use clap::{Args, ValueEnum};
use rust_htslib::bam::{self, Format, Header, IndexedReader, Read as BamRead, Record};
use url::Url;

use crate::bed::{self, Region};
use crate::crai;
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

    /// Write every pass-1 IRR candidate as-is, skipping mate retrieval
    /// (pass 2) and the passing-anchor requirement entirely.
    #[arg(long)]
    pub irr_only: bool,

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
    // in so htslib only seeks/decompresses each slice once, rather than once
    // per (possibly much narrower) BED/mate region. This can pull in extra
    // reads that fall inside a shared slice but outside any originally
    // requested region, so callers must re-check membership against the
    // pre-expansion region set (via `bed::contains`) while iterating.
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
        Some(slices) => bed::merge_regions(&crai::expand_to_slices(&bed_regions, slices)),
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
    let mut candidate_keys: HashSet<(i32, i64, Vec<u8>, u16)> = HashSet::new();
    let mut mate_targets: Vec<Region> = Vec::new();
    // (qname, is_first_in_template) identifying the specific mate we still need.
    let mut wanted_mates: HashSet<(Vec<u8>, bool)> = HashSet::new();

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
            if irr::classify_in_repeat_read(
                &record.seq().as_bytes(),
                record.qual(),
                args.motif_min_len,
                args.motif_max_len,
            )
            .is_none()
            {
                continue;
            }

            pass1_count += 1;

            let key = (
                record.tid(),
                record.pos(),
                record.qname().to_vec(),
                record.flags(),
            );
            if candidate_keys.insert(key) {
                candidate_records.push(record.clone());
            }

            if !args.irr_only {
                mate_targets.push(Region {
                    tid: record.mtid(),
                    start: record.mpos(),
                    end: record.mpos() + 1,
                });
                wanted_mates.insert((record.qname().to_vec(), !record.is_first_in_template()));
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

    let mut merged_mate_region_count = 0;
    let mut mate_fetch_region_count = 0;
    let mut pass2_count: u64 = 0;

    if !args.irr_only {
        let merged_mate_regions = bed::merge_regions(&mate_targets);
        let mate_fetch_regions = match &crai_slices {
            Some(slices) => {
                bed::merge_regions(&crai::expand_to_slices(&merged_mate_regions, slices))
            }
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
             {merged_mate_region_count} merged mate regions ({mate_fetch_region_count} mate fetches), \
             {pass2_count} mate reads written, \
             {pass1_dropped_count} IRR reads dropped for lacking a passing anchor, {pass1_written_count} IRR reads written, \
             {} total records written",
            bed_fetch_regions.len(),
            bed_regions.len(),
            written_keys.len(),
        );
    }

    Ok(())
}

fn is_cram_path(path: &str) -> bool {
    let path_no_query = path.split(['?', '#']).next().unwrap_or(path);
    path_no_query.to_ascii_lowercase().ends_with(".cram")
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
    fn is_cram_path_detects_extension_case_insensitively() {
        assert!(is_cram_path("foo.cram"));
        assert!(is_cram_path("foo.CRAM"));
        assert!(is_cram_path("s3://bucket/foo.cram"));
        assert!(is_cram_path("s3://bucket/foo.cram?X-Amz-Signature=abc"));
        assert!(!is_cram_path("foo.bam"));
    }

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
