//! Scans an entire CRAM/BAM (no BED-restricted fetch, unlike `profile`) for
//! in-repeat reads (IRRs) -- low-MAPQ reads with at least one qualifying
//! motif -- and reports where they cluster, as a BED file: each row is a
//! same-motif group of IRR-read alignment spans merged within
//! `--merge-distance` bp of each other, with that motif and the number of
//! IRR reads contributing to it. Regions are merged separately per motif,
//! so the same coordinates can appear on more than one row if reads there
//! qualify under more than one motif.
//!
//! Intended to build a `--sink-bed` / `--exclude-bed` input for `profile`:
//! loci that are themselves saturated with IRRs are exactly the ones whose
//! anchor-mate evidence isn't informative and should be excluded from its
//! `--summary`.

use std::collections::HashMap;
use std::io::Write as IoWrite;
use std::path::PathBuf;

use anyhow::{Context, Result, bail};
use clap::Args;
use rust_htslib::bam::{self, Read as BamRead, Record};
use url::Url;

use crate::bed::{self, Region};
use crate::hts_io::is_cram_path;
use crate::irr;

#[derive(Args, Debug)]
pub struct SinksArgs {
    /// Input CRAM/BAM: local path or s3:// / gs:// / https:// URL. Scanned
    /// in its entirety; no index is required.
    #[arg(short = 'i', long)]
    pub input: String,

    /// BED output path: one row per merged same-motif group of IRR-read
    /// regions, as `chrom<TAB>start<TAB>end<TAB>motif<TAB>irr_count`
    /// (0-based, half-open coordinates).
    #[arg(short = 'o', long)]
    pub output: PathBuf,

    /// Reads with MAPQ below this are candidates for IRR classification.
    #[arg(long, default_value_t = 40)]
    pub max_irr_mapq: u8,

    /// Shortest repeat-unit (motif) length considered when checking whether
    /// a candidate read is an in-repeat read (IRR).
    #[arg(long, default_value_t = irr::DEFAULT_MOTIF_MIN_LEN)]
    pub motif_min_len: u32,

    /// Longest repeat-unit (motif) length considered.
    #[arg(long, default_value_t = irr::DEFAULT_MOTIF_MAX_LEN)]
    pub motif_max_len: u32,

    /// Maximum number of IUPAC-ambiguous (non-A/C/G/T) positions allowed in
    /// a mononucleotide (1bp) motif; motifs exceeding this are rejected.
    #[arg(long, default_value_t = irr::DEFAULT_MAX_DEGENERATE_MONONUCLEOTIDE)]
    pub max_degenerate_mononucleotide: u32,

    /// Same, for a dinucleotide (2bp) motif.
    #[arg(long, default_value_t = irr::DEFAULT_MAX_DEGENERATE_DINUCLEOTIDE)]
    pub max_degenerate_dinucleotide: u32,

    /// Same, for a trinucleotide (3bp) motif.
    #[arg(long, default_value_t = irr::DEFAULT_MAX_DEGENERATE_TRINUCLEOTIDE)]
    pub max_degenerate_trinucleotide: u32,

    /// Same, for any motif of 4bp or longer.
    #[arg(long, default_value_t = irr::DEFAULT_MAX_DEGENERATE_OTHER)]
    pub max_degenerate_other: u32,

    /// Merge same-motif IRR-read regions within this many bp of each other
    /// into one output region.
    #[arg(long, default_value_t = 0)]
    pub merge_distance: i64,

    /// Reference FASTA. Required when the input uses CRAM.
    #[arg(short = 'r', long)]
    pub reference: Option<PathBuf>,

    /// Number of htslib I/O threads to use for reading.
    #[arg(long, default_value_t = 1)]
    pub threads: usize,
}

pub fn run(args: SinksArgs) -> Result<()> {
    if is_cram_path(&args.input) && args.reference.is_none() {
        bail!("--reference is required when the input is CRAM (input={})", args.input);
    }

    let mut reader = match Url::parse(&args.input) {
        Ok(url) => bam::Reader::from_url(&url)
            .with_context(|| format!("failed to open input {}", args.input))?,
        Err(_) => bam::Reader::from_path(&args.input)
            .with_context(|| format!("failed to open input {}", args.input))?,
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

    let degenerate_limits = irr::DegenerateLimits {
        mononucleotide: args.max_degenerate_mononucleotide,
        dinucleotide: args.max_degenerate_dinucleotide,
        trinucleotide: args.max_degenerate_trinucleotide,
        other: args.max_degenerate_other,
    };

    // Each IRR read's alignment span is filed under every motif it
    // qualifies under, so regions are only ever merged with other regions
    // of the *same* motif.
    let mut regions_by_motif: HashMap<Vec<u8>, Vec<Region>> = HashMap::new();
    let mut scanned_count: u64 = 0;
    let mut irr_count: u64 = 0;
    let mut record = Record::new();
    while let Some(result) = reader.read(&mut record) {
        result.context("failed to read record")?;
        scanned_count += 1;

        if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
            continue;
        }
        if record.mapq() > args.max_irr_mapq {
            continue;
        }
        let motifs = irr::identify_repeat_motifs(
            &record.seq().as_bytes(),
            record.qual(),
            args.motif_min_len,
            args.motif_max_len,
            degenerate_limits,
        );
        if motifs.is_empty() {
            continue;
        }
        irr_count += 1;

        let region = Region {
            tid: record.tid(),
            start: record.pos(),
            end: record.cigar().end_pos(),
        };
        for motif in motifs {
            regions_by_motif.entry(motif).or_default().push(region);
        }
    }

    // Merged (region, motif, irr_count) rows, one per same-motif cluster;
    // sorted below for deterministic BED output.
    let mut output_rows: Vec<(Region, Vec<u8>, usize)> = Vec::new();
    for (motif, regions) in &regions_by_motif {
        let clusters = bed::merge_within(regions, args.merge_distance);
        let mut counts = vec![0usize; clusters.len()];
        for region in regions {
            let idx = bed::locate(&clusters, region.tid, region.start)
                .expect("every region must fall within its own merged cluster");
            counts[idx] += 1;
        }
        output_rows.extend(
            clusters
                .into_iter()
                .zip(counts)
                .map(|(cluster, count)| (cluster, motif.clone(), count)),
        );
    }
    output_rows.sort_by(|(a_region, a_motif, _), (b_region, b_motif, _)| {
        (a_region.tid, a_region.start, a_region.end, a_motif).cmp(&(
            b_region.tid,
            b_region.start,
            b_region.end,
            b_motif,
        ))
    });

    let mut writer = std::io::BufWriter::new(
        std::fs::File::create(&args.output)
            .with_context(|| format!("failed to create output BED {:?}", args.output))?,
    );
    for (region, motif, count) in &output_rows {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}",
            String::from_utf8_lossy(reader.header().tid2name(region.tid as u32)),
            region.start,
            region.end,
            String::from_utf8_lossy(motif),
            count
        )
        .context("failed to write output BED row")?;
    }
    writer
        .flush()
        .with_context(|| format!("failed to flush output BED {:?}", args.output))?;

    log::info!(
        "sinks: {scanned_count} records scanned, {irr_count} IRR reads found, {} merged regions \
         written to {:?}",
        output_rows.len(),
        args.output,
    );

    Ok(())
}
