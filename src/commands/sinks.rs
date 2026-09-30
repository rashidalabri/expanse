//! Scans an entire CRAM/BAM (no BED-restricted fetch, unlike `profile`) for
//! in-repeat reads (IRRs) -- low-MAPQ reads with at least one qualifying
//! motif -- and reports where they cluster, as a BED file: each row is a
//! same-motif group of IRR-read alignment spans merged within
//! `--merge-distance` bp of each other, with that motif and the number of
//! IRR reads contributing to it. Each IRR read is attributed to its single
//! best-scoring motif (see `irr::best_repeat_motif`), not every motif it
//! happens to also satisfy, so it contributes to at most one row.
//!
//! An IRR read is only counted if it pairs with a mapped mate whose MAPQ
//! clears `--min-anchor-mapq` (an "anchored IRR pair") or whose mate is
//! itself an IRR ("paired IRRs" -- both reads counted); an IRR read whose
//! mate is neither is dropped. A read is an IRR candidate if it's unmapped
//! *or* its MAPQ clears `--max-irr-mapq` -- reads entirely consumed by a
//! long repeat often fail to map at all, not just map with low confidence
//! -- so a paired IRR's mate may itself be unmapped, as long as it still
//! independently qualifies via its motif content. An anchor, in contrast,
//! must be mapped by definition; an anchored pair whose IRR side happens to
//! be unmapped is still correctly recognized as anchored, but contributes
//! no output row, since `sinks` only ever reports an IRR read's own
//! position and an unmapped read has none.
//!
//! Pairing is resolved in a single pass over the file: every
//! non-secondary/non-supplementary read is classified (IRR candidate,
//! anchor, or neither) and, the first time either mate of a pair is seen,
//! cached by qname until its mate is read; the pair is resolved the moment
//! the second mate arrives, and both cache entries are dropped
//! immediately, so memory stays proportional to how many mate pairs are
//! simultaneously "in flight" in the file, not to the file's total size.
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

    /// Reads with MAPQ at or below this, or unmapped entirely, are
    /// candidates for IRR classification.
    #[arg(long, default_value_t = 40)]
    pub max_irr_mapq: u8,

    /// Minimum MAPQ for a mapped mate to count as an anchor. An IRR
    /// candidate is only counted if its mate clears this MAPQ (an
    /// "anchored IRR pair") or is itself an IRR ("paired IRRs"); an IRR
    /// candidate whose mate is neither is dropped.
    #[arg(long, default_value_t = 50)]
    pub min_anchor_mapq: u8,

    /// Shortest repeat-unit (motif) length considered when checking whether
    /// a candidate read is an in-repeat read (IRR).
    #[arg(long, default_value_t = irr::DEFAULT_MOTIF_MIN_LEN)]
    pub motif_min_len: u32,

    /// Longest repeat-unit (motif) length considered.
    #[arg(long, default_value_t = irr::DEFAULT_MOTIF_MAX_LEN)]
    pub motif_max_len: u32,

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

/// A read's role for IRR/anchor pairing (see [`run`]): either a high-MAPQ
/// anchor, or an IRR candidate with its best-scoring motif and, if mapped,
/// its own alignment span (`None` for an unmapped IRR read, which has no
/// position to report). Cached by qname until its mate is seen, so a pair
/// can be resolved the moment both members have been read once each.
enum Classification {
    Anchor,
    Irr {
        motif: Vec<u8>,
        region: Option<Region>,
    },
}

pub fn run(args: SinksArgs) -> Result<()> {
    if is_cram_path(&args.input) && args.reference.is_none() {
        bail!(
            "--reference is required when the input is CRAM (input={})",
            args.input
        );
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

    // Each counted IRR read's alignment span is filed under its single
    // best-scoring motif only, so regions are only ever merged with other
    // regions of the *same* motif.
    let mut regions_by_motif: HashMap<Vec<u8>, Vec<Region>> = HashMap::new();
    // Reads seen so far whose mate hasn't been reached yet, keyed by
    // qname; resolved (and removed) the moment their mate is read. See the
    // module docs for why this stays small in practice despite having no
    // explicit eviction.
    let mut pending: HashMap<Vec<u8>, Classification> = HashMap::new();
    let mut scanned_count: u64 = 0;
    let mut anchored_irr_count: u64 = 0;
    let mut paired_irr_count: u64 = 0;
    let mut record = Record::new();
    while let Some(result) = reader.read(&mut record) {
        result.context("failed to read record")?;
        scanned_count += 1;

        if record.is_secondary() || record.is_supplementary() {
            continue;
        }

        let mapq = record.mapq();
        // A read entirely consumed by a long repeat often fails to map at
        // all, not just map with low confidence, so an unmapped read is
        // always an IRR candidate regardless of its (otherwise
        // meaningless) MAPQ -- matching the low-MAPQ gate below with `||`.
        let classification = if record.is_unmapped() || mapq <= args.max_irr_mapq {
            let Some(motif) = irr::best_repeat_motif(
                &record.seq().as_bytes(),
                record.qual(),
                args.motif_min_len,
                args.motif_max_len,
            ) else {
                continue;
            };
            let region = (!record.is_unmapped()).then(|| Region {
                tid: record.tid(),
                start: record.pos(),
                end: record.cigar().end_pos(),
            });
            Classification::Irr { motif, region }
        } else if mapq >= args.min_anchor_mapq {
            Classification::Anchor
        } else {
            continue;
        };

        match (pending.remove(record.qname()), classification) {
            (
                Some(Classification::Irr { motif, region }),
                Classification::Irr {
                    motif: mate_motif,
                    region: mate_region,
                },
            ) => {
                // Paired IRRs: each side that has a real alignment
                // position contributes its own region -- an unmapped IRR
                // read has none, but still validly pairs with (and
                // confirms) its mate.
                if let Some(region) = region {
                    regions_by_motif.entry(motif).or_default().push(region);
                    paired_irr_count += 1;
                }
                if let Some(mate_region) = mate_region {
                    regions_by_motif
                        .entry(mate_motif)
                        .or_default()
                        .push(mate_region);
                    paired_irr_count += 1;
                }
            }
            (Some(Classification::Irr { motif, region }), Classification::Anchor)
            | (Some(Classification::Anchor), Classification::Irr { motif, region }) => {
                // Anchored IRR pair: only the IRR side has a locus to
                // report -- if it's unmapped, there's nothing to push, but
                // the pairing is still resolved (both cache entries
                // consumed).
                if let Some(region) = region {
                    regions_by_motif.entry(motif).or_default().push(region);
                    anchored_irr_count += 1;
                }
            }
            (Some(Classification::Anchor), Classification::Anchor) => {
                // Neither read is an IRR: nothing to count.
            }
            (None, classification) => {
                pending.insert(record.qname().to_vec(), classification);
            }
        }
    }
    // Any reads still pending never had their mate resolved in this scan
    // (e.g. the mate was itself dropped -- neither an IRR candidate nor
    // anchor-tier), so they were never anchored or paired; purely
    // informational for the log message below.
    let unresolved_count = pending.len() as u64;
    drop(pending);

    let irr_count = anchored_irr_count + paired_irr_count;

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
        "sinks: {scanned_count} records scanned, {irr_count} IRR reads counted \
         ({anchored_irr_count} anchored, {paired_irr_count} paired, {unresolved_count} \
         unresolved mates dropped), {} merged regions written to {:?}",
        output_rows.len(),
        args.output,
    );

    Ok(())
}
