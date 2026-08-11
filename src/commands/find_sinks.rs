//! `find-sinks`: a full linear scan (via a plain `Reader`, no index) for
//! in-repeat reads (IRRs) and the "sink" regions they cluster into.
//!
//! Unlike `profile` (which uses an `IndexedReader` to fetch candidate BED
//! regions and then seek to each candidate's mate), this command never
//! seeks: it reads the whole file once, front to back. That means it can
//! never randomly jump to a read's mate to check the mate's mapping
//! quality. Instead each primary mapped read is classified as it streams by
//! and cached by QNAME; when its mate streams by later (in either order --
//! pairs aren't necessarily read1-then-read2), the two are resolved against
//! each other from the cache and evicted. Because coordinate-sorted BAMs
//! place most mates within one insert-size window of each other, the
//! cache's live size tracks local read density rather than genome size,
//! even though the scan itself covers the whole file.
//!
//! Two further optimizations:
//!   - Unmapped reads are skipped entirely, before any decoding. A read's
//!     own flags already say whether *its* mate is unmapped
//!     (`is_mate_unmapped`), so we never need to observe the unmapped mate
//!     record itself to resolve pairing -- it would only ever be classified
//!     as "other" and cached for no benefit.
//!   - Base/quality decoding (and the purity test) only runs for reads that
//!     already pass the cheap mapped+mapq check, rather than for every
//!     primary read.

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use anyhow::{Context, Result, bail};
use clap::Args;
use rust_htslib::bam::{self, HeaderView, Read as BamRead, Record};
use url::Url;

use crate::bed::Region;
use crate::hts_io::is_cram_path;
use crate::irr;

#[derive(Args, Debug)]
pub struct FindSinksArgs {
    /// Input CRAM/BAM/SAM: local path or s3:// / gs:// / https:// URL.
    /// Scanned front-to-back with a plain (non-indexed) reader.
    #[arg(short = 'i', long)]
    pub input: String,

    /// Output BED path: chrom, start, end, motif, irr_count, anchored_count.
    #[arg(short = 'o', long)]
    pub output: PathBuf,

    /// Reads with MAPQ at or below this are in-repeat-read (IRR) candidates.
    #[arg(long, default_value_t = 40)]
    pub max_irr_mapq: u8,

    /// Minimum MAPQ required of a mapped mate for an IRR to count as
    /// anchored (rather than unanchored).
    #[arg(long, default_value_t = 50)]
    pub min_anchor_mapq: u8,

    /// Shortest repeat-unit (motif) length considered when checking whether
    /// a candidate read is an in-repeat read (IRR).
    #[arg(long, default_value_t = irr::DEFAULT_MOTIF_MIN_LEN)]
    pub motif_min_len: u32,

    /// Longest repeat-unit (motif) length considered.
    #[arg(long, default_value_t = irr::DEFAULT_MOTIF_MAX_LEN)]
    pub motif_max_len: u32,

    /// Reference FASTA. Required when the input is CRAM.
    #[arg(short = 'r', long)]
    pub reference: Option<PathBuf>,

    /// Number of htslib I/O threads to use for decompression.
    #[arg(long, default_value_t = 1)]
    pub threads: usize,
}

/// One candidate IRR read pulled out of the scan: its own reference span,
/// canonical motif, and whether it resolved as anchored.
type Candidate = (Region, Vec<u8>, bool);

/// How a single mapped primary read is classified, independent of pairing.
enum ReadKind {
    /// A read that passed the mapq + purity test, carrying its own
    /// reference span and canonical motif so it can be turned into a
    /// `Candidate` once its anchored/unanchored status is known.
    Irr(Region, Vec<u8>),
    /// A read whose own MAPQ is high enough to anchor a mate IRR.
    Anchor,
    /// Neither of the above: mapq strictly between the two thresholds, or a
    /// low-mapq read that failed the IRR purity test.
    Other,
}

/// The pairing-cache counterpart of `ReadKind` (same shape, but owned and
/// detached from any `Record`, since it may outlive many subsequent reads).
enum CachedKind {
    Irr(Region, Vec<u8>),
    Anchor,
    Other,
}

impl From<ReadKind> for CachedKind {
    fn from(kind: ReadKind) -> Self {
        match kind {
            ReadKind::Irr(region, motif) => CachedKind::Irr(region, motif),
            ReadKind::Anchor => CachedKind::Anchor,
            ReadKind::Other => CachedKind::Other,
        }
    }
}

/// A tiny FNV-1a hasher for the QNAME pairing cache. QNAMEs are short ASCII
/// byte strings hashed on nearly every primary mapped read in the file, so
/// a cheap non-cryptographic hash meaningfully beats the default SipHash on
/// this hot path.
struct FnvHasher(u64);

impl Default for FnvHasher {
    fn default() -> Self {
        FnvHasher(0xcbf2_9ce4_8422_2325)
    }
}

impl Hasher for FnvHasher {
    fn write(&mut self, bytes: &[u8]) {
        const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;
        for &byte in bytes {
            self.0 ^= byte as u64;
            self.0 = self.0.wrapping_mul(FNV_PRIME);
        }
    }

    fn finish(&self) -> u64 {
        self.0
    }
}

type QnameCache = HashMap<Vec<u8>, CachedKind, BuildHasherDefault<FnvHasher>>;

/// A merged sink region: a run of same-motif IRR reads whose spans overlap
/// or touch, with support counts.
struct SinkRegion {
    tid: i32,
    start: i64,
    end: i64,
    motif: Vec<u8>,
    irr_count: usize,
    anchored_count: usize,
}

pub fn run(args: FindSinksArgs) -> Result<()> {
    if is_cram_path(&args.input) && args.reference.is_none() {
        bail!("--reference is required for CRAM input (input={})", args.input);
    }

    let mut reader = match Url::parse(&args.input) {
        Ok(url) => bam::Reader::from_url(&url).with_context(|| format!("failed to open input {}", args.input))?,
        Err(_) => bam::Reader::from_path(&args.input).with_context(|| format!("failed to open input {}", args.input))?,
    };
    if let Some(reference) = &args.reference {
        reader
            .set_reference(reference)
            .with_context(|| format!("failed to set reference {reference:?} on reader"))?;
    }
    if args.threads > 1 {
        reader.set_threads(args.threads).context("failed to set reader thread count")?;
    }

    let mut cache: QnameCache = QnameCache::default();
    let mut candidates: Vec<Candidate> = Vec::new();

    let mut scanned_count: u64 = 0;
    let mut irr_count: u64 = 0;
    let mut anchored_count: u64 = 0;
    let mut unresolved_at_eof: u64 = 0;

    let mut record = Record::new();
    while let Some(result) = reader.read(&mut record) {
        result.context("failed to read record")?;

        // Unmapped reads are never IRR candidates (must be mapped) and
        // never anchors (must be mapped), and any pairing decision that
        // depends on "is my mate unmapped" is already answered by the
        // *mapped* mate's own `is_mate_unmapped` flag -- so there's nothing
        // useful to learn by decoding an unmapped record at all.
        if record.is_secondary() || record.is_supplementary() || record.is_unmapped() {
            continue;
        }
        scanned_count += 1;

        let kind = classify(&record, &args);
        let is_pairable = record.is_paired() && !record.is_mate_unmapped();

        if !is_pairable {
            // No mate to wait for (unpaired read, or a mate we already know
            // is unmapped): resolve immediately without touching the cache.
            if let ReadKind::Irr(region, motif) = kind {
                irr_count += 1;
                candidates.push((region, motif, false));
            }
            continue;
        }

        match cache.remove(record.qname()) {
            Some(cached) => {
                // Both `kind` and `cached` are consumed by their own `if
                // let` below, so snapshot the "is this an anchor" bit each
                // needs from the *other* side up front.
                let current_is_anchor = matches!(&kind, ReadKind::Anchor);
                let cached_is_anchor = matches!(&cached, CachedKind::Anchor);

                if let CachedKind::Irr(region, motif) = cached {
                    irr_count += 1;
                    if current_is_anchor {
                        anchored_count += 1;
                    }
                    candidates.push((region, motif, current_is_anchor));
                }
                if let ReadKind::Irr(region, motif) = kind {
                    irr_count += 1;
                    if cached_is_anchor {
                        anchored_count += 1;
                    }
                    candidates.push((region, motif, cached_is_anchor));
                }
            }
            None => {
                cache.insert(record.qname().to_vec(), kind.into());
            }
        }
    }

    // Anything still cached at EOF is a read whose mate's flags claimed a
    // mapped mate that never actually turned up as a primary alignment
    // (truncated input, unusual filtering upstream, ...). We have no
    // evidence of a qualifying anchor, so treat any leftover IRR entries as
    // unanchored rather than silently dropping them.
    for cached in cache.into_values() {
        if let CachedKind::Irr(region, motif) = cached {
            irr_count += 1;
            unresolved_at_eof += 1;
            candidates.push((region, motif, false));
        }
    }

    let sinks = merge_into_sinks(candidates);
    write_bed(&args.output, &sinks, reader.header())?;

    log::info!(
        "find-sinks: {scanned_count} primary mapped reads scanned, {irr_count} IRR reads \
         ({anchored_count} anchored, {} unanchored, {unresolved_at_eof} unresolved at EOF), \
         {} sink regions written to {:?}",
        irr_count - anchored_count,
        sinks.len(),
        args.output,
    );

    Ok(())
}

/// Classifies a single mapped primary read. Base/quality decoding and the
/// purity test only happen for reads that already pass the mapq gate, since
/// that's the expensive path and most reads in a typical BAM are high-mapq.
fn classify(record: &Record, args: &FindSinksArgs) -> ReadKind {
    let mapq = record.mapq();

    if mapq <= args.max_irr_mapq
        && let Some(motif) =
            irr::classify_in_repeat_read(&record.seq().as_bytes(), record.qual(), args.motif_min_len, args.motif_max_len)
    {
        // `.max(pos + 1)` guards against a degenerate empty CIGAR (no
        // ref-consuming ops), which would otherwise produce an inverted
        // empty region.
        let end = record.cigar().end_pos().max(record.pos() + 1);
        let region = Region { tid: record.tid(), start: record.pos(), end };
        return ReadKind::Irr(region, motif);
    }

    if mapq >= args.min_anchor_mapq {
        return ReadKind::Anchor;
    }

    ReadKind::Other
}

/// Groups candidates by motif, then within each motif group merges spans
/// that overlap or touch into a single `SinkRegion`, accumulating support
/// counts along the way.
fn merge_into_sinks(candidates: Vec<Candidate>) -> Vec<SinkRegion> {
    let mut by_motif: HashMap<Vec<u8>, Vec<(Region, bool)>> = HashMap::new();
    for (region, motif, anchored) in candidates {
        by_motif.entry(motif).or_default().push((region, anchored));
    }

    let mut sinks: Vec<SinkRegion> = Vec::new();
    for (motif, mut items) in by_motif {
        items.sort_by_key(|(region, _)| (region.tid, region.start, region.end));

        for (region, anchored) in items {
            // `sinks` accumulates across motif groups, so the motif check
            // here matters: without it, the first region of a new group
            // could wrongly merge into the previous group's last entry.
            let merges_into_last =
                sinks.last().is_some_and(|last| last.motif == motif && last.tid == region.tid && region.start <= last.end);

            if merges_into_last {
                let last = sinks.last_mut().expect("checked by merges_into_last");
                last.end = last.end.max(region.end);
                last.irr_count += 1;
                if anchored {
                    last.anchored_count += 1;
                }
            } else {
                sinks.push(SinkRegion {
                    tid: region.tid,
                    start: region.start,
                    end: region.end,
                    motif: motif.clone(),
                    irr_count: 1,
                    anchored_count: anchored as usize,
                });
            }
        }
    }

    sinks.sort_by(|a, b| (a.tid, a.start, a.end, &a.motif).cmp(&(b.tid, b.start, b.end, &b.motif)));
    sinks
}

fn write_bed(path: &PathBuf, sinks: &[SinkRegion], header: &HeaderView) -> Result<()> {
    let file = std::fs::File::create(path).with_context(|| format!("failed to create output BED {path:?}"))?;
    let mut writer = BufWriter::new(file);

    for sink in sinks {
        let chrom = String::from_utf8_lossy(header.tid2name(sink.tid as u32));
        let motif = String::from_utf8_lossy(&sink.motif);
        writeln!(writer, "{chrom}\t{}\t{}\t{motif}\t{}\t{}", sink.start, sink.end, sink.irr_count, sink.anchored_count)
            .with_context(|| format!("failed to write output BED {path:?}"))?;
    }

    writer.flush().with_context(|| format!("failed to flush output BED {path:?}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn r(tid: i32, start: i64, end: i64) -> Region {
        Region { tid, start, end }
    }

    fn sink_tuples(sinks: &[SinkRegion]) -> Vec<(i32, i64, i64, String, usize, usize)> {
        sinks
            .iter()
            .map(|s| (s.tid, s.start, s.end, String::from_utf8_lossy(&s.motif).into_owned(), s.irr_count, s.anchored_count))
            .collect()
    }

    #[test]
    fn merge_into_sinks_merges_overlapping_same_motif_candidates() {
        let candidates = vec![
            (r(0, 100, 150), b"CAG".to_vec(), true),
            (r(0, 140, 190), b"CAG".to_vec(), false),
            (r(0, 500, 550), b"CAG".to_vec(), false),
        ];

        let sinks = merge_into_sinks(candidates);
        assert_eq!(sink_tuples(&sinks), vec![(0, 100, 190, "CAG".to_string(), 2, 1), (0, 500, 550, "CAG".to_string(), 1, 0),]);
    }

    #[test]
    fn merge_into_sinks_touching_spans_merge() {
        let candidates = vec![(r(0, 100, 150), b"CAG".to_vec(), false), (r(0, 150, 200), b"CAG".to_vec(), true)];

        let sinks = merge_into_sinks(candidates);
        assert_eq!(sink_tuples(&sinks), vec![(0, 100, 200, "CAG".to_string(), 2, 1)]);
    }

    #[test]
    fn merge_into_sinks_keeps_different_motifs_at_same_coordinates_separate() {
        let candidates = vec![(r(0, 100, 150), b"CAG".to_vec(), false), (r(0, 100, 150), b"GATA".to_vec(), false)];

        let sinks = merge_into_sinks(candidates);
        let mut tuples = sink_tuples(&sinks);
        tuples.sort();
        assert_eq!(
            tuples,
            vec![(0, 100, 150, "CAG".to_string(), 1, 0), (0, 100, 150, "GATA".to_string(), 1, 0)]
        );
    }

    #[test]
    fn merge_into_sinks_keeps_different_contigs_separate() {
        let candidates = vec![(r(0, 100, 150), b"CAG".to_vec(), false), (r(1, 100, 150), b"CAG".to_vec(), false)];

        let sinks = merge_into_sinks(candidates);
        let mut tuples = sink_tuples(&sinks);
        tuples.sort();
        assert_eq!(tuples, vec![(0, 100, 150, "CAG".to_string(), 1, 0), (1, 100, 150, "CAG".to_string(), 1, 0)]);
    }

    #[test]
    fn merge_into_sinks_non_touching_spans_stay_separate() {
        let candidates = vec![(r(0, 100, 150), b"CAG".to_vec(), false), (r(0, 151, 200), b"CAG".to_vec(), false)];

        let sinks = merge_into_sinks(candidates);
        let mut tuples = sink_tuples(&sinks);
        tuples.sort();
        assert_eq!(tuples, vec![(0, 100, 150, "CAG".to_string(), 1, 0), (0, 151, 200, "CAG".to_string(), 1, 0)]);
    }

    #[test]
    fn fnv_hasher_is_deterministic_and_sensitive_to_input() {
        use std::hash::{BuildHasher, BuildHasherDefault};

        let build: BuildHasherDefault<FnvHasher> = BuildHasherDefault::default();
        let h1 = build.hash_one(b"read1");
        let h2 = build.hash_one(b"read1");
        let h3 = build.hash_one(b"read2");
        assert_eq!(h1, h2);
        assert_ne!(h1, h3);
    }
}
