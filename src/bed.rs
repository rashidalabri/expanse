use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{bail, Context, Result};
use rust_htslib::bam::HeaderView;

/// A half-open genomic interval `[start, end)` on a contig identified by its
/// numeric target id (`tid`) in a BAM/CRAM header.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Region {
    pub tid: i32,
    pub start: i64,
    pub end: i64,
}

/// Parse a BED file (BED3+, extra columns ignored) into `Region`s, resolving
/// each contig name to its `tid` via the given BAM/CRAM header.
pub fn parse_bed(path: &Path, header: &HeaderView) -> Result<Vec<Region>> {
    let file = File::open(path).with_context(|| format!("failed to open BED file {path:?}"))?;
    let reader = BufReader::new(file);

    let mut regions = Vec::new();
    for (line_no, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("failed to read {path:?} at line {}", line_no + 1))?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
            continue;
        }

        let mut fields = line.split('\t');
        let chrom = fields
            .next()
            .with_context(|| format!("{path:?}:{}: missing chrom column", line_no + 1))?;
        let start: i64 = fields
            .next()
            .with_context(|| format!("{path:?}:{}: missing start column", line_no + 1))?
            .parse()
            .with_context(|| format!("{path:?}:{}: invalid start coordinate", line_no + 1))?;
        let end: i64 = fields
            .next()
            .with_context(|| format!("{path:?}:{}: missing end column", line_no + 1))?
            .parse()
            .with_context(|| format!("{path:?}:{}: invalid end coordinate", line_no + 1))?;

        let tid = header
            .tid(chrom.as_bytes())
            .with_context(|| format!("{path:?}:{}: contig {chrom:?} not found in BAM/CRAM header", line_no + 1))?;

        if end <= start {
            bail!("{path:?}:{}: empty or invalid interval {chrom}:{start}-{end}", line_no + 1);
        }

        regions.push(Region {
            tid: tid as i32,
            start,
            end,
        });
    }

    Ok(regions)
}

/// Find the index of the region in `regions` that contains `(tid, pos)`,
/// or `None` if there isn't one. `regions` must be sorted and
/// non-overlapping per contig (i.e. the output of [`merge_regions`] or
/// [`merge_within`]).
pub fn locate(regions: &[Region], tid: i32, pos: i64) -> Option<usize> {
    let idx = regions.partition_point(|r| (r.tid, r.start) <= (tid, pos));
    if idx == 0 {
        return None;
    }
    let region = regions[idx - 1];
    (region.tid == tid && region.start <= pos && pos < region.end).then_some(idx - 1)
}

/// Total bp of overlap between `[start, end)` on `tid` and `regions`.
/// `regions` must be sorted and non-overlapping per contig (i.e. the output
/// of [`merge_regions`] or [`merge_within`]).
pub fn overlap_length(regions: &[Region], tid: i32, start: i64, end: i64) -> i64 {
    let idx = regions.partition_point(|r| (r.tid, r.end) <= (tid, start));

    let mut total = 0i64;
    for region in &regions[idx..] {
        if region.tid != tid || region.start >= end {
            break;
        }
        total += region.end.min(end) - region.start.max(start);
    }
    total
}

/// Sort and merge overlapping/touching regions, per contig.
pub fn merge_regions(regions: &[Region]) -> Vec<Region> {
    merge_within(regions, 0)
}

/// Sort and merge regions that overlap or are within `distance` bp of each
/// other, per contig. `distance = 0` is equivalent to [`merge_regions`]
/// (overlapping/touching only).
pub fn merge_within(regions: &[Region], distance: i64) -> Vec<Region> {
    let mut sorted: Vec<Region> = regions.to_vec();
    sorted.sort_by_key(|r| (r.tid, r.start, r.end));

    let mut merged: Vec<Region> = Vec::with_capacity(sorted.len());
    for region in sorted {
        match merged.last_mut() {
            Some(last) if last.tid == region.tid && region.start <= last.end + distance => {
                last.end = last.end.max(region.end);
            }
            _ => merged.push(region),
        }
    }
    merged
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn r(tid: i32, start: i64, end: i64) -> Region {
        Region { tid, start, end }
    }

    fn test_header() -> HeaderView {
        let sam_header = b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000000\n@SQ\tSN:chr2\tLN:2000000\n";
        HeaderView::from_bytes(sam_header)
    }

    fn write_bed(contents: &str) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(format!("expanse-test-{}.bed", uuid_like()));
        let mut f = File::create(&path).unwrap();
        f.write_all(contents.as_bytes()).unwrap();
        path
    }

    fn uuid_like() -> u64 {
        use std::time::{SystemTime, UNIX_EPOCH};
        SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos() as u64
    }

    #[test]
    fn parse_bed_basic() {
        let header = test_header();
        let path = write_bed("chr1\t10\t20\nchr2\t100\t200\textra\tcolumns\n# a comment\n\nchr1\t30\t40\n");
        let regions = parse_bed(&path, &header).unwrap();
        std::fs::remove_file(&path).ok();
        assert_eq!(
            regions,
            vec![r(0, 10, 20), r(1, 100, 200), r(0, 30, 40)]
        );
    }

    #[test]
    fn parse_bed_unknown_contig_errors() {
        let header = test_header();
        let path = write_bed("chrUnknown\t10\t20\n");
        let result = parse_bed(&path, &header);
        std::fs::remove_file(&path).ok();
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("chrUnknown"));
    }

    #[test]
    fn parse_bed_invalid_interval_errors() {
        let header = test_header();
        let path = write_bed("chr1\t20\t10\n");
        let result = parse_bed(&path, &header);
        std::fs::remove_file(&path).ok();
        assert!(result.is_err());
    }

    #[test]
    fn parse_bed_malformed_line_errors() {
        let header = test_header();
        let path = write_bed("chr1\t10\n");
        let result = parse_bed(&path, &header);
        std::fs::remove_file(&path).ok();
        assert!(result.is_err());
    }

    #[test]
    fn merge_overlapping() {
        let regions = vec![r(0, 10, 20), r(0, 15, 25), r(0, 30, 40)];
        let merged = merge_regions(&regions);
        assert_eq!(merged, vec![r(0, 10, 25), r(0, 30, 40)]);
    }

    #[test]
    fn merge_touching() {
        // half-open intervals [10,20) and [20,30) touch at 20 and should merge
        let regions = vec![r(0, 10, 20), r(0, 20, 30)];
        let merged = merge_regions(&regions);
        assert_eq!(merged, vec![r(0, 10, 30)]);
    }

    #[test]
    fn merge_disjoint() {
        let regions = vec![r(0, 10, 20), r(0, 21, 30)];
        let merged = merge_regions(&regions);
        assert_eq!(merged, vec![r(0, 10, 20), r(0, 21, 30)]);
    }

    #[test]
    fn merge_multi_contig() {
        let regions = vec![r(1, 10, 20), r(0, 10, 20), r(1, 15, 25), r(0, 25, 35)];
        let merged = merge_regions(&regions);
        assert_eq!(merged, vec![r(0, 10, 20), r(0, 25, 35), r(1, 10, 25)]);
    }

    #[test]
    fn merge_unsorted_input() {
        let regions = vec![r(0, 30, 40), r(0, 10, 20)];
        let merged = merge_regions(&regions);
        assert_eq!(merged, vec![r(0, 10, 20), r(0, 30, 40)]);
    }

    #[test]
    fn merge_empty() {
        let merged = merge_regions(&[]);
        assert!(merged.is_empty());
    }

    #[test]
    fn locate_returns_index_of_containing_region() {
        let regions = merge_regions(&[r(0, 10, 20), r(0, 30, 40), r(1, 5, 15)]);
        assert_eq!(locate(&regions, 0, 10), Some(0));
        assert_eq!(locate(&regions, 0, 35), Some(1));
        assert_eq!(locate(&regions, 1, 10), Some(2));
        assert_eq!(locate(&regions, 0, 20), None);
        assert_eq!(locate(&regions, 2, 10), None);
    }

    #[test]
    fn merge_within_merges_regions_closer_than_distance() {
        // 20..30 and gap of 5 (30..35) before 35..40: within distance 5.
        let regions = vec![r(0, 20, 30), r(0, 35, 40)];
        assert_eq!(merge_within(&regions, 5), vec![r(0, 20, 40)]);
        assert_eq!(merge_within(&regions, 4), vec![r(0, 20, 30), r(0, 35, 40)]);
    }

    #[test]
    fn merge_within_zero_distance_matches_merge_regions() {
        let regions = vec![r(0, 10, 20), r(0, 20, 30), r(0, 32, 40)];
        assert_eq!(merge_within(&regions, 0), merge_regions(&regions));
    }

    #[test]
    fn merge_within_does_not_bridge_across_contigs() {
        let regions = vec![r(0, 10, 20), r(1, 25, 35)];
        assert_eq!(
            merge_within(&regions, 1_000_000),
            vec![r(0, 10, 20), r(1, 25, 35)]
        );
    }

    #[test]
    fn overlap_length_no_regions() {
        assert_eq!(overlap_length(&[], 0, 10, 20), 0);
    }

    #[test]
    fn overlap_length_partial_and_full_overlap() {
        let regions = merge_regions(&[r(0, 10, 20), r(0, 30, 50)]);
        // Partial overlap with the first region, none with the second.
        assert_eq!(overlap_length(&regions, 0, 15, 25), 5);
        // Fully contains the second region.
        assert_eq!(overlap_length(&regions, 0, 25, 60), 20);
        // Spans both regions plus the gap between them.
        assert_eq!(overlap_length(&regions, 0, 0, 60), 30);
    }

    #[test]
    fn overlap_length_no_overlap() {
        let regions = merge_regions(&[r(0, 10, 20)]);
        assert_eq!(overlap_length(&regions, 0, 20, 30), 0);
        assert_eq!(overlap_length(&regions, 0, 0, 10), 0);
    }

    #[test]
    fn overlap_length_ignores_other_contigs() {
        let regions = merge_regions(&[r(0, 10, 20), r(1, 10, 20)]);
        assert_eq!(overlap_length(&regions, 2, 10, 20), 0);
        assert_eq!(overlap_length(&regions, 1, 10, 20), 10);
    }

    #[test]
    fn merge_within_chains_transitively() {
        // Each gap is within distance 10 of its neighbor, chaining all three
        // into one region even though the first and last are 40 bp apart.
        let regions = vec![r(0, 0, 10), r(0, 20, 30), r(0, 40, 50)];
        assert_eq!(merge_within(&regions, 10), vec![r(0, 0, 50)]);
    }
}
