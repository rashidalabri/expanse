//! Minimal reader for CRAM `.crai` index files.
//!
//! htslib already consults the CRAI when it opens a CRAM (to answer
//! `fetch()`), but it doesn't expose the parsed slice boundaries through
//! rust-htslib's safe API. We read and parse the CRAI ourselves so callers
//! can group BED regions by the underlying CRAM slice they land in, and
//! issue one `fetch()` per slice instead of one per region.

use std::ffi::CString;
use std::io::Read;
use std::os::raw::c_void;

use anyhow::{Context, Result, bail};
use rust_htslib::htslib;

use crate::bed::Region;

/// A safety ceiling on CRAI size: real indexes are a tiny fraction of the
/// CRAM they describe, so hitting this means something is wrong rather than
/// that we should silently truncate the index.
const MAX_CRAI_BYTES: usize = 512 * 1024 * 1024;

/// One slice entry from a CRAI: the genomic span it covers (0-based,
/// half-open) and its byte location, which we don't need beyond using it to
/// identify distinct slices while parsing.
#[derive(Debug, Clone, Copy)]
pub struct Slice {
    pub tid: i32,
    pub start: i64,
    pub end: i64,
}

/// Load and parse the CRAI index for `input` (a local path or s3:// / gs:// /
/// https:// URL). Slices covering only unmapped reads (refID < 0) are
/// dropped since they can never overlap a BED region.
pub fn load(input: &str) -> Result<Vec<Slice>> {
    let idx_path = crai_path(input);
    let raw = unsafe { read_all(&idx_path) }
        .with_context(|| format!("failed to read CRAI index {idx_path}"))?;

    let decompressed = if raw.len() >= 2 && raw[0] == 0x1f && raw[1] == 0x8b {
        let mut out = Vec::new();
        flate2::read::MultiGzDecoder::new(&raw[..])
            .read_to_end(&mut out)
            .with_context(|| format!("failed to gunzip CRAI index {idx_path}"))?;
        out
    } else {
        raw
    };

    parse(&decompressed).with_context(|| format!("failed to parse CRAI index {idx_path}"))
}

/// Expand each region to the full span of every CRAM slice it overlaps.
/// Regions that don't land on any known slice (e.g. an empty/stale index)
/// are passed through unchanged, so callers always cover at least as much as
/// the original region.
pub fn expand_to_slices(regions: &[Region], slices: &[Slice]) -> Vec<Region> {
    let mut sorted: Vec<&Slice> = slices.iter().collect();
    sorted.sort_by_key(|s| (s.tid, s.start, s.end));

    regions
        .iter()
        .flat_map(|region| {
            let lo = sorted.partition_point(|s| s.tid < region.tid);
            let hi = sorted.partition_point(|s| s.tid <= region.tid);

            let mut matches = Vec::new();
            for slice in &sorted[lo..hi] {
                if slice.start >= region.end {
                    break;
                }
                if slice.end > region.start {
                    matches.push(Region {
                        tid: slice.tid,
                        start: slice.start,
                        end: slice.end,
                    });
                }
            }

            if matches.is_empty() {
                vec![*region]
            } else {
                matches
            }
        })
        .collect()
}

fn crai_path(input: &str) -> String {
    match input.find(['?', '#']) {
        Some(idx) => format!("{}.crai{}", &input[..idx], &input[idx..]),
        None => format!("{input}.crai"),
    }
}

fn parse(data: &[u8]) -> Result<Vec<Slice>> {
    let text = std::str::from_utf8(data).context("CRAI index is not valid UTF-8 text")?;

    let mut slices = Vec::new();
    for (line_no, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        let mut fields = line.split_whitespace();
        let mut next_i64 = |name: &str| -> Result<i64> {
            fields
                .next()
                .with_context(|| format!("CRAI line {}: missing {name}", line_no + 1))?
                .parse::<i64>()
                .with_context(|| format!("CRAI line {}: invalid {name}", line_no + 1))
        };

        let ref_id = next_i64("refID")?;
        let start = next_i64("start")?;
        let span = next_i64("span")?;
        let _container_offset = next_i64("containerOffset")?;
        let _slice_offset = next_i64("sliceOffset")?;
        let _slice_size = next_i64("sliceSize")?;

        // Unmapped-reads-only slices (refID == -1) never overlap a BED
        // region, which are always on a real contig.
        if ref_id < 0 || span <= 0 {
            continue;
        }

        slices.push(Slice {
            tid: ref_id as i32,
            // CRAI start is 1-based; convert to the 0-based half-open
            // convention used by `Region` elsewhere in this crate.
            start: start - 1,
            end: start - 1 + span,
        });
    }

    Ok(slices)
}

unsafe fn read_all(path: &str) -> Result<Vec<u8>> {
    let c_path = CString::new(path).with_context(|| format!("invalid CRAI path {path:?}"))?;
    let c_mode = c"r";

    let fp = unsafe { htslib::hopen(c_path.as_ptr(), c_mode.as_ptr()) };
    if fp.is_null() {
        bail!("htslib could not open {path}");
    }

    let mut buf = Vec::new();
    let mut chunk = [0u8; 65536];
    let read_result = (|| -> Result<()> {
        loop {
            let n = unsafe { hread(fp, &mut chunk)? };
            if n == 0 {
                return Ok(());
            }
            buf.extend_from_slice(&chunk[..n]);
            if buf.len() > MAX_CRAI_BYTES {
                bail!("CRAI index exceeds {MAX_CRAI_BYTES}-byte safety limit");
            }
        }
    })();

    unsafe {
        htslib::hclose(fp);
    }
    read_result?;
    Ok(buf)
}

unsafe extern "C" {
    // htslib's `hread()` is `static inline` in hfile.h, so bindgen doesn't
    // generate a binding for it, but the `hread2` helper it delegates to
    // (declared as a local `extern` inside that inline function) is a
    // normal exported symbol in libhts that we can link against directly.
    fn hread2(fp: *mut htslib::hFILE, buffer: *mut c_void, nbytes: usize, nread: usize) -> isize;
}

/// Reimplementation of htslib's `hread()`: drain whatever is already
/// buffered in `fp`, then pull the remainder from the backend via
/// `hread2`. Safe to call repeatedly on the same handle, exactly like the
/// real `hread()` (this is the same pattern htslib's own CRAI loader uses).
unsafe fn hread(fp: *mut htslib::hFILE, buf: &mut [u8]) -> Result<usize> {
    unsafe {
        let begin = (*fp).begin;
        let end = (*fp).end;
        let avail = (end as usize).saturating_sub(begin as usize);
        let n = avail.min(buf.len());
        if n > 0 {
            std::ptr::copy_nonoverlapping(begin as *const u8, buf.as_mut_ptr(), n);
            (*fp).begin = begin.add(n);
        }

        if n == buf.len() || htslib::hFILE::mobile_raw(fp) == 0 {
            Ok(n)
        } else {
            let r = hread2(fp, buf.as_mut_ptr() as *mut c_void, buf.len(), n);
            if r < 0 {
                bail!("htslib read error");
            }
            Ok(r as usize)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn s(tid: i32, start: i64, end: i64) -> Slice {
        Slice { tid, start, end }
    }

    fn r(tid: i32, start: i64, end: i64) -> Region {
        Region { tid, start, end }
    }

    #[test]
    fn crai_path_appends_suffix() {
        assert_eq!(crai_path("foo.cram"), "foo.cram.crai");
        assert_eq!(
            crai_path("s3://bucket/foo.cram?X-Amz-Signature=abc"),
            "s3://bucket/foo.cram.crai?X-Amz-Signature=abc"
        );
    }

    #[test]
    fn parse_reads_slice_rows_and_converts_to_half_open() {
        // refid start span containerOffset sliceOffset sliceSize
        let data = b"0\t1\t100\t0\t0\t50\n0\t101\t50\t1000\t0\t20\n1\t1\t10\t2000\t0\t5\n-1\t0\t0\t3000\t0\t5\n";
        let slices = parse(data).unwrap();
        assert_eq!(
            slices
                .iter()
                .map(|s| (s.tid, s.start, s.end))
                .collect::<Vec<_>>(),
            vec![(0, 0, 100), (0, 100, 150), (1, 0, 10)]
        );
    }

    #[test]
    fn parse_skips_blank_lines() {
        let data = b"0\t1\t100\t0\t0\t50\n\n  \n0\t101\t50\t1000\t0\t20\n";
        let slices = parse(data).unwrap();
        assert_eq!(slices.len(), 2);
    }

    #[test]
    fn expand_to_slices_replaces_region_with_overlapping_slice_spans() {
        let slices = vec![s(0, 0, 100), s(0, 100, 200), s(1, 0, 50)];
        let regions = vec![r(0, 10, 20), r(0, 90, 110)];
        let expanded = expand_to_slices(&regions, &slices);
        assert_eq!(expanded, vec![r(0, 0, 100), r(0, 0, 100), r(0, 100, 200)]);
    }

    #[test]
    fn expand_to_slices_falls_back_to_original_region_when_unmatched() {
        let slices = vec![s(0, 0, 100)];
        let regions = vec![r(1, 10, 20)];
        let expanded = expand_to_slices(&regions, &slices);
        assert_eq!(expanded, vec![r(1, 10, 20)]);
    }

    #[test]
    fn expand_to_slices_empty_index_is_noop() {
        let regions = vec![r(0, 10, 20), r(1, 5, 15)];
        let expanded = expand_to_slices(&regions, &[]);
        assert_eq!(expanded, regions);
    }
}
