//! Small helpers shared by commands that open BAM/CRAM/SAM inputs.

/// Detects a CRAM path/URL by its `.cram` extension (case-insensitive),
/// ignoring any trailing query string or fragment (e.g. S3 presigned URLs).
pub fn is_cram_path(path: &str) -> bool {
    let path_no_query = path.split(['?', '#']).next().unwrap_or(path);
    path_no_query.to_ascii_lowercase().ends_with(".cram")
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
}
