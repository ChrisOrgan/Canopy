//! Checks GitHub for a newer Canopy release.

use anyhow::Result;
use std::time::Duration;

/// The latest release on GitHub.
const LATEST: &str = "https://api.github.com/repos/ChrisOrgan/Canopy/releases/latest";

/// A release newer than the running version.
#[derive(Clone, Debug)]
pub struct Release {
    /// Version without the leading "v", e.g. "1.7.0".
    pub version: String,
    /// The release page.
    pub url: String,
}

/// Numeric parts of a version tag ("v1.6.5" → [1, 6, 5]); pre-release and
/// build suffixes are ignored.
fn parts(v: &str) -> Vec<u64> {
    let v = v.trim().trim_start_matches(['v', 'V']);
    let v = v.split(['-', '+']).next().unwrap_or("");
    v.split('.').map(|p| p.parse().unwrap_or(0)).collect()
}

/// True when version `a` is newer than `b` (missing parts count as 0).
pub fn newer(a: &str, b: &str) -> bool {
    let (a, b) = (parts(a), parts(b));
    for i in 0..a.len().max(b.len()) {
        let (x, y) = (a.get(i).copied().unwrap_or(0), b.get(i).copied().unwrap_or(0));
        if x != y {
            return x > y;
        }
    }
    false
}

/// The latest GitHub release if it is newer than this build, None if not.
pub fn check() -> Result<Option<Release>> {
    let text = ureq::get(LATEST)
        .set("User-Agent", concat!("Canopy/", env!("CARGO_PKG_VERSION")))
        .set("Accept", "application/vnd.github+json")
        .timeout(Duration::from_secs(10))
        .call()?
        .into_string()?;
    let body: serde_json::Value = serde_json::from_str(&text)?;
    let tag = body["tag_name"].as_str().unwrap_or("");
    if tag.is_empty() || !newer(tag, env!("CARGO_PKG_VERSION")) {
        return Ok(None);
    }
    let url = body["html_url"].as_str().unwrap_or("https://github.com/ChrisOrgan/Canopy/releases/latest").to_string();
    Ok(Some(Release { version: tag.trim_start_matches(['v', 'V']).to_string(), url }))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compares_versions() {
        assert!(newer("v1.7.0", "1.6.5"));
        assert!(newer("1.6.10", "1.6.9"));
        assert!(newer("v2", "1.9.9"));
        assert!(!newer("v1.6.5", "1.6.5"));
        assert!(!newer("1.6", "1.6.0"));
        assert!(!newer("v1.6.4", "1.6.5"));
        assert!(!newer("v1.7.0-beta", "1.7.0"));
    }
}
