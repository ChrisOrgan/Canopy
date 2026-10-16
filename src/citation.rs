//! How to cite Canopy: one source for the README text, the About box and
//! `canopy --cite`. The DOI is Zenodo's concept DOI, which always resolves to the
//! latest version, so citations to every version are counted together.

/// Zenodo concept DOI (all versions).
pub const DOI: &str = "10.5281/zenodo.23425860";

/// The year of the cited release.
pub const YEAR: &str = "2026";

/// The citation in text form, with this build's version.
pub fn text() -> String {
    format!("Organ, C. ({}) Canopy (version {}) [Computer software]. Zenodo. https://doi.org/{}", YEAR, env!("CARGO_PKG_VERSION"), DOI)
}

/// The citation as a BibTeX entry, with this build's version.
pub fn bibtex() -> String {
    format!(
        "@software{{organ_canopy_{year},\n  author    = {{Organ, Chris}},\n  title     = {{Canopy}},\n  year      = {{{year}}},\n  version   = {{{version}}},\n  publisher = {{Zenodo}},\n  doi       = {{{doi}}},\n  url       = {{https://doi.org/{doi}}}\n}}",
        year = YEAR,
        version = env!("CARGO_PKG_VERSION"),
        doi = DOI
    )
}

#[cfg(test)]
mod tests {
    #[test]
    fn citation_names_this_version() {
        let v = env!("CARGO_PKG_VERSION");
        assert!(super::text().contains(v) && super::text().contains(super::DOI));
        assert!(super::bibtex().contains(&format!("version   = {{{}}}", v)));
    }
}
