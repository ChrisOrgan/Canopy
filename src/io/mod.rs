//! File input/output.

pub mod data;
pub mod hyphy;
pub mod newick;
pub mod nexus;
pub mod phyloxml;

use crate::tree::Tree;
use anyhow::{Context, Result};
use std::path::Path;

/// Read every tree in a Newick, NEXUS or phyloXML file (format detected from content).
pub fn read_trees_str(text: &str) -> Result<Vec<Tree>> {
    if hyphy::is_hyphy_json(text) {
        return Ok(vec![hyphy::parse(text)?.0]);
    }
    if phyloxml::is_phyloxml(text) {
        phyloxml::parse_phyloxml(text)
    } else if nexus::is_nexus(text) {
        Ok(nexus::parse_nexus(text)?.trees)
    } else {
        newick::parse_newick_multi(text)
    }
}

pub fn read_trees(path: &Path) -> Result<Vec<Tree>> {
    let bytes = std::fs::read(path).with_context(|| format!("reading {}", path.display()))?;
    let text = String::from_utf8_lossy(&bytes);
    read_trees_str(&text).with_context(|| format!("parsing {}", path.display()))
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FileKind {
    Tree,
    Data,
    /// BayesTraits output: a VarRates file (variable rates or Fabric) or an MCMC log.
    VarRates,
    BayesLog,
    Image,
    Unknown,
}

pub fn classify(path: &Path) -> FileKind {
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("").to_ascii_lowercase();
    match ext.as_str() {
        "nwk" | "newick" | "tre" | "tree" | "trees" | "nex" | "nexus" | "nxs" | "t" | "con" | "treefile"
        | "contree" | "phy" | "xml" | "phyloxml" => FileKind::Tree,
        "csv" | "tsv" | "tab" => FileKind::Data,
        // HyPhy results (aBSREL, BUSTED, MEME, FEL, RELAX, …).
        "json" => match std::fs::read_to_string(path) {
            Ok(s) if hyphy::is_hyphy_json(&s) => FileKind::Tree,
            _ => FileKind::Unknown,
        },
        "png" | "jpg" | "jpeg" => FileKind::Image,
        "txt" => {
            // Sniff: a tree file starts with '(' or #NEXUS.
            match std::fs::read_to_string(path) {
                Ok(s) if s.trim_start().starts_with('(') || nexus::is_nexus(&s) => FileKind::Tree,
                Ok(s) if crate::bayestraits::is_varrates(&s) => FileKind::VarRates,
                Ok(s) if crate::bayestraits::is_log(&s) => FileKind::BayesLog,
                Ok(_) => FileKind::Data,
                Err(_) => FileKind::Unknown,
            }
        }
        _ => FileKind::Unknown,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn recognizes_bayestraits_output() {
        let ex = |f: &str| Path::new(env!("CARGO_MANIFEST_DIR")).join("examples").join(f);
        assert_eq!(classify(&ex("primates.MCC.VR.VarRates.txt")), FileKind::VarRates);
        assert_eq!(classify(&ex("primates.MCC.FAB.VarRates.txt")), FileKind::VarRates);
        assert_eq!(classify(&ex("primates.MCC.VR.Log.txt")), FileKind::BayesLog);
        assert_eq!(classify(&ex("primate_traits.txt")), FileKind::Data);
        for f in ["primates_100_NeocortexSA_EndocastVol.txt", "ToyTreeData.tsv", "adapid_fossil_ranges.tsv"] {
            assert_eq!(classify(&ex(f)), FileKind::Data, "{}", f);
        }
        // A trait table with columns named like BayesTraits output is still data.
        let dir = std::env::temp_dir().join("canopy-classify-test");
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("traits.txt");
        std::fs::write(&p, "Iteration	Node ID	mass
A	1	2
").unwrap();
        assert_eq!(classify(&p), FileKind::Data);
    }
}
