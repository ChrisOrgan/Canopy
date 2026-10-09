//! File input/output.

pub mod data;
pub mod newick;
pub mod nexus;

use crate::tree::Tree;
use anyhow::{Context, Result};
use std::path::Path;

/// Read every tree in a Newick or NEXUS file (format detected from content).
pub fn read_trees_str(text: &str) -> Result<Vec<Tree>> {
    if nexus::is_nexus(text) {
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
    Image,
    Project,
    Unknown,
}

pub fn classify(path: &Path) -> FileKind {
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("").to_ascii_lowercase();
    match ext.as_str() {
        "nwk" | "newick" | "tre" | "tree" | "trees" | "nex" | "nexus" | "nxs" | "t" | "con" | "treefile"
        | "contree" | "phy" => FileKind::Tree,
        "csv" | "tsv" | "tab" => FileKind::Data,
        "png" | "jpg" | "jpeg" => FileKind::Image,
        "canopy" => FileKind::Project,
        "txt" => {
            // Sniff: a tree file starts with '(' or #NEXUS.
            match std::fs::read_to_string(path) {
                Ok(s) if s.trim_start().starts_with('(') || nexus::is_nexus(&s) => FileKind::Tree,
                Ok(_) => FileKind::Data,
                Err(_) => FileKind::Unknown,
            }
        }
        _ => FileKind::Unknown,
    }
}
