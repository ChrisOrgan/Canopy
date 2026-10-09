//! Comparative (trait) data tables joined to tree tips by label.

use crate::tree::{Attr, Tree};
use anyhow::{bail, Result};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DataTable {
    pub name: String,
    pub columns: Vec<String>,
    pub rows: Vec<Vec<String>>,
    /// Column holding taxon names.
    pub key: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ColumnKind {
    Numeric,
    Categorical,
}

/// Taxon name as used for matching: case ignored, any run of whitespace
/// (including tabs and non-breaking spaces) or underscores becomes one underscore.
pub fn normalize_name(s: &str) -> String {
    s.split(|c: char| c.is_whitespace() || c == '_').filter(|w| !w.is_empty()).collect::<Vec<_>>().join("_").to_lowercase()
}

impl DataTable {
    pub fn from_str(text: &str, name: &str, delimiter: Option<u8>) -> Result<Self> {
        let delim = delimiter.unwrap_or_else(|| {
            let first = text.lines().next().unwrap_or("");
            let tabs = first.matches('\t').count();
            let commas = first.matches(',').count();
            let semis = first.matches(';').count();
            if tabs >= commas && tabs >= semis && tabs > 0 {
                b'\t'
            } else if semis > commas {
                b';'
            } else {
                b','
            }
        });
        let mut rdr = csv::ReaderBuilder::new().delimiter(delim).flexible(true).has_headers(true).from_reader(text.as_bytes());
        let columns: Vec<String> = rdr.headers()?.iter().map(|s| s.trim().to_string()).collect();
        if columns.is_empty() {
            bail!("data table has no columns");
        }
        let mut rows = Vec::new();
        for rec in rdr.records() {
            let rec = rec?;
            let mut row: Vec<String> = rec.iter().map(|s| s.trim().to_string()).collect();
            row.resize(columns.len(), String::new());
            rows.push(row);
        }
        let key = columns
            .iter()
            .position(|c| matches!(c.to_ascii_lowercase().as_str(), "label" | "taxon" | "taxa" | "species" | "tip" | "name" | "id"))
            .unwrap_or(0);
        // Each taxon may appear once: a second row would silently overwrite the first.
        let mut first_row: HashMap<String, usize> = HashMap::new();
        let mut dups: Vec<String> = Vec::new();
        for (i, r) in rows.iter().enumerate() {
            let k = normalize_name(&r[key]);
            if k.is_empty() {
                continue;
            }
            // Row numbers as in a spreadsheet: the header is row 1.
            match first_row.get(&k) {
                Some(&j) => dups.push(format!("{} (rows {} and {})", r[key], j + 2, i + 2)),
                None => {
                    first_row.insert(k, i);
                }
            }
        }
        if !dups.is_empty() {
            let more = if dups.len() > 10 { format!(" and {} more", dups.len() - 10) } else { String::new() };
            bail!(
                "{} taxa appear in more than one row of column \"{}\": {}{}. Keep one row per taxon (e.g. average the values) and load it again.",
                dups.len(),
                columns[key],
                dups.iter().take(10).cloned().collect::<Vec<_>>().join(", "),
                more
            );
        }
        Ok(DataTable { name: name.to_string(), columns, rows, key })
    }

    pub fn from_path(path: &Path) -> Result<Self> {
        let text = std::fs::read_to_string(path)?;
        let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("").to_ascii_lowercase();
        let delim = match ext.as_str() {
            "tsv" | "tab" => Some(b'\t'),
            "csv" => None,
            _ => None,
        };
        let name = path.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
        Self::from_str(&text, &name, delim)
    }

    pub fn kind(&self, col: usize) -> ColumnKind {
        let mut any = false;
        for r in &self.rows {
            let v = r[col].trim();
            if v.is_empty() || v.eq_ignore_ascii_case("na") || v == "?" || v == "-" {
                continue;
            }
            any = true;
            if v.parse::<f64>().is_err() {
                return ColumnKind::Categorical;
            }
        }
        if any {
            ColumnKind::Numeric
        } else {
            ColumnKind::Categorical
        }
    }

    /// Copy table columns onto matching tips as attributes. Returns (matched, unmatched rows).
    pub fn join_to_tree(&self, tree: &mut Tree) -> (usize, Vec<String>) {
        let mut tips: HashMap<String, usize> = HashMap::new();
        for n in tree.tips() {
            tips.insert(normalize_name(tree.label(n)), n);
        }
        let kinds: Vec<ColumnKind> = (0..self.columns.len()).map(|c| self.kind(c)).collect();
        let mut matched = 0;
        let mut unmatched = Vec::new();
        for row in &self.rows {
            let key = normalize_name(&row[self.key]);
            let Some(&n) = tips.get(&key) else {
                unmatched.push(row[self.key].clone());
                continue;
            };
            matched += 1;
            for (c, col) in self.columns.iter().enumerate() {
                if c == self.key {
                    continue;
                }
                let v = row[c].trim();
                if v.is_empty() || v.eq_ignore_ascii_case("na") || v == "?" {
                    continue;
                }
                let attr = match kinds[c] {
                    ColumnKind::Numeric => Attr::Num(v.parse().unwrap_or(f64::NAN)),
                    ColumnKind::Categorical => Attr::Text(v.to_string()),
                };
                tree.nodes[n].attrs.insert(col.clone(), attr);
            }
        }
        (matched, unmatched)
    }

    pub fn trait_columns(&self) -> Vec<(String, ColumnKind)> {
        (0..self.columns.len()).filter(|&c| c != self.key).map(|c| (self.columns[c].clone(), self.kind(c))).collect()
    }
}

/// Table of the given tips: label, distance to the root, number of nodes on
/// the path to the root (internal nodes including the root, i.e. the number
/// of splits the lineage passed through), then every attribute any of them
/// carries. Numbers are rounded to 3 decimals.
pub fn tips_table(tree: &Tree, tips: &[usize]) -> (Vec<String>, Vec<Vec<String>>) {
    let depths = tree.depths(true);
    let mut keys: Vec<String> = Vec::new();
    for &n in tips {
        for k in tree.nodes[n].attrs.keys() {
            if !keys.contains(k) {
                keys.push(k.clone());
            }
        }
    }
    // Without a `height` annotation (e.g. a single non-BEAST tree), compute it
    // from branch lengths: time before the youngest tip.
    let computed_height = !keys.iter().any(|k| k == "height") && tree.has_lengths();
    let heights = if computed_height { tree.heights() } else { Vec::new() };
    let mut header = vec!["label".to_string(), "root_distance".to_string(), "nodes_to_root".to_string()];
    if computed_height {
        header.push("height".to_string());
    }
    header.extend(keys.iter().cloned());
    let rows = tips
        .iter()
        .map(|&n| {
            let mut r = vec![
                tree.label(n).to_string(),
                crate::tree::format_num(depths[n], 3),
                (tree.path_to_root(n).len() - 1).to_string(),
            ];
            if computed_height {
                r.push(table_value("height", &Attr::Num(heights[n])));
            }
            r.extend(keys.iter().map(|k| tree.nodes[n].attrs.get(k).map(|a| table_value(k, a)).unwrap_or_default()));
            r
        })
        .collect();
    (header, rows)
}

/// One row per internal node (root excluded) for bipartition / clade support
/// reports, as from RAxML: node id, number of tips, each support value the
/// tree carries (posterior, bootstrap, ...), branch length, height when the
/// tree has lengths, and the taxa in the clade. Returns (header, rows, node ids).
pub fn node_report(tree: &Tree) -> (Vec<String>, Vec<Vec<String>>, Vec<usize>) {
    let keys = crate::ops::support_keys(tree);
    let lengths = tree.has_lengths();
    let heights = if lengths { tree.heights() } else { Vec::new() };
    let mut header = vec!["node".to_string(), "tips".to_string()];
    header.extend(keys.iter().map(|k| k.to_string()));
    header.push("branch_length".into());
    if lengths {
        header.push("height".into());
    }
    header.push("taxa".into());
    let mut rows = Vec::new();
    let mut ids = Vec::new();
    for n in tree.preorder() {
        if tree.is_tip(n) || n == tree.root {
            continue;
        }
        let tips = tree.tips_below(n);
        let mut r = vec![n.to_string(), tips.len().to_string()];
        r.extend(keys.iter().map(|k| tree.nodes[n].attrs.get(*k).map(|a| table_value(k, a)).unwrap_or_default()));
        r.push(tree.nodes[n].length.map(|l| crate::tree::format_num(l, 3)).unwrap_or_default());
        if lengths {
            r.push(table_value("height", &Attr::Num(heights[n])));
        }
        r.push(tips.iter().map(|&t| tree.label(t)).collect::<Vec<_>>().join(" "));
        rows.push(r);
        ids.push(n);
    }
    (header, rows, ids)
}

/// Variance–covariance matrix of the tips below `n` (`ops::vcv`) as a table:
/// a blank corner, then one column and one row per tip. `digits` rounds the
/// values for display (values below 0.0001 show as 0); None keeps full precision for export.
pub fn vcv_table(tree: &Tree, n: usize, digits: Option<usize>) -> (Vec<String>, Vec<Vec<String>>) {
    let (tips, m) = crate::ops::vcv(tree, n);
    let mut header = vec![String::new()];
    header.extend(tips.iter().map(|&t| tree.label(t).to_string()));
    let cell = |x: f64| match digits {
        Some(d) => crate::tree::format_num(if x.abs() < 1e-4 { 0.0 } else { x }, d),
        None => x.to_string(),
    };
    let rows = tips.iter().zip(&m).map(|(&t, row)| std::iter::once(tree.label(t).to_string()).chain(row.iter().map(|&x| cell(x))).collect()).collect();
    (header, rows)
}

/// A cell value rounded to 3 decimals. Heights below 0.0001 (rounding noise
/// for tips at the present) are shown as 0.
pub fn table_value(key: &str, a: &Attr) -> String {
    let clean = |x: f64| if key.starts_with("height") && x.abs() < 1e-4 { 0.0 } else { x };
    match a {
        Attr::Num(x) => crate::tree::format_num(clean(*x), 3),
        Attr::Range(lo, hi) => Attr::Range(clean(*lo), clean(*hi)).display(3),
        other => other.display(3),
    }
}

/// Render a table as delimited text (tab for spreadsheets, comma for CSV).
pub fn table_to_text(header: &[String], rows: &[Vec<String>], delimiter: u8) -> String {
    let mut w = csv::WriterBuilder::new().delimiter(delimiter).from_writer(Vec::new());
    let _ = w.write_record(header);
    for r in rows {
        let _ = w.write_record(r);
    }
    String::from_utf8(w.into_inner().unwrap_or_default()).unwrap_or_default()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::newick::parse_newick;

    #[test]
    fn table_columns() {
        let mut t = parse_newick("((A:1.23456,B:1):2,C:3);").unwrap();
        let a = t.find_label("A").unwrap();
        t.nodes[a].attrs.insert("x".into(), Attr::Num(0.123456));
        let (h, rows) = tips_table(&t, &[a, t.find_label("C").unwrap()]);
        assert_eq!(h[..5], ["label", "root_distance", "nodes_to_root", "height", "x"]);
        // A is the deepest tip (3.235), so its height is 0; C is 0.235 above it.
        assert_eq!(rows[0], vec!["A", "3.235", "2", "0", "0.123"]);
        assert_eq!(rows[1][3], "0.235");
        assert_eq!(rows[1][2], "1");
        t.nodes[a].attrs.insert("height".into(), Attr::Num(3.2e-7));
        t.nodes[a].attrs.insert("height_95%_HPD".into(), Attr::Range(-1e-6, 0.5));
        let (h, rows) = tips_table(&t, &[a]);
        let col = |name: &str| h.iter().position(|x| x == name).unwrap();
        assert_eq!(rows[0][col("height")], "0");
        assert_eq!(rows[0][col("height_95%_HPD")], "[0, 0.5]");
    }

    #[test]
    fn join() {
        let mut t = parse_newick("((Homo_sapiens,Pan),Gorilla);").unwrap();
        let d = DataTable::from_str("species,mass,diet\nHomo sapiens,70,omni\nPan,45,frugi\nFoo,1,x\n", "d", None).unwrap();
        assert_eq!(d.kind(1), ColumnKind::Numeric);
        assert_eq!(d.kind(2), ColumnKind::Categorical);
        let (m, un) = d.join_to_tree(&mut t);
        assert_eq!(m, 2);
        assert_eq!(un, vec!["Foo".to_string()]);
        let h = t.find_label("Homo_sapiens").unwrap();
        assert_eq!(t.nodes[h].attrs["mass"], Attr::Num(70.0));
    }

    #[test]
    fn node_report_lists_support() {
        let t = parse_newick("((A:0.1,B:0.2):0.3[100],(C:0.1,D:0.1):0.2[75],E:0.4);").unwrap();
        let (h, rows, ids) = node_report(&t);
        assert_eq!(h, ["node", "tips", "support", "branch_length", "height", "taxa"]);
        assert_eq!(rows.len(), 2);
        assert_eq!(ids.len(), 2);
        assert_eq!(rows[0][2], "100");
        assert_eq!(rows[1][2], "75");
        assert_eq!(rows[1][5], "C D");
    }

    #[test]
    fn duplicate_taxa_rejected() {
        let e = DataTable::from_str("species,mass\nHomo sapiens,70\nPan,45\nhomo_sapiens,65\n", "d", None).unwrap_err().to_string();
        assert!(e.contains("homo_sapiens (rows 2 and 4)"), "{}", e);
        // Non-breaking and doubled spaces (common in spreadsheet exports) still count.
        assert!(DataTable::from_str("species,mass\nHomo\u{a0}sapiens,70\nHomo  sapiens ,65\n", "d", None).is_err());
        assert!(DataTable::from_str("species,mass\nA,1\n,2\n,3\n", "d", None).is_ok(), "blank names are not duplicates");
    }
}
