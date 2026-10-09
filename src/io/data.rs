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

pub fn normalize_name(s: &str) -> String {
    s.trim().replace(' ', "_").to_ascii_lowercase()
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
}
