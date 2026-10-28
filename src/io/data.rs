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
    /// Column types chosen by the user (by column name), overriding detection.
    #[serde(default)]
    pub kinds: std::collections::BTreeMap<String, ColumnKind>,
}

/// New columns made from a numeric column (`DataTable::add_transform`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Transform {
    /// x − mean.
    Center,
    /// (x − mean) / sd: z-scores.
    Standardize,
    /// Natural logarithm.
    Ln,
    /// Base-10 logarithm.
    Log10,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum ColumnKind {
    Numeric,
    Categorical,
}

/// Quantile `p` of sorted values, interpolating between order statistics (R type 7).
pub fn quantile(sorted: &[f64], p: f64) -> f64 {
    let h = (sorted.len() as f64 - 1.0) * p;
    let (lo, hi) = (h.floor() as usize, h.ceil() as usize);
    sorted[lo] + (h - lo as f64) * (sorted[hi] - sorted[lo])
}

/// Histogram of `xs` scaled to a density (bar areas sum to 1): bin edges and
/// heights. Bin width by the Freedman–Diaconis rule, 5 to 50 bins.
pub fn histogram(xs: &[f64]) -> (Vec<f64>, Vec<f64>) {
    histogram_bins(xs, 0)
}

/// As `histogram`, with `bins` equal-width bins (0 = Freedman–Diaconis).
pub fn histogram_bins(xs: &[f64], bins: usize) -> (Vec<f64>, Vec<f64>) {
    let mut s = xs.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let (lo, hi) = (s[0], s[s.len() - 1]);
    if hi <= lo {
        return (vec![lo - 0.5, lo + 0.5], vec![1.0]);
    }
    let iqr = quantile(&s, 0.75) - quantile(&s, 0.25);
    let width = 2.0 * iqr / (s.len() as f64).cbrt();
    let auto = if width > 0.0 { ((hi - lo) / width).ceil() as usize } else { 10 }.clamp(5, 50);
    let bins = if bins > 0 { bins } else { auto };
    let w = (hi - lo) / bins as f64;
    let mut counts = vec![0usize; bins];
    for &x in &s {
        counts[(((x - lo) / w) as usize).min(bins - 1)] += 1;
    }
    let edges = (0..=bins).map(|i| lo + i as f64 * w).collect();
    (edges, counts.iter().map(|&c| c as f64 / (s.len() as f64 * w)).collect())
}

/// Gaussian kernel density of `xs` at `points` evenly spaced values over the
/// data range padded by three bandwidths; bandwidth by Silverman's rule of
/// thumb (R's bw.nrd0).
pub fn kernel_density(xs: &[f64], points: usize) -> Vec<(f64, f64)> {
    let mut s = xs.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = s.len() as f64;
    let mean = s.iter().sum::<f64>() / n;
    let sd = if s.len() > 1 { (s.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0)).sqrt() } else { 0.0 };
    let iqr = quantile(&s, 0.75) - quantile(&s, 0.25);
    let spread = if iqr > 0.0 { sd.min(iqr / 1.34) } else { sd };
    // As R's bw.nrd0 when the spread is 0: |x₁|, else 1.
    let fallback = if s[0] != 0.0 { s[0].abs() } else { 1.0 };
    let bw = (0.9 * if spread > 0.0 { spread } else { fallback } * n.powf(-0.2)).max(1e-12);
    let (lo, hi) = (s[0] - 3.0 * bw, s[s.len() - 1] + 3.0 * bw);
    let norm = 1.0 / (n * bw * (2.0 * std::f64::consts::PI).sqrt());
    (0..points)
        .map(|i| {
            let x = lo + (hi - lo) * i as f64 / (points - 1).max(1) as f64;
            (x, norm * s.iter().map(|v| (-0.5 * ((x - v) / bw).powi(2)).exp()).sum::<f64>())
        })
        .collect()
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
        // No header row (e.g. BayesTraits data files): the first line is a taxon
        // and its values. Columns are then named taxon, trait_1, trait_2, ...
        let first: Vec<&str> = text.lines().find(|l| !l.trim().is_empty()).unwrap_or("").split(delim as char).map(str::trim).collect();
        let missing = |s: &str| matches!(s, "" | "-" | "?" | "NA" | "na" | "N/A");
        let headerless = first.len() >= 2 && first[1..].iter().all(|s| s.parse::<f64>().is_ok() || missing(s)) && first[1..].iter().any(|s| s.parse::<f64>().is_ok());
        let named;
        let text = if headerless {
            let d = delim as char;
            let names: Vec<String> = (1..first.len()).map(|i| if first.len() == 2 { "trait".to_string() } else { format!("trait_{}", i) }).collect();
            named = format!("taxon{}{}
{}", d, names.join(&d.to_string()), text);
            named.as_str()
        } else {
            text
        };
        // Row numbers as in a spreadsheet (row 1 is the header, if there is one).
        let row0 = if headerless { 1 } else { 2 };
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
            match first_row.get(&k) {
                Some(&j) => dups.push(format!("{} (rows {} and {})", r[key], j + row0, i + row0)),
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
        Ok(DataTable { name: name.to_string(), columns, rows, key, kinds: Default::default() })
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

    /// Column type: the user's choice if any, else numeric when every value parses.
    pub fn kind(&self, col: usize) -> ColumnKind {
        if let Some(&k) = self.kinds.get(&self.columns[col]) {
            return k;
        }
        self.detected_kind(col)
    }

    /// Numeric if every non-missing value parses as a number.
    pub fn detected_kind(&self, col: usize) -> ColumnKind {
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
                    // A column set to numeric treats values that aren't numbers as missing.
                    ColumnKind::Numeric => match v.parse() {
                        Ok(x) => Attr::Num(x),
                        Err(_) => continue,
                    },
                    ColumnKind::Categorical => Attr::Text(v.to_string()),
                };
                tree.nodes[n].attrs.insert(col.clone(), attr);
            }
        }
        (matched, unmatched)
    }

    /// Remove a trait column (not the taxon column).
    pub fn remove_column(&mut self, col: usize) {
        if col == self.key || col >= self.columns.len() {
            return;
        }
        let name = self.columns.remove(col);
        for r in self.rows.iter_mut() {
            if col < r.len() {
                r.remove(col);
            }
        }
        if col < self.key {
            self.key -= 1;
        }
        self.kinds.remove(&name);
    }

    /// Add a column of `col` transformed (centered, standardized or logged),
    /// named after it (e.g. "mass_z"). Missing values, values that aren't
    /// numbers and (for logs) values ≤ 0 are left empty. Returns the new
    /// column's name and how many values were left empty that weren't before.
    pub fn add_transform(&mut self, col: usize, kind: Transform) -> Result<(String, usize)> {
        let parse = |s: &str| s.trim().parse::<f64>().ok().filter(|x| x.is_finite());
        let xs: Vec<f64> = self.rows.iter().filter_map(|r| parse(&r[col])).collect();
        if xs.is_empty() {
            bail!("\"{}\" has no numbers to transform", self.columns[col]);
        }
        let n = xs.len() as f64;
        let mean = xs.iter().sum::<f64>() / n;
        // Sample standard deviation (n − 1), as in R's scale().
        let sd = if xs.len() > 1 { (xs.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0)).sqrt() } else { 0.0 };
        if kind == Transform::Standardize && sd <= 0.0 {
            bail!("\"{}\" does not vary, so it cannot be standardized", self.columns[col]);
        }
        let suffix = match kind {
            Transform::Center => "centered",
            Transform::Standardize => "z",
            Transform::Ln => "ln",
            Transform::Log10 => "log10",
        };
        let base = format!("{}_{}", self.columns[col], suffix);
        let mut name = base.clone();
        let mut i = 2;
        while self.columns.contains(&name) {
            name = format!("{}_{}", base, i);
            i += 1;
        }
        let mut lost = 0;
        for r in self.rows.iter_mut() {
            let v = match parse(&r[col]) {
                Some(x) => match kind {
                    Transform::Center => Some(x - mean),
                    Transform::Standardize => Some((x - mean) / sd),
                    Transform::Ln if x > 0.0 => Some(x.ln()),
                    Transform::Log10 if x > 0.0 => Some(x.log10()),
                    _ => {
                        lost += 1;
                        None
                    }
                },
                None => None,
            };
            r.push(v.map(|x| format!("{}", (x * 1e10).round() / 1e10)).unwrap_or_default());
        }
        self.columns.push(name.clone());
        Ok((name, lost))
    }

    /// Summary statistics of each trait column, as a table: N, missing, mean,
    /// SD (sample), min, 2.5%, 5%, Q1, median, Q3, 95%, 97.5% and max for numeric columns (quartiles
    /// interpolated as in R's default); the number of states and the most common
    /// one for categorical columns. With `only`, just the rows of those taxa
    /// (normalized names, e.g. the tips of the tree).
    pub fn summary(&self, only: Option<&std::collections::HashSet<String>>) -> (Vec<String>, Vec<Vec<String>>) {
        let header: Vec<String> = ["Column", "Type", "N", "Missing", "Mean", "SD", "Min", "2.5%", "5%", "Q1", "Median", "Q3", "95%", "97.5%", "Max", "States", "Most common"].iter().map(|s| s.to_string()).collect();
        let rows: Vec<&Vec<String>> = self.rows.iter().filter(|r| only.map(|s| s.contains(&normalize_name(&r[self.key]))).unwrap_or(true)).collect();
        let missing = |v: &str| v.is_empty() || v.eq_ignore_ascii_case("na") || v == "?" || v == "-";
        let f = |x: f64| crate::tree::format_num(x, 3);
        let mut out = Vec::new();
        for (c, (name, kind)) in (0..self.columns.len()).filter(|&c| c != self.key).map(|c| (c, (self.columns[c].clone(), self.kind(c)))) {
            let vals: Vec<&str> = rows.iter().map(|r| r[c].trim()).collect();
            let present: Vec<&str> = vals.iter().copied().filter(|v| !missing(v)).collect();
            let mut row = vec![name, String::new(), String::new(), (vals.len() - present.len()).to_string()];
            match kind {
                ColumnKind::Numeric => {
                    let mut xs: Vec<f64> = present.iter().filter_map(|v| v.parse::<f64>().ok()).filter(|x| x.is_finite()).collect();
                    xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
                    row[1] = "numeric".into();
                    row[2] = xs.len().to_string();
                    row[3] = (vals.len() - xs.len()).to_string();
                    if xs.is_empty() {
                        row.extend(std::iter::repeat(String::new()).take(13));
                    } else {
                        let n = xs.len() as f64;
                        let mean = xs.iter().sum::<f64>() / n;
                        let sd = if xs.len() > 1 { (xs.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0)).sqrt() } else { f64::NAN };
                        // Quantile with linear interpolation between order statistics (R type 7).
                        let q = |p: f64| quantile(&xs, p);
                        row.extend([f(mean), if sd.is_finite() { f(sd) } else { String::new() }, f(xs[0]), f(q(0.025)), f(q(0.05)), f(q(0.25)), f(q(0.5)), f(q(0.75)), f(q(0.95)), f(q(0.975)), f(xs[xs.len() - 1]), String::new(), String::new()]);
                    }
                }
                ColumnKind::Categorical => {
                    let mut counts: std::collections::BTreeMap<&str, usize> = std::collections::BTreeMap::new();
                    for v in &present {
                        *counts.entry(v).or_default() += 1;
                    }
                    row[1] = "categorical".into();
                    row[2] = present.len().to_string();
                    row.extend(std::iter::repeat(String::new()).take(11));
                    row.push(counts.len().to_string());
                    // Most common state (ties: the first alphabetically).
                    let top = counts.iter().fold(None, |best: Option<(&str, usize)>, (k, &n)| match best {
                        Some((_, b)) if b >= n => best,
                        _ => Some((k, n)),
                    });
                    row.push(top.map(|(k, n)| format!("{} ({})", k, n)).unwrap_or_default());
                }
            }
            out.push(row);
        }
        (header, out)
    }

    /// The numbers in column `col` (optionally only the rows of the taxa in `only`).
    pub fn numbers(&self, col: usize, only: Option<&std::collections::HashSet<String>>) -> Vec<f64> {
        self.rows
            .iter()
            .filter(|r| only.map(|s| s.contains(&normalize_name(&r[self.key]))).unwrap_or(true))
            .filter_map(|r| r[col].trim().parse::<f64>().ok())
            .filter(|x| x.is_finite())
            .collect()
    }

    /// The states of column `col` with their counts, most common first.
    pub fn state_counts(&self, col: usize, only: Option<&std::collections::HashSet<String>>) -> Vec<(String, usize)> {
        let mut counts: std::collections::BTreeMap<String, usize> = std::collections::BTreeMap::new();
        for r in self.rows.iter().filter(|r| only.map(|s| s.contains(&normalize_name(&r[self.key]))).unwrap_or(true)) {
            let v = r[col].trim();
            if !(v.is_empty() || v.eq_ignore_ascii_case("na") || v == "?" || v == "-") {
                *counts.entry(v.to_string()).or_default() += 1;
            }
        }
        let mut v: Vec<(String, usize)> = counts.into_iter().collect();
        v.sort_by(|a, b| b.1.cmp(&a.1));
        v
    }

    /// The table without the rows for the named taxa (matched as tip labels are).
    pub fn without(&self, names: &[String]) -> DataTable {
        let drop: std::collections::HashSet<String> = names.iter().map(|n| normalize_name(n)).collect();
        let mut d = self.clone();
        d.rows.retain(|r| !drop.contains(&normalize_name(&r[d.key])));
        d
    }

    /// Taxa in the tree with no row in the table, and rows naming no taxon in the tree.
    pub fn mismatches(&self, tree: &Tree) -> (Vec<String>, Vec<String>) {
        let rows: std::collections::HashSet<String> = self.rows.iter().map(|r| normalize_name(&r[self.key])).filter(|k| !k.is_empty()).collect();
        let tips: Vec<String> = tree.tips().into_iter().map(|n| tree.label(n).to_string()).collect();
        let tip_keys: std::collections::HashSet<String> = tips.iter().map(|t| normalize_name(t)).collect();
        let no_row = tips.into_iter().filter(|t| !rows.contains(&normalize_name(t))).collect();
        let no_tip = self.rows.iter().map(|r| r[self.key].clone()).filter(|k| !k.trim().is_empty() && !tip_keys.contains(&normalize_name(k))).collect();
        (no_row, no_tip)
    }

    pub fn trait_columns(&self) -> Vec<(String, ColumnKind)> {
        (0..self.columns.len()).filter(|&c| c != self.key).map(|c| (self.columns[c].clone(), self.kind(c))).collect()
    }
}

/// True if tip `n` has no usable value for column `col`: no row in the table,
/// an empty/NA cell, or a non-numeric entry in a numeric column.
pub fn is_missing(tree: &Tree, n: usize, col: &str) -> bool {
    match tree.nodes[n].attrs.get(col) {
        None => true,
        Some(Attr::Num(x)) => !x.is_finite(),
        Some(Attr::Text(s)) => s.trim().is_empty(),
        Some(_) => false,
    }
}

/// Tips missing a value for any of `cols`.
pub fn tips_missing(tree: &Tree, cols: &[String]) -> Vec<usize> {
    tree.tips().into_iter().filter(|&n| cols.iter().any(|c| is_missing(tree, n, c))).collect()
}

/// Table of the given tips: label, distance to the root, number of nodes on
/// the path to the root (internal nodes including the root, i.e. the number
/// of splits the lineage passed through), the tip's branch (`tip_length`, then
/// its mean, median and 95% HPD across a tree set), its height (then mean,
/// median and HPD), then every other attribute. Numbers are rounded to 3 decimals.
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
    // Summaries of a quantity in the order value, mean, median, 95% HPD.
    let group = |keys: &mut Vec<String>, stem: &str| -> Vec<String> {
        let names = ["_mean", "_median", "_95%_HPD"].map(|s| format!("{}{}", stem, s));
        let found: Vec<String> = names.into_iter().filter(|k| keys.contains(k)).collect();
        keys.retain(|k| !found.contains(k));
        found
    };
    let tip_spread = group(&mut keys, "length");
    let height_spread = group(&mut keys, "height");
    let height_attr = keys.iter().position(|k| k == "height").map(|i| keys.remove(i));
    // Any other quantity with summaries: its mean, median and HPD together, in that order.
    let stems: Vec<String> = keys.iter().filter_map(|k| k.strip_suffix("_95%_HPD").or(k.strip_suffix("_median")).or(k.strip_suffix("_mean")).map(str::to_string)).collect();
    for stem in stems {
        let found = group(&mut keys, &stem);
        if found.is_empty() {
            continue;
        }
        let at = keys.iter().position(|k| *k == stem).map(|i| i + 1).unwrap_or(keys.len());
        for (j, k) in found.into_iter().enumerate() {
            keys.insert(at + j, k);
        }
    }
    let lengths = tree.has_lengths();
    // Without a `height` annotation, the tree's own: time before the youngest tip.
    let computed_height = height_attr.is_none() && lengths;
    let heights = if computed_height { tree.heights() } else { Vec::new() };
    let mut header = vec!["label".to_string(), "root_distance".to_string(), "nodes_to_root".to_string()];
    if lengths {
        header.push("tip_length".to_string());
    }
    header.extend(tip_spread.iter().map(|k| format!("tip_{}", k)));
    if computed_height || height_attr.is_some() {
        header.push("height".to_string());
    }
    header.extend(height_spread.iter().cloned());
    header.extend(keys.iter().cloned());
    let value = |n: usize, k: &str| tree.nodes[n].attrs.get(k).map(|a| table_value(k, a)).unwrap_or_default();
    let rows = tips
        .iter()
        .map(|&n| {
            let mut r = vec![tree.label(n).to_string(), crate::tree::format_num(depths[n], 3), (tree.path_to_root(n).len() - 1).to_string()];
            if lengths {
                r.push(tree.nodes[n].length.map(|l| crate::tree::format_num(l, 3)).unwrap_or_default());
            }
            r.extend(tip_spread.iter().map(|k| value(n, k)));
            if computed_height {
                r.push(table_value("height", &Attr::Num(heights[n])));
            } else if height_attr.is_some() {
                r.push(value(n, "height"));
            }
            r.extend(height_spread.iter().map(|k| value(n, k)));
            r.extend(keys.iter().map(|k| value(n, k)));
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
    fn reads_tables_without_a_header() {
        let t = DataTable::from_str("Homo_sapiens	62
Pan_troglodytes	45
", "x", None).unwrap();
        assert_eq!(t.columns, vec!["taxon", "trait"]);
        assert_eq!(t.rows.len(), 2);
        assert_eq!(t.rows[0], vec!["Homo_sapiens", "62"]);
        let t = DataTable::from_str("A 1 -
B 2 3
", "x", Some(b' ')).unwrap();
        assert_eq!(t.columns, vec!["taxon", "trait_1", "trait_2"]);
        // A header row stays a header.
        let t = DataTable::from_str("species,mass
A,1
", "x", None).unwrap();
        assert_eq!(t.columns, vec!["species", "mass"]);
        assert_eq!(t.rows.len(), 1);
    }

    #[test]
    fn table_columns() {
        let mut t = parse_newick("((A:1.23456,B:1):2,C:3);").unwrap();
        let a = t.find_label("A").unwrap();
        t.nodes[a].attrs.insert("x".into(), Attr::Num(0.123456));
        let (h, rows) = tips_table(&t, &[a, t.find_label("C").unwrap()]);
        assert_eq!(h[..6], ["label", "root_distance", "nodes_to_root", "tip_length", "height", "x"]);
        // A is the deepest tip (3.235), so its height is 0; C is 0.235 above it.
        assert_eq!(rows[0], vec!["A", "3.235", "2", "1.235", "0", "0.123"]);
        assert_eq!(rows[1][4], "0.235");
        assert_eq!(rows[1][2], "1");
        t.nodes[a].attrs.insert("height".into(), Attr::Num(3.2e-7));
        t.nodes[a].attrs.insert("height_95%_HPD".into(), Attr::Range(-1e-6, 0.5));
        let (h, rows) = tips_table(&t, &[a]);
        let col = |name: &str| h.iter().position(|x| x == name).unwrap();
        assert_eq!(rows[0][col("height")], "0");
        assert_eq!(rows[0][col("height_95%_HPD")], "[0, 0.5]");
    }

    #[test]
    fn missing_data() {
        let mut t = parse_newick("((A,B),(C,D));").unwrap();
        let d = DataTable::from_str("species,mass,diet\nA,70,omni\nB,NA,frugi\nC,12,\n", "d", None).unwrap();
        d.join_to_tree(&mut t);
        let names = |v: Vec<usize>| v.into_iter().map(|n| t.label(n).to_string()).collect::<Vec<_>>();
        // D has no row; B lacks mass; C lacks diet.
        assert_eq!(names(tips_missing(&t, &["mass".into()])), ["B", "D"]);
        assert_eq!(names(tips_missing(&t, &["diet".into()])), ["C", "D"]);
        assert_eq!(names(tips_missing(&t, &["mass".into(), "diet".into()])), ["B", "C", "D"]);
        assert!(tips_missing(&t, &[]).is_empty());
    }

    #[test]
    fn summary_columns_in_order() {
        let t = crate::io::newick::parse_newick("(A[&height=1,height_95%_HPD={0.5,2},height_median=0.9,height_mean=1.1,length_95%_HPD={1,3},length_mean=1.5,length_median=1.4,mass=3]:1,B[&height=1]:1);").unwrap();
        let (header, rows) = tips_table(&t, &t.tips());
        // Each quantity: value, mean, median, 95% HPD; the tip's branch first, renamed.
        let want = ["tip_length", "tip_length_mean", "tip_length_median", "tip_length_95%_HPD", "height", "height_mean", "height_median", "height_95%_HPD", "mass"];
        assert_eq!(header[3..], want, "{:?}", header);
        assert_eq!(rows[0][3..8], ["1", "1.5", "1.4", "[1, 3]", "1"]);
    }

    #[test]
    fn transforms_and_deleting_columns() {
        let mut d = DataTable::from_str("species,x,y
A,1,a
B,2,b
C,3,c
D,NA,d
", "d", None).unwrap();
        // Mean 2, sample sd 1: centered and z-scores are both -1, 0, 1.
        assert_eq!(d.add_transform(1, Transform::Center).unwrap(), ("x_centered".to_string(), 0));
        assert_eq!(d.add_transform(1, Transform::Standardize).unwrap().0, "x_z");
        let col = |d: &DataTable, name: &str| { let c = d.columns.iter().position(|x| x == name).unwrap(); d.rows.iter().map(|r| r[c].clone()).collect::<Vec<_>>() };
        assert_eq!(col(&d, "x_centered"), ["-1", "0", "1", ""]);
        assert_eq!(col(&d, "x_z"), ["-1", "0", "1", ""]);
        let (name, _) = d.add_transform(1, Transform::Log10).unwrap();
        assert_eq!(col(&d, &name)[0], "0");
        // A second centered column gets a new name; log of 0 is left empty.
        assert_eq!(d.add_transform(1, Transform::Center).unwrap().0, "x_centered_2");
        let mut z = DataTable::from_str("species,v
A,0
B,1
", "z", None).unwrap();
        assert_eq!(z.add_transform(1, Transform::Ln).unwrap().1, 1);
        assert!(d.add_transform(2, Transform::Center).is_err());
        // Deleting a column.
        let n = d.columns.len();
        d.remove_column(2);
        assert_eq!(d.columns.len(), n - 1);
        assert!(!d.columns.contains(&"y".to_string()) && d.rows.iter().all(|r| r.len() == n - 1));
        d.remove_column(d.key);
        assert_eq!(d.columns.len(), n - 1);
    }

    #[test]
    fn summary_statistics() {
        let d = DataTable::from_str("species,x,diet
A,1,omni
B,2,frugi
C,3,omni
D,4,
E,NA,omni
", "d", None).unwrap();
        let (h, rows) = d.summary(None);
        let get = |r: usize, k: &str| rows[r][h.iter().position(|x| x == k).unwrap()].clone();
        // x: 1-4, mean 2.5, sample SD 1.291, quartiles 1.75 / 2.5 / 3.25 (R type 7).
        assert_eq!((get(0, "N"), get(0, "Missing"), get(0, "Mean"), get(0, "SD")), ("4".into(), "1".into(), "2.5".into(), "1.291".into()));
        assert_eq!((get(0, "Q1"), get(0, "Median"), get(0, "Q3"), get(0, "Max")), ("1.75".into(), "2.5".into(), "3.25".into(), "4".into()));
        assert_eq!((get(1, "Type"), get(1, "N"), get(1, "Missing"), get(1, "States"), get(1, "Most common")), ("categorical".into(), "4".into(), "1".into(), "2".into(), "omni (3)".into()));
        // Only some taxa.
        let only: std::collections::HashSet<String> = ["a", "b"].iter().map(|s| s.to_string()).collect();
        let (_, r2) = d.summary(Some(&only));
        assert_eq!(r2[0][4], "1.5");
    }

    #[test]
    fn quantiles_histogram_and_density() {
        let xs: Vec<f64> = (1..=100).map(|i| i as f64).collect();
        // R: quantile(1:100, c(.025, .975)) = 3.475, 97.525.
        assert!((quantile(&xs, 0.025) - 3.475).abs() < 1e-9 && (quantile(&xs, 0.975) - 97.525).abs() < 1e-9);
        let (edges, h) = histogram(&xs);
        let area: f64 = h.iter().zip(edges.windows(2)).map(|(d, e)| d * (e[1] - e[0])).sum();
        assert!((area - 1.0).abs() < 1e-9 && h.len() >= 5);
        let k = kernel_density(&xs, 400);
        let step = k[1].0 - k[0].0;
        let integral: f64 = k.iter().map(|p| p.1 * step).sum();
        assert!((integral - 1.0).abs() < 0.02, "{}", integral);
        let peak = k.iter().cloned().fold((0.0, 0.0), |a, b| if b.1 > a.1 { b } else { a });
        assert!((peak.0 - 50.5).abs() < 10.0);
    }

    #[test]
    fn mismatches_and_without() {
        let t = crate::io::newick::parse_newick("((A,B),C);").unwrap();
        let d = DataTable::from_str("species,x
A,1
C,2
Z,3
", "d", None).unwrap();
        assert_eq!(d.mismatches(&t), (vec!["B".to_string()], vec!["Z".to_string()]));
        let d2 = d.without(&["a".to_string()]);
        assert_eq!(d2.rows.len(), 2);
        assert_eq!(d2.mismatches(&t).0, ["A", "B"]);
    }

    #[test]
    fn chosen_column_kinds() {
        let mut t = parse_newick("((A,B),C);").unwrap();
        let mut d = DataTable::from_str("species,state,size\nA,0,big\nB,1,3\nC,1,4\n", "d", None).unwrap();
        assert_eq!(d.kind(1), ColumnKind::Numeric);
        assert_eq!(d.kind(2), ColumnKind::Categorical);
        // 0/1 states read as categories; a mostly numeric column read as numbers.
        d.kinds.insert("state".into(), ColumnKind::Categorical);
        d.kinds.insert("size".into(), ColumnKind::Numeric);
        assert_eq!(d.kind(1), ColumnKind::Categorical);
        d.join_to_tree(&mut t);
        let a = t.find_label("A").unwrap();
        let b = t.find_label("B").unwrap();
        assert_eq!(t.nodes[a].attrs["state"], Attr::Text("0".into()));
        assert_eq!(t.nodes[b].attrs["size"], Attr::Num(3.0));
        // "big" is not a number, so A's size counts as missing.
        assert!(!t.nodes[a].attrs.contains_key("size"));
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
