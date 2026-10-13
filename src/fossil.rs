//! Fossil (stratigraphic) ranges: the oldest and youngest ages at which a taxon
//! is known, drawn as wide bars along its branch. Stored on tips as the
//! attributes `fossil_max_age` (first appearance) and `fossil_min_age` (last).

use crate::io::data::{normalize_name, DataTable};
use crate::tree::{Attr, NodeId, Tree};
use anyhow::{bail, Result};
use std::collections::HashMap;

pub const MIN_KEY: &str = "fossil_min_age";
pub const MAX_KEY: &str = "fossil_max_age";

/// A node's fossil range (youngest, oldest) in Ma, if it has one.
pub fn range(tree: &Tree, n: NodeId) -> Option<(f64, f64)> {
    let a = &tree.nodes[n].attrs;
    Some((a.get(MIN_KEY)?.as_f64()?, a.get(MAX_KEY)?.as_f64()?))
}

/// Set a node's fossil range; the two ages may be given in either order.
pub fn set_range(tree: &mut Tree, n: NodeId, a: f64, b: f64) {
    let attrs = &mut tree.nodes[n].attrs;
    attrs.insert(MIN_KEY.into(), Attr::Num(a.min(b)));
    attrs.insert(MAX_KEY.into(), Attr::Num(a.max(b)));
}

pub fn clear_range(tree: &mut Tree, n: NodeId) {
    tree.nodes[n].attrs.remove(MIN_KEY);
    tree.nodes[n].attrs.remove(MAX_KEY);
}

/// Read a table of fossil ranges (comma- or tab-delimited, with a header such as
/// `Species,minimum_age,maximum_age`) onto the matching tips. The age columns
/// are the ones whose names contain "min" and "max", else the two after the
/// taxon column. Returns the number of taxa matched and the names not found.
pub fn import(text: &str, name: &str, tree: &mut Tree) -> Result<(usize, Vec<String>)> {
    // A file without a header starts straight away with a taxon and two ages:
    // give it one, so its first row isn't taken for column names.
    let first = text.lines().find(|l| !l.trim().is_empty()).unwrap_or("");
    let delim = if first.contains('\t') { '\t' } else if first.contains(';') { ';' } else { ',' };
    let fields: Vec<&str> = first.split(delim).map(str::trim).collect();
    let headerless = fields.len() >= 3 && fields[1].parse::<f64>().is_ok() && fields[2].parse::<f64>().is_ok();
    let text = if headerless { format!("taxon{d}age_1{d}age_2\n{}", text, d = delim) } else { text.to_string() };
    let table = DataTable::from_str(&text, name, None)?;
    let col = |word: &str| (0..table.columns.len()).find(|&c| c != table.key && table.columns[c].to_ascii_lowercase().contains(word));
    let others: Vec<usize> = (0..table.columns.len()).filter(|&c| c != table.key).collect();
    let (Some(lo), Some(hi)) = (col("min").or(others.first().copied()), col("max").or(others.get(1).copied())) else {
        bail!("expected a taxon column and two age columns (e.g. Species,minimum_age,maximum_age)");
    };
    if lo == hi {
        bail!("could not tell the minimum and maximum age columns apart");
    }
    let tips: HashMap<String, NodeId> = tree.tips().into_iter().map(|n| (normalize_name(tree.label(n)), n)).collect();
    let (mut matched, mut unmatched) = (0, Vec::new());
    for row in &table.rows {
        let (Ok(a), Ok(b)) = (row[lo].trim().parse::<f64>(), row[hi].trim().parse::<f64>()) else {
            if !row[table.key].trim().is_empty() {
                unmatched.push(format!("{} (ages not numbers)", row[table.key]));
            }
            continue;
        };
        match tips.get(&normalize_name(&row[table.key])) {
            Some(&n) => {
                set_range(tree, n, a, b);
                matched += 1;
            }
            None => unmatched.push(row[table.key].clone()),
        }
    }
    Ok((matched, unmatched))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::newick::parse_newick;

    #[test]
    fn imports_comma_or_tab_tables() {
        let mut t = parse_newick("((Homo_sapiens:2,Pan:2):3,Gorilla:5);").unwrap();
        let (m, un) = import("Species,minimum_age,maximum_age\nHomo sapiens,0,0.3\nGorilla,8,2\nFoo,1,2\n", "f", &mut t).unwrap();
        assert_eq!((m, un), (2, vec!["Foo".to_string()]));
        assert_eq!(range(&t, t.find_label("Homo_sapiens").unwrap()), Some((0.0, 0.3)));
        // Ages given in the wrong order are swapped.
        assert_eq!(range(&t, t.find_label("Gorilla").unwrap()), Some((2.0, 8.0)));
        let mut t2 = parse_newick("(A:1,B:1);").unwrap();
        let (m, _) = import("taxon\tmax\tmin\nA\t5\t1\n", "f", &mut t2).unwrap();
        assert_eq!(m, 1);
        assert_eq!(range(&t2, t2.find_label("A").unwrap()), Some((1.0, 5.0)));
        let a = t2.find_label("A").unwrap();
        clear_range(&mut t2, a);
        assert_eq!(range(&t2, t2.find_label("A").unwrap()), None);
        // No header: the first line is data, not column names.
        let mut t3 = parse_newick("((Smilodectes:1,Notharctus:1):1,Adapis:2);").unwrap();
        let (m, un) = import("Smilodectes	50.3	48.62
Notharctus	50.3	48.62
Adapis	40	36.6
", "f", &mut t3).unwrap();
        assert_eq!((m, un.len()), (3, 0));
        assert_eq!(range(&t3, t3.find_label("Smilodectes").unwrap()), Some((48.62, 50.3)));
    }
}
