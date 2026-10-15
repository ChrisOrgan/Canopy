//! HyPhy results files (JSON from aBSREL, BUSTED, MEME, FEL, RELAX, …): the
//! tree they analysed with each branch's results attached as attributes.
//!
//! Layout (HyPhy 2.5): `input.trees` holds one Newick tree per partition, with
//! internal nodes named (`Node1`, …); `branch attributes` has, per partition,
//! an object keyed by node name with that branch's values (e.g. "Corrected
//! P-value", "Baseline MG94xREV omega ratio"); `tested` marks the branches
//! tested ("test") and the rest ("background"). Site-level results (MEME, FEL)
//! are not branch values and are left out.

use crate::io::newick::parse_newick;
use crate::tree::{Attr, Tree};
use anyhow::{bail, Context, Result};
use serde_json::Value;

/// Whether `text` looks like a HyPhy results file.
pub fn is_hyphy_json(text: &str) -> bool {
    let t = text.trim_start();
    t.starts_with('{') && t.contains("\"branch attributes\"") && t.contains("\"input\"")
}

/// Name of the method (e.g. "aBSREL") from the `analysis` block, if given.
fn method(v: &Value) -> Option<String> {
    let info = v["analysis"]["info"].as_str()?;
    // The info text usually starts with the method's name.
    let first = info.split(|c: char| c == '(' || c == '\n' || c == ':' || c == '.').next()?.trim();
    (!first.is_empty() && first.len() < 60).then(|| first.to_string())
}

/// A HyPhy results file's tree (the first partition) with branch results as
/// attributes; the method's name if the file gives it.
pub fn parse(text: &str) -> Result<(Tree, Option<String>)> {
    let v: Value = serde_json::from_str(text).context("not valid JSON")?;
    let trees = &v["input"]["trees"];
    let newick = match trees {
        Value::Object(m) => m.get("0").or_else(|| m.values().next()).and_then(Value::as_str),
        Value::Array(a) => a.first().and_then(Value::as_str),
        Value::String(s) => Some(s.as_str()),
        _ => None,
    };
    let Some(newick) = newick else { bail!("no tree under \"input\" › \"trees\"") };
    let newick = if newick.trim_end().ends_with(';') { newick.to_string() } else { format!("{};", newick.trim_end()) };
    let mut tree = parse_newick(&newick).context("reading the tree in the HyPhy file")?;
    let branches = match &v["branch attributes"] {
        Value::Object(m) => m.get("0").or_else(|| m.iter().find(|(k, _)| *k != "attributes").map(|(_, x)| x)),
        _ => None,
    };
    let tested = match &v["tested"] {
        Value::Object(m) => m.get("0").or_else(|| m.values().next()),
        _ => None,
    };
    for n in tree.preorder() {
        let Some(name) = tree.nodes[n].label.clone() else { continue };
        if let Some(Value::Object(vals)) = branches.and_then(|b| b.get(&name)) {
            for (k, x) in vals {
                let attr = match x {
                    Value::Number(num) => num.as_f64().map(Attr::Num),
                    Value::String(s) if k != "original name" => Some(Attr::Text(s.clone())),
                    Value::Bool(b) => Some(Attr::Text(b.to_string())),
                    _ => None,
                };
                if let Some(a) = attr {
                    tree.nodes[n].attrs.insert(k.clone(), a);
                }
            }
            // HyPhy renames taxa it cannot use; show the original name.
            if let Some(Value::String(orig)) = vals.get("original name") {
                if tree.is_tip(n) && !orig.is_empty() {
                    tree.nodes[n].label = Some(orig.clone());
                }
            }
        }
        if let Some(Value::String(t)) = tested.and_then(|x| x.get(&name)) {
            tree.nodes[n].attrs.insert("tested".into(), Attr::Text(t.clone()));
        }
        // HyPhy's internal names (Node1, …) only key the results: keep them as
        // an attribute rather than a label drawn on the tree.
        if !tree.is_tip(n) {
            tree.nodes[n].attrs.insert("hyphy_node".into(), Attr::Text(name));
            tree.nodes[n].label = None;
        }
    }
    tree.name = method(&v);
    Ok((tree, method(&v)))
}

#[cfg(test)]
mod tests {
    use super::*;

    const ABSREL: &str = r#"{
      "analysis": {"info": "aBSREL (Adaptive branch-site random effects likelihood) uses an adaptive random effects branch-site model"},
      "input": {"file name": "x.fas", "number of sequences": 3, "trees": {"0": "((Homo_sapiens:0.1,Pan_troglodytes:0.12)Node2:0.05,Gorilla:0.2)"}},
      "branch attributes": {
        "0": {
          "Homo_sapiens": {"Corrected P-value": 0.012, "Baseline MG94xREV omega ratio": 0.8, "Rate classes": 2, "original name": "Homo sapiens", "Rate Distributions": [[0.1, 0.9], [5.0, 0.1]]},
          "Pan_troglodytes": {"Corrected P-value": 1, "Rate classes": 1, "original name": "Pan_troglodytes"},
          "Node2": {"Corrected P-value": 0.4, "Rate classes": 1}
        },
        "attributes": {"Corrected P-value": {"attribute type": "branch label"}}
      },
      "tested": {"0": {"Homo_sapiens": "test", "Pan_troglodytes": "test", "Node2": "background", "Gorilla": "background"}}
    }"#;

    #[test]
    fn reads_absrel_results() {
        assert!(is_hyphy_json(ABSREL));
        let (t, m) = parse(ABSREL).unwrap();
        assert_eq!(m.as_deref(), Some("aBSREL"));
        // Original names restored; values attached; arrays left out.
        let h = t.find_label("Homo sapiens").unwrap();
        assert_eq!(t.nodes[h].attrs["Corrected P-value"], Attr::Num(0.012));
        assert_eq!(t.nodes[h].attrs["tested"], Attr::Text("test".into()));
        assert!(!t.nodes[h].attrs.contains_key("Rate Distributions"));
        let n2 = t.nodes[h].parent.unwrap();
        assert_eq!(t.nodes[n2].label, None);
        assert_eq!(t.nodes[n2].attrs["hyphy_node"], Attr::Text("Node2".into()));
        assert_eq!(t.nodes[n2].attrs["Corrected P-value"], Attr::Num(0.4));
        assert!(!is_hyphy_json("{\"a\": 1}"));
    }
}
