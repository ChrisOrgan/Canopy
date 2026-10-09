//! phyloXML (<http://www.phyloxml.org>, schema 1.20): reading and writing.
//!
//! Clade `name` (or the taxonomy's scientific name) becomes the label,
//! `branch_length` (element or attribute) the length, `confidence` a support
//! attribute (type "bootstrap" -> `bootstrap`, "probability" -> `posterior`,
//! other types keep their name), and taxonomy fields, `color`, `date` and
//! `property` elements become attributes. Canopy's other attributes are
//! written as `canopy:` properties, so a save and reload keeps them.

use crate::io::newick::WriteOptions;
use crate::tree::{Attr, NodeId, Tree};
use anyhow::{bail, Context, Result};
use std::fmt::Write as _;

/// True if `text` looks like a phyloXML document.
pub fn is_phyloxml(text: &str) -> bool {
    let t = text.trim_start_matches('\u{feff}').trim_start();
    t.starts_with('<') && (t.contains("<phyloxml") || t.contains("<phylogeny"))
}

/// Every `<phylogeny>` in the document.
pub fn parse_phyloxml(text: &str) -> Result<Vec<Tree>> {
    let doc = roxmltree::Document::parse(text.trim_start_matches('\u{feff}')).context("invalid XML")?;
    let mut trees = Vec::new();
    for ph in doc.descendants().filter(|n| n.has_tag_name_local("phylogeny")) {
        let mut tree = Tree::new();
        tree.rooted = match ph.attribute("rooted") {
            Some("true") => Some(true),
            Some("false") => Some(false),
            _ => None,
        };
        tree.name = child(ph, "name").and_then(|n| n.text()).map(|s| s.trim().to_string()).filter(|s| !s.is_empty());
        let Some(clade) = child(ph, "clade") else { continue };
        let root = tree.root;
        read_clade(clade, &mut tree, root)?;
        trees.push(tree);
    }
    if trees.is_empty() {
        bail!("no <phylogeny> with a <clade> found");
    }
    Ok(trees)
}

trait LocalName {
    fn has_tag_name_local(&self, name: &str) -> bool;
}

impl LocalName for roxmltree::Node<'_, '_> {
    fn has_tag_name_local(&self, name: &str) -> bool {
        self.is_element() && self.tag_name().name() == name
    }
}

fn child<'a, 'i>(n: roxmltree::Node<'a, 'i>, name: &str) -> Option<roxmltree::Node<'a, 'i>> {
    n.children().find(|c| c.has_tag_name_local(name))
}

fn text_of(n: roxmltree::Node) -> String {
    n.text().unwrap_or("").trim().to_string()
}

fn value(s: &str) -> Attr {
    match s.trim().parse::<f64>() {
        Ok(v) => Attr::Num(v),
        Err(_) => Attr::Text(s.trim().to_string()),
    }
}

fn read_clade(el: roxmltree::Node, tree: &mut Tree, id: NodeId) -> Result<()> {
    if let Some(bl) = el.attribute("branch_length") {
        tree.nodes[id].length = Some(bl.trim().parse().with_context(|| format!("bad branch_length '{}'", bl))?);
    }
    for c in el.children().filter(|c| c.is_element()) {
        let attrs = &mut tree.nodes[id].attrs;
        match c.tag_name().name() {
            "name" => {
                let s = text_of(c);
                if !s.is_empty() {
                    tree.nodes[id].label = Some(s);
                }
            }
            "branch_length" => {
                let s = text_of(c);
                tree.nodes[id].length = Some(s.parse().with_context(|| format!("bad branch_length '{}'", s))?);
            }
            "confidence" => {
                let key = match c.attribute("type").unwrap_or("").to_ascii_lowercase().as_str() {
                    "bootstrap" => "bootstrap".to_string(),
                    "probability" | "posterior" | "bayesian" | "bayesian_posterior" => "posterior".to_string(),
                    "" => "support".to_string(),
                    other => other.to_string(),
                };
                attrs.insert(key, value(&text_of(c)));
            }
            "taxonomy" => {
                for t in c.children().filter(|t| t.is_element()) {
                    let key = match t.tag_name().name() {
                        "scientific_name" => "scientific_name",
                        "common_name" => "common_name",
                        "rank" => "rank",
                        "code" => "taxonomy_code",
                        "id" => "taxonomy_id",
                        _ => continue,
                    };
                    let s = text_of(t);
                    if !s.is_empty() {
                        attrs.insert(key.into(), Attr::Text(s));
                    }
                }
            }
            "color" => {
                let ch = |n: &str| child(c, n).map(text_of).and_then(|s| s.parse::<u8>().ok()).unwrap_or(0);
                attrs.insert("color".into(), Attr::Text(format!("#{:02X}{:02X}{:02X}", ch("red"), ch("green"), ch("blue"))));
            }
            "date" => {
                if let Some(v) = child(c, "value") {
                    attrs.insert("date".into(), value(&text_of(v)));
                }
            }
            "property" => {
                let Some(r) = c.attribute("ref") else { continue };
                let s = text_of(c);
                let attr = match r.strip_prefix("canopy:") {
                    // Canopy's own: ranges and lists were written comma-separated.
                    Some(k) => {
                        let parts: Vec<&str> = s.split(',').collect();
                        let nums: Option<Vec<f64>> = parts.iter().map(|p| p.trim().parse().ok()).collect();
                        let a = match nums {
                            Some(v) if v.len() == 2 && c.attribute("datatype") == Some("xsd:string") => Attr::Range(v[0], v[1]),
                            Some(v) if v.len() > 2 => Attr::List(v.into_iter().map(Attr::Num).collect()),
                            _ => value(&s),
                        };
                        tree.nodes[id].attrs.insert(decode_key(k), a);
                        continue;
                    }
                    None => value(&s),
                };
                attrs.insert(r.to_string(), attr);
            }
            "clade" => {
                let k = tree.add_node(Some(id));
                read_clade(c, tree, k)?;
            }
            _ => {}
        }
    }
    // Tips named only through their taxonomy.
    if tree.nodes[id].label.is_none() && tree.nodes[id].children.is_empty() {
        if let Some(Attr::Text(s)) = tree.nodes[id].attrs.get("scientific_name") {
            tree.nodes[id].label = Some(s.clone());
        }
    }
    Ok(())
}

/// phyloXML property names must match `[a-zA-Z0-9_]+`: other characters are
/// written as `_xHH_` (UTF-8 bytes).
fn encode_key(k: &str) -> String {
    let mut out = String::new();
    for b in k.bytes() {
        if b.is_ascii_alphanumeric() || b == b'_' {
            out.push(b as char);
        } else {
            let _ = write!(out, "_x{:02X}_", b);
        }
    }
    out
}

fn decode_key(k: &str) -> String {
    let b = k.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    while i < b.len() {
        if b[i] == b'_' && i + 4 < b.len() && b[i + 1] == b'x' && b[i + 4] == b'_' {
            if let Ok(v) = u8::from_str_radix(&k[i + 2..i + 4], 16) {
                out.push(v);
                i += 5;
                continue;
            }
        }
        out.push(b[i]);
        i += 1;
    }
    String::from_utf8_lossy(&out).into_owned()
}

fn esc(s: &str) -> String {
    s.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;").replace('"', "&quot;").replace('\'', "&apos;")
}

/// Write trees as one phyloXML document. `opts` decides whether branch
/// lengths, internal node names and annotations (confidence, taxonomy, color,
/// date, properties) are included.
pub fn write_phyloxml(trees: &[&Tree], opts: &WriteOptions) -> String {
    let mut s = String::from(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<phyloxml xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" \
         xsi:schemaLocation=\"http://www.phyloxml.org http://www.phyloxml.org/1.20/phyloxml.xsd\" xmlns=\"http://www.phyloxml.org\">\n",
    );
    for t in trees {
        let rooted = t.rooted != Some(false);
        let _ = writeln!(s, "  <phylogeny rooted=\"{}\">", rooted);
        if let Some(n) = &t.name {
            let _ = writeln!(s, "    <name>{}</name>", esc(n));
        }
        write_clade(t, t.root, 2, opts, &mut s);
        s.push_str("  </phylogeny>\n");
    }
    s.push_str("</phyloxml>\n");
    s
}

/// Attributes written as phyloXML elements rather than `canopy:` properties.
const NATIVE: &[&str] = &["scientific_name", "common_name", "rank", "taxonomy_code", "taxonomy_id", "color", "date"];

fn write_clade(t: &Tree, n: NodeId, depth: usize, opts: &WriteOptions, s: &mut String) {
    let pad = "  ".repeat(depth);
    let node = &t.nodes[n];
    let _ = writeln!(s, "{}<clade>", pad);
    let inner = "  ".repeat(depth + 1);
    let support = crate::ops::SUPPORT_KEYS;
    // A numeric internal label is the support value, written as confidence.
    let label_is_support = !t.is_tip(n) && node.label.as_deref().is_some_and(|l| l.trim().parse::<f64>().is_ok()) && node.attrs.contains_key("support");
    let show_name = t.is_tip(n) || (opts.internal_labels && !(label_is_support && opts.annotations));
    if let Some(l) = node.label.as_deref().filter(|_| show_name) {
        let _ = writeln!(s, "{}<name>{}</name>", inner, esc(l));
    }
    if let Some(len) = node.length.filter(|_| opts.lengths) {
        let _ = writeln!(s, "{}<branch_length>{}</branch_length>", inner, len);
    }
    if opts.annotations {
        for k in support {
            if let Some(v) = node.attrs.get(k).and_then(Attr::as_f64) {
                let ty = if k == "posterior" { "probability" } else { k };
                let _ = writeln!(s, "{}<confidence type=\"{}\">{}</confidence>", inner, ty, v);
            }
        }
        if let Some(Attr::Text(c)) = node.attrs.get("color") {
            let h = c.trim_start_matches('#');
            let p = |i: usize| u8::from_str_radix(h.get(i..i + 2).unwrap_or("00"), 16).unwrap_or(0);
            let _ = writeln!(s, "{}<color><red>{}</red><green>{}</green><blue>{}</blue></color>", inner, p(0), p(2), p(4));
        }
        let tax: Vec<(&str, &str)> = [("taxonomy_id", "id"), ("taxonomy_code", "code"), ("scientific_name", "scientific_name"), ("common_name", "common_name"), ("rank", "rank")]
            .into_iter()
            .filter(|(k, _)| node.attrs.contains_key(*k))
            .collect();
        if !tax.is_empty() {
            let _ = writeln!(s, "{}<taxonomy>", inner);
            for (k, el) in tax {
                let _ = writeln!(s, "{}  <{}>{}</{}>", inner, el, esc(&node.attrs[k].display(12)), el);
            }
            let _ = writeln!(s, "{}</taxonomy>", inner);
        }
        if let Some(d) = node.attrs.get("date") {
            let _ = writeln!(s, "{}<date><value>{}</value></date>", inner, esc(&d.display(12)));
        }
        for (k, v) in &node.attrs {
            if support.contains(&k.as_str()) || NATIVE.contains(&k.as_str()) {
                continue;
            }
            let (datatype, text) = match v {
                Attr::Num(x) => ("xsd:double", x.to_string()),
                Attr::Text(t) => ("xsd:string", t.clone()),
                Attr::Range(a, b) => ("xsd:string", format!("{},{}", a, b)),
                Attr::List(items) => ("xsd:string", items.iter().map(|a| a.display(12)).collect::<Vec<_>>().join(",")),
            };
            let _ = writeln!(s, "{}<property ref=\"canopy:{}\" datatype=\"{}\" applies_to=\"clade\">{}</property>", inner, encode_key(k), datatype, esc(&text));
        }
    }
    for &c in &node.children {
        write_clade(t, c, depth + 1, opts, s);
    }
    let _ = writeln!(s, "{}</clade>", pad);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::newick::{parse_newick, write_newick, WriteOptions};

    #[test]
    fn reads_phyloxml_example() {
        let x = r#"<?xml version="1.0" encoding="UTF-8"?>
<phyloxml xmlns="http://www.phyloxml.org">
  <phylogeny rooted="true">
    <name>example</name>
    <clade>
      <clade branch_length="0.06">
        <confidence type="bootstrap">89</confidence>
        <clade><name>A</name><branch_length>0.102</branch_length></clade>
        <clade><taxonomy><scientific_name>Homo sapiens</scientific_name></taxonomy><branch_length>0.23</branch_length></clade>
      </clade>
      <clade><name>C</name><branch_length>0.4</branch_length>
        <property ref="NOAA:depth" datatype="xsd:double" applies_to="clade">1200</property></clade>
    </clade>
  </phylogeny>
</phyloxml>"#;
        assert!(is_phyloxml(x));
        let t = parse_phyloxml(x).unwrap().remove(0);
        assert_eq!(t.name.as_deref(), Some("example"));
        assert_eq!(t.rooted, Some(true));
        assert_eq!(t.num_tips(), 3);
        let a = t.find_label("A").unwrap();
        assert_eq!(t.nodes[a].length, Some(0.102));
        let p = t.nodes[a].parent.unwrap();
        assert_eq!(t.nodes[p].attrs["bootstrap"], Attr::Num(89.0));
        assert_eq!(t.nodes[p].length, Some(0.06));
        assert!(t.find_label("Homo sapiens").is_some());
        let c = t.find_label("C").unwrap();
        assert_eq!(t.nodes[c].attrs["NOAA:depth"], Attr::Num(1200.0));
    }

    #[test]
    fn round_trip_keeps_annotations() {
        let t = parse_newick("((A:1,B:2)[&posterior=0.95,height_95%_HPD={0.5,1.5},rate=0.01]:1,(C<x>:3,D:3)88:1);").unwrap();
        let x = write_phyloxml(&[&t], &WriteOptions { annotations: true, ..Default::default() });
        let back = parse_phyloxml(&x).unwrap().remove(0);
        let plain = WriteOptions::default();
        assert_eq!(write_newick(&back, &plain), "((A:1,B:2):1,(C<x>:3,D:3):1);");
        let a = back.find_label("A").unwrap();
        let ab = back.nodes[a].parent.unwrap();
        assert_eq!(back.nodes[ab].attrs["height_95%_HPD"], Attr::Range(0.5, 1.5));
        assert_eq!(back.nodes[ab].attrs["posterior"], Attr::Num(0.95));
        assert_eq!(back.nodes[ab].attrs["rate"], Attr::Num(0.01));
        let cd = back.nodes[back.find_label("D").unwrap()].parent.unwrap();
        assert_eq!(back.nodes[cd].attrs["support"], Attr::Num(88.0));
        assert_eq!(decode_key(&encode_key("height_95%_HPD")), "height_95%_HPD");
    }
}
