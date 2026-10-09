//! NEXUS reader/writer (TAXA and TREES blocks, TRANSLATE tables).

use super::newick::{parse_newick, quote_label, write_newick_with, WriteOptions};
use crate::tree::Tree;
use anyhow::{bail, Context, Result};
use std::collections::HashMap;

pub struct NexusFile {
    pub taxa: Vec<String>,
    pub trees: Vec<Tree>,
}

/// Split NEXUS text into `;`-terminated commands, respecting quotes and
/// bracketed comments (which are preserved, since `[&R]` and `[&...]`
/// annotations matter inside TREE commands).
fn commands(s: &str) -> Vec<&str> {
    let b = s.as_bytes();
    let mut out = Vec::new();
    let mut start = 0;
    let mut depth = 0;
    let mut quote: Option<u8> = None;
    for (i, &c) in b.iter().enumerate() {
        match (quote, c) {
            (Some(q), _) if c == q => quote = None,
            (Some(_), _) => {}
            (None, b'\'') | (None, b'"') if depth == 0 => quote = Some(c),
            (None, b'[') => depth += 1,
            (None, b']') => depth -= 1,
            (None, b';') if depth == 0 => {
                out.push(&s[start..i]);
                start = i + 1;
            }
            _ => {}
        }
    }
    if start < s.len() && !s[start..].trim().is_empty() {
        out.push(&s[start..]);
    }
    out
}

/// Remove leading whitespace and plain (non-`&`) comments.
fn strip_leading_comments(mut s: &str) -> &str {
    loop {
        s = s.trim_start();
        if s.starts_with('[') && !s.starts_with("[&") {
            match s.find(']') {
                Some(e) => s = &s[e + 1..],
                None => return "",
            }
        } else {
            return s;
        }
    }
}

/// Remove every plain comment (keeping `[&...]` ones).
fn strip_plain_comments(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut depth = 0;
    let mut keep = false;
    let chars: Vec<char> = s.chars().collect();
    for (i, &c) in chars.iter().enumerate() {
        if c == '[' {
            if depth == 0 {
                keep = chars.get(i + 1) == Some(&'&');
            }
            depth += 1;
            if keep {
                out.push(c);
            }
        } else if c == ']' && depth > 0 {
            depth -= 1;
            if keep {
                out.push(c);
            }
        } else if depth == 0 || keep {
            out.push(c);
        }
    }
    out
}

fn split_word(s: &str) -> (&str, &str) {
    let s = s.trim_start();
    match s.find(|c: char| c.is_whitespace()) {
        Some(i) => (&s[..i], &s[i..]),
        None => (s, ""),
    }
}

fn unquote(s: &str) -> String {
    let s = s.trim();
    if s.len() >= 2 && ((s.starts_with('\'') && s.ends_with('\'')) || (s.starts_with('"') && s.ends_with('"'))) {
        s[1..s.len() - 1].replace("''", "'")
    } else {
        s.to_string()
    }
}

/// Split on whitespace outside quotes.
fn tokens(s: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut quote: Option<char> = None;
    for c in s.chars() {
        match quote {
            Some(q) if c == q => {
                quote = None;
                cur.push(c);
            }
            Some(_) => cur.push(c),
            None if c == '\'' || c == '"' => {
                quote = Some(c);
                cur.push(c);
            }
            None if c.is_whitespace() || c == ',' => {
                if !cur.is_empty() {
                    out.push(unquote(&cur));
                    cur.clear();
                }
            }
            None => cur.push(c),
        }
    }
    if !cur.is_empty() {
        out.push(unquote(&cur));
    }
    out
}

pub fn is_nexus(s: &str) -> bool {
    s.trim_start().get(..6).map(|h| h.eq_ignore_ascii_case("#nexus")).unwrap_or(false)
}

pub fn parse_nexus(s: &str) -> Result<NexusFile> {
    let body = s.trim_start();
    let body = if is_nexus(body) { &body[6..] } else { body };
    let mut block = String::new();
    let mut taxa = Vec::new();
    let mut translate: HashMap<String, String> = HashMap::new();
    let mut trees = Vec::new();

    for cmd in commands(body) {
        let cmd = strip_leading_comments(cmd);
        if cmd.is_empty() {
            continue;
        }
        let (word, rest) = split_word(cmd);
        let word = word.to_ascii_lowercase();
        match word.as_str() {
            "begin" => block = rest.trim().to_ascii_lowercase(),
            "end" | "endblock" => block.clear(),
            "taxlabels" if block == "taxa" || block == "data" => {
                taxa = tokens(&strip_plain_comments(rest));
            }
            "translate" if block == "trees" => {
                let rest = strip_plain_comments(rest);
                let toks = tokens(&rest);
                for pair in toks.chunks(2) {
                    if let [k, v] = pair {
                        translate.insert(k.clone(), v.clone());
                    }
                }
            }
            "tree" | "utree" if block == "trees" => {
                // BEAST writes `tree STATE_0 [&lnP=-1.0] = [&R] (...)`: find '=' outside comments.
                let mut depth = 0;
                let eq = rest.char_indices().find_map(|(i, c)| {
                    match c {
                        '[' => depth += 1,
                        ']' => depth -= 1,
                        '=' if depth == 0 => return Some(i),
                        _ => {}
                    }
                    None
                });
                let Some(eq) = eq else { bail!("TREE command without '='") };
                let raw_name = strip_plain_comments(&rest[..eq]);
                let raw_name = raw_name.split('[').next().unwrap_or("");
                let name = unquote(raw_name.trim().trim_start_matches('*').trim());
                let newick = format!("{};", &rest[eq + 1..]);
                let mut t = parse_newick(&newick).with_context(|| format!("in tree '{}'", name))?;
                t.name = Some(name);
                if word == "utree" && t.rooted.is_none() {
                    t.rooted = Some(false);
                }
                if !translate.is_empty() {
                    for n in t.tips() {
                        if let Some(l) = &t.nodes[n].label {
                            if let Some(real) = translate.get(l) {
                                t.nodes[n].label = Some(real.clone());
                            }
                        }
                    }
                }
                trees.push(t);
            }
            _ => {}
        }
    }
    if trees.is_empty() {
        bail!("no TREES block / TREE commands found in NEXUS file");
    }
    Ok(NexusFile { taxa, trees })
}

/// Write trees as NEXUS with a TAXA block and a TRANSLATE table.
pub fn write_nexus(trees: &[&Tree], opts: &WriteOptions) -> String {
    let mut taxa: Vec<String> = Vec::new();
    let mut index: HashMap<String, usize> = HashMap::new();
    for t in trees {
        for n in t.tips() {
            let l = t.label(n).to_string();
            if !index.contains_key(&l) {
                index.insert(l.clone(), taxa.len() + 1);
                taxa.push(l);
            }
        }
    }
    let mut s = String::from("#NEXUS\n\nBEGIN TAXA;\n");
    s += &format!("\tDIMENSIONS NTAX={};\n\tTAXLABELS\n", taxa.len());
    for t in &taxa {
        s += &format!("\t\t{}\n", quote_label(t));
    }
    s += "\t;\nEND;\n\nBEGIN TREES;\n\tTRANSLATE\n";
    for (i, t) in taxa.iter().enumerate() {
        let sep = if i + 1 == taxa.len() { "" } else { "," };
        s += &format!("\t\t{} {}{}\n", i + 1, quote_label(t), sep);
    }
    s += "\t;\n";
    for (i, t) in trees.iter().enumerate() {
        let name = t.name.clone().unwrap_or_else(|| format!("tree_{}", i + 1));
        let rooted = match t.rooted {
            Some(false) => "[&U] ",
            _ => "[&R] ",
        };
        let nwk = write_newick_with(t, opts, &|n| {
            if t.is_tip(n) {
                t.nodes[n].label.as_ref().and_then(|l| index.get(l)).map(|i| i.to_string())
            } else {
                None
            }
        });
        s += &format!("\tTREE {} = {}{}\n", quote_label(&name), rooted, nwk);
    }
    s += "END;\n";
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    const MRBAYES: &str = "#NEXUS\n[comment; with semicolon]\nbegin trees;\n  translate\n    1 Homo_sapiens,\n    2 'Pan troglodytes',\n    3 Gorilla;\n  tree gen.1 = [&U] ((1:0.1,2:0.2)[&prob=0.9]:0.05,3:0.3);\n  tree gen.2 = [&U] ((1:0.1,3:0.2):0.05,2:0.3);\nend;\n";

    #[test]
    fn translate_and_comments() {
        let f = parse_nexus(MRBAYES).unwrap();
        assert_eq!(f.trees.len(), 2);
        let t = &f.trees[0];
        assert_eq!(t.rooted, Some(false));
        assert!(t.find_label("Pan troglodytes").is_some());
        assert!(t.find_label("Homo_sapiens").is_some());
        assert_eq!(t.name.as_deref(), Some("gen.1"));
    }

    #[test]
    fn roundtrip() {
        let f = parse_nexus(MRBAYES).unwrap();
        let refs: Vec<&Tree> = f.trees.iter().collect();
        let out = write_nexus(&refs, &WriteOptions { annotations: true, ..Default::default() });
        let g = parse_nexus(&out).unwrap();
        assert_eq!(g.trees.len(), 2);
        assert!(g.trees[1].find_label("Pan troglodytes").is_some());
        assert_eq!(g.taxa.len(), 3);
    }
}
