//! Newick reader/writer, including BEAST/MrBayes `[&key=value]` and NHX
//! `[&&NHX:key=value]` annotations.

use crate::tree::{format_num, Attr, NodeId, Tree};
use anyhow::{bail, Result};
use std::collections::BTreeMap;

struct Parser<'a> {
    s: &'a [u8],
    i: usize,
}

impl<'a> Parser<'a> {
    fn peek(&self) -> Option<u8> {
        self.s.get(self.i).copied()
    }

    fn skip_ws(&mut self) {
        while let Some(c) = self.peek() {
            if c.is_ascii_whitespace() {
                self.i += 1;
            } else {
                break;
            }
        }
    }

    fn err<T>(&self, msg: &str) -> Result<T> {
        let start = self.i.saturating_sub(20);
        let end = (self.i + 20).min(self.s.len());
        bail!(
            "Newick parse error at byte {}: {} (near \"{}\")",
            self.i,
            msg,
            String::from_utf8_lossy(&self.s[start..end])
        )
    }

    /// Reads a bracketed comment, returning its inner text.
    fn comment(&mut self) -> Result<String> {
        debug_assert_eq!(self.peek(), Some(b'['));
        self.i += 1;
        let start = self.i;
        let mut depth = 1;
        let mut quote: Option<u8> = None;
        while let Some(c) = self.peek() {
            match (quote, c) {
                (Some(q), _) if c == q => quote = None,
                (Some(_), _) => {}
                (None, b'"') | (None, b'\'') => quote = Some(c),
                (None, b'[') => depth += 1,
                (None, b']') => {
                    depth -= 1;
                    if depth == 0 {
                        let text = String::from_utf8_lossy(&self.s[start..self.i]).into_owned();
                        self.i += 1;
                        return Ok(text);
                    }
                }
                _ => {}
            }
            self.i += 1;
        }
        self.err("unterminated comment")
    }

    fn label(&mut self) -> Result<Option<String>> {
        self.skip_ws();
        match self.peek() {
            Some(b'\'') | Some(b'"') => {
                let q = self.peek().unwrap();
                self.i += 1;
                let mut out = Vec::new();
                loop {
                    match self.peek() {
                        None => return self.err("unterminated quoted label"),
                        Some(c) if c == q => {
                            if self.s.get(self.i + 1) == Some(&q) {
                                out.push(q);
                                self.i += 2;
                            } else {
                                self.i += 1;
                                break;
                            }
                        }
                        Some(c) => {
                            out.push(c);
                            self.i += 1;
                        }
                    }
                }
                Ok(Some(String::from_utf8_lossy(&out).into_owned()))
            }
            _ => {
                let start = self.i;
                while let Some(c) = self.peek() {
                    if matches!(c, b'(' | b')' | b',' | b':' | b';' | b'[') {
                        break;
                    }
                    self.i += 1;
                }
                let s = String::from_utf8_lossy(&self.s[start..self.i]).trim().to_string();
                Ok(if s.is_empty() { None } else { Some(s) })
            }
        }
    }

    fn number(&mut self) -> Result<f64> {
        self.skip_ws();
        let start = self.i;
        while let Some(c) = self.peek() {
            if c.is_ascii_digit() || matches!(c, b'.' | b'-' | b'+' | b'e' | b'E') {
                self.i += 1;
            } else {
                break;
            }
        }
        let txt = std::str::from_utf8(&self.s[start..self.i]).unwrap_or("");
        match txt.parse::<f64>() {
            Ok(v) => Ok(v),
            Err(_) => self.err("invalid branch length"),
        }
    }

    /// Label, comments and branch length following a node.
    fn node_suffix(&mut self, tree: &mut Tree, id: NodeId) -> Result<()> {
        if let Some(l) = self.label()? {
            tree.nodes[id].label = Some(l);
        }
        loop {
            self.skip_ws();
            match self.peek() {
                Some(b'[') => {
                    let c = self.comment()?;
                    parse_annotation(&c, &mut tree.nodes[id].attrs);
                }
                Some(b':') => {
                    self.i += 1;
                    tree.nodes[id].length = Some(self.number()?);
                }
                _ => break,
            }
        }
        Ok(())
    }

    fn tree(&mut self) -> Result<Option<Tree>> {
        let mut tree = Tree::new();
        // Leading comments: [&R], [&U], or anything else.
        loop {
            self.skip_ws();
            if self.peek() == Some(b'[') {
                let c = self.comment()?;
                let c = c.trim();
                if c.eq_ignore_ascii_case("&R") {
                    tree.rooted = Some(true);
                } else if c.eq_ignore_ascii_case("&U") {
                    tree.rooted = Some(false);
                }
            } else {
                break;
            }
        }
        if self.peek().is_none() {
            return Ok(None);
        }
        let mut stack: Vec<NodeId> = Vec::new();
        let mut cur = tree.root;
        loop {
            self.skip_ws();
            match self.peek() {
                Some(b'(') => {
                    self.i += 1;
                    stack.push(cur);
                    cur = tree.add_node(Some(cur));
                }
                Some(b',') => {
                    self.i += 1;
                    let Some(&p) = stack.last() else { return self.err("unexpected ','") };
                    cur = tree.add_node(Some(p));
                }
                Some(b')') => {
                    self.i += 1;
                    let Some(p) = stack.pop() else { return self.err("unbalanced ')'") };
                    cur = p;
                    self.node_suffix(&mut tree, cur)?;
                }
                Some(b';') => {
                    self.i += 1;
                    break;
                }
                None => {
                    if stack.is_empty() {
                        break;
                    }
                    return self.err("unexpected end of input");
                }
                Some(_) => {
                    let before = self.i;
                    self.node_suffix(&mut tree, cur)?;
                    if self.i == before {
                        return self.err("unexpected character");
                    }
                }
            }
        }
        if !stack.is_empty() {
            return self.err("unbalanced '('");
        }
        post_process(&mut tree);
        Ok(Some(tree))
    }
}

/// Numeric internal node labels become a `support` attribute.
fn post_process(tree: &mut Tree) {
    for n in tree.preorder() {
        if !tree.is_tip(n) {
            if let Some(v) = tree.nodes[n].label.as_deref().and_then(|l| l.trim().parse::<f64>().ok()) {
                tree.nodes[n].attrs.entry("support".into()).or_insert(Attr::Num(v));
            }
        }
    }
}

/// Parse `&key=value,...` (BEAST) or `&&NHX:key=value:...` comments into attributes.
pub fn parse_annotation(comment: &str, attrs: &mut BTreeMap<String, Attr>) {
    let c = comment.trim();
    if !c.starts_with('&') {
        return;
    }
    let c = c.trim_start_matches('&');
    if let Some(nhx) = c.strip_prefix("NHX:") {
        for kv in nhx.split(':') {
            if let Some((k, v)) = kv.split_once('=') {
                attrs.insert(k.trim().to_string(), parse_value(v.trim()));
            }
        }
        return;
    }
    for kv in split_top(c, b',') {
        if let Some((k, v)) = kv.split_once('=') {
            let key = k.trim().trim_matches('"').to_string();
            let mut val = parse_value(v.trim());
            if let Attr::List(items) = &val {
                let lk = key.to_ascii_lowercase();
                if items.len() == 2 && (lk.contains("hpd") || lk.contains("range")) {
                    if let (Some(a), Some(b)) = (items.first().and_then(Attr::as_f64), items.get(1).and_then(Attr::as_f64)) {
                        val = Attr::Range(a.min(b), a.max(b));
                    }
                }
            }
            attrs.insert(key, val);
        }
    }
}

fn parse_value(v: &str) -> Attr {
    if v.starts_with('{') && v.ends_with('}') {
        let inner = &v[1..v.len() - 1];
        return Attr::List(split_top(inner, b',').into_iter().map(|s| parse_value(s.trim())).collect());
    }
    let unq = v.trim_matches('"').trim_matches('\'');
    if unq.len() != v.len() {
        return Attr::Text(unq.to_string());
    }
    match v.parse::<f64>() {
        Ok(x) => Attr::Num(x),
        Err(_) => Attr::Text(v.to_string()),
    }
}

/// Split at `sep` outside of braces and quotes.
fn split_top(s: &str, sep: u8) -> Vec<&str> {
    let b = s.as_bytes();
    let mut out = Vec::new();
    let mut depth = 0i32;
    let mut quote: Option<u8> = None;
    let mut start = 0;
    for (i, &c) in b.iter().enumerate() {
        match (quote, c) {
            (Some(q), _) if c == q => quote = None,
            (Some(_), _) => {}
            (None, b'"') | (None, b'\'') => quote = Some(c),
            (None, b'{') => depth += 1,
            (None, b'}') => depth -= 1,
            (None, _) if c == sep && depth == 0 => {
                out.push(&s[start..i]);
                start = i + 1;
            }
            _ => {}
        }
    }
    if start < s.len() {
        out.push(&s[start..]);
    }
    out
}

/// Parse a single Newick tree.
pub fn parse_newick(s: &str) -> Result<Tree> {
    let mut p = Parser { s: s.as_bytes(), i: 0 };
    match p.tree()? {
        Some(t) => Ok(t),
        None => bail!("no tree found"),
    }
}

/// Parse every `;`-terminated tree in the input.
pub fn parse_newick_multi(s: &str) -> Result<Vec<Tree>> {
    let mut p = Parser { s: s.as_bytes(), i: 0 };
    let mut out = Vec::new();
    while let Some(t) = p.tree()? {
        out.push(t);
        p.skip_ws();
        if p.peek().is_none() {
            break;
        }
    }
    if out.is_empty() {
        bail!("no tree found");
    }
    Ok(out)
}

#[derive(Clone, Debug)]
pub struct WriteOptions {
    pub lengths: bool,
    pub internal_labels: bool,
    /// Write node attributes as `[&k=v,...]` comments.
    pub annotations: bool,
    pub digits: usize,
}

impl Default for WriteOptions {
    fn default() -> Self {
        WriteOptions { lengths: true, internal_labels: true, annotations: false, digits: 10 }
    }
}

pub fn quote_label(l: &str) -> String {
    let needs = l.is_empty() || l.chars().any(|c| "()[]':;, \t\n".contains(c));
    if needs {
        format!("'{}'", l.replace('\'', "''"))
    } else {
        l.to_string()
    }
}

fn attr_to_nexus(a: &Attr, digits: usize) -> String {
    match a {
        Attr::Num(v) => format_num(*v, digits),
        Attr::Text(s) => format!("\"{}\"", s.replace('"', "'")),
        Attr::Range(x, y) => format!("{{{},{}}}", format_num(*x, digits), format_num(*y, digits)),
        Attr::List(v) => format!("{{{}}}", v.iter().map(|x| attr_to_nexus(x, digits)).collect::<Vec<_>>().join(",")),
    }
}

/// Write a tree as Newick. `rename` optionally maps tip labels (used for Nexus translate tables).
pub fn write_newick_with(tree: &Tree, opts: &WriteOptions, rename: &dyn Fn(NodeId) -> Option<String>) -> String {
    fn rec(t: &Tree, n: NodeId, o: &WriteOptions, rename: &dyn Fn(NodeId) -> Option<String>, out: &mut String) {
        let node = &t.nodes[n];
        if !node.children.is_empty() {
            out.push('(');
            for (i, &c) in node.children.iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                rec(t, c, o, rename, out);
            }
            out.push(')');
        }
        let is_tip = node.children.is_empty();
        if is_tip || o.internal_labels {
            if let Some(name) = rename(n) {
                out.push_str(&name);
            } else if let Some(l) = &node.label {
                out.push_str(&quote_label(l));
            }
        }
        if o.annotations && !node.attrs.is_empty() {
            let parts: Vec<String> =
                node.attrs.iter().map(|(k, v)| format!("{}={}", k, attr_to_nexus(v, o.digits))).collect();
            out.push_str(&format!("[&{}]", parts.join(",")));
        }
        if o.lengths {
            if let Some(l) = node.length {
                out.push(':');
                out.push_str(&format_num(l, o.digits));
            }
        }
    }
    let mut s = String::new();
    rec(tree, tree.root, opts, rename, &mut s);
    s.push(';');
    s
}

pub fn write_newick(tree: &Tree, opts: &WriteOptions) -> String {
    write_newick_with(tree, opts, &|_| None)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic() {
        let t = parse_newick("((A:0.1,B:0.2)90:0.3,'C d':0.4)root;").unwrap();
        assert_eq!(t.num_tips(), 3);
        let n = t.find_label("C d").unwrap();
        assert_eq!(t.nodes[n].length, Some(0.4));
        let ab = t.mrca(&[t.find_label("A").unwrap(), t.find_label("B").unwrap()]).unwrap();
        assert_eq!(t.nodes[ab].attrs.get("support"), Some(&Attr::Num(90.0)));
        assert_eq!(t.label(t.root), "root");
    }

    #[test]
    fn beast_annotations() {
        let s = "[&R] ((A[&rate=1.5]:1.0,B[&rate=2]:1.0)[&posterior=0.98,height_95%_HPD={0.8,1.3}]:1.0,C:2.0);";
        let t = parse_newick(s).unwrap();
        assert_eq!(t.rooted, Some(true));
        let a = t.find_label("A").unwrap();
        assert_eq!(t.nodes[a].attrs["rate"], Attr::Num(1.5));
        let ab = t.nodes[a].parent.unwrap();
        assert_eq!(t.nodes[ab].attrs["posterior"], Attr::Num(0.98));
        assert_eq!(t.nodes[ab].attrs["height_95%_HPD"].as_range(), Some((0.8, 1.3)));
    }

    #[test]
    fn nhx_and_roundtrip() {
        let t = parse_newick("((A:1[&&NHX:S=human:E=1.1],B:2):3,C:4);").unwrap();
        let a = t.find_label("A").unwrap();
        assert_eq!(t.nodes[a].attrs["S"], Attr::Text("human".into()));
        let s = write_newick(&t, &WriteOptions::default());
        assert_eq!(s, "((A:1,B:2):3,C:4);");
        let t2 = parse_newick(&s).unwrap();
        assert_eq!(t2.num_tips(), 3);
    }

    #[test]
    fn multi() {
        let v = parse_newick_multi("(A,B);\n(A,(B,C));\n").unwrap();
        assert_eq!(v.len(), 2);
        assert_eq!(v[1].num_tips(), 3);
    }
}
