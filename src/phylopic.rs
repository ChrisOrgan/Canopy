//! PhyloPic (https://www.phylopic.org/) silhouette lookup.
//!
//! Requests run on background threads; results are cached on disk
//! (including misses) so each taxon is only looked up once. Attribution
//! and license information is kept for the credits dialog, since many
//! silhouettes are CC BY and must be credited in publications.

use crate::render::ImageMap;
use anyhow::{anyhow, Context, Result};
use image::RgbaImage;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::io::Read;
use std::path::PathBuf;
use std::sync::mpsc::{channel, Receiver, Sender};
use std::sync::{Arc, Mutex, OnceLock};
use std::time::Duration;

const API: &str = "https://api.phylopic.org";
const USER_AGENT: &str = concat!("Canopy/", env!("CARGO_PKG_VERSION"), " (phylogenetic tree viewer)");

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PicMeta {
    pub query: String,
    pub matched_name: String,
    pub attribution: Option<String>,
    pub license: Option<String>,
    pub image_uuid: String,
    pub image_url: String,
}

impl PicMeta {
    pub fn page_url(&self) -> String {
        format!("https://www.phylopic.org/images/{}", self.image_uuid)
    }

    pub fn license_name(&self) -> String {
        let l = self.license.as_deref().unwrap_or("");
        if l.contains("publicdomain/zero") {
            "CC0 1.0".into()
        } else if l.contains("publicdomain/mark") {
            "Public Domain Mark".into()
        } else if l.contains("by-nc-sa") {
            "CC BY-NC-SA".into()
        } else if l.contains("by-nc") {
            "CC BY-NC".into()
        } else if l.contains("by-sa") {
            "CC BY-SA".into()
        } else if l.contains("/by/") {
            "CC BY".into()
        } else if l.is_empty() {
            "unknown".into()
        } else {
            l.to_string()
        }
    }
}

#[derive(Clone)]
pub enum PicStatus {
    Pending,
    Missing,
    Failed(String),
    Ready(Arc<RgbaImage>, Option<PicMeta>),
}

type Outcome = (String, Result<Option<(RgbaImage, PicMeta)>, String>);

pub struct PhyloPic {
    tx: Sender<String>,
    rx: Receiver<Outcome>,
    pub entries: HashMap<String, PicStatus>,
    pub enabled: bool,
}

impl PhyloPic {
    /// `repaint` is called from worker threads when a result arrives.
    pub fn new(repaint: Arc<dyn Fn() + Send + Sync>) -> Self {
        let (tx, job_rx) = channel::<String>();
        let (res_tx, rx) = channel::<Outcome>();
        let job_rx = Arc::new(Mutex::new(job_rx));
        for _ in 0..3 {
            let job_rx = job_rx.clone();
            let res_tx = res_tx.clone();
            let repaint = repaint.clone();
            std::thread::spawn(move || loop {
                let job = match job_rx.lock() {
                    Ok(rx) => rx.recv(),
                    Err(_) => return,
                };
                let Ok(key) = job else { return };
                let r = fetch(&key).map_err(|e| format!("{:#}", e));
                if res_tx.send((key, r)).is_err() {
                    return;
                }
                repaint();
            });
        }
        PhyloPic { tx, rx, entries: HashMap::new(), enabled: true }
    }

    pub fn request(&mut self, key: &str) {
        if !self.enabled || key.is_empty() || self.entries.contains_key(key) {
            return;
        }
        if query_names(key).is_empty() {
            self.entries.insert(key.to_string(), PicStatus::Missing);
            return;
        }
        self.entries.insert(key.to_string(), PicStatus::Pending);
        let _ = self.tx.send(key.to_string());
    }

    /// Forget a result (e.g. to retry after a network failure).
    pub fn retry_failed(&mut self) {
        let failed: Vec<String> =
            self.entries.iter().filter(|(_, s)| matches!(s, PicStatus::Failed(_))).map(|(k, _)| k.clone()).collect();
        for k in failed {
            self.entries.remove(&k);
            self.request(&k);
        }
    }

    /// Drain finished jobs; returns true if anything changed.
    pub fn poll(&mut self) -> bool {
        let mut changed = false;
        while let Ok((key, r)) = self.rx.try_recv() {
            let st = match r {
                Ok(Some((img, meta))) => PicStatus::Ready(Arc::new(img), Some(meta)),
                Ok(None) => PicStatus::Missing,
                Err(e) => PicStatus::Failed(e),
            };
            self.entries.insert(key, st);
            changed = true;
        }
        changed
    }

    pub fn set_custom(&mut self, key: &str, img: RgbaImage) {
        self.entries.insert(key.to_string(), PicStatus::Ready(Arc::new(img), None));
    }

    pub fn images(&self) -> ImageMap {
        self.entries
            .iter()
            .filter_map(|(k, s)| match s {
                PicStatus::Ready(img, _) => Some((k.clone(), img.clone())),
                _ => None,
            })
            .collect()
    }

    pub fn credits(&self, keys: &[String]) -> Vec<PicMeta> {
        let mut v: Vec<PicMeta> = keys
            .iter()
            .filter_map(|k| match self.entries.get(k) {
                Some(PicStatus::Ready(_, Some(m))) => Some(m.clone()),
                _ => None,
            })
            .collect();
        v.sort_by(|a, b| a.query.cmp(&b.query));
        v.dedup_by(|a, b| a.image_uuid == b.image_uuid);
        v
    }

    pub fn counts(&self) -> (usize, usize, usize, usize) {
        let mut c = (0, 0, 0, 0);
        for s in self.entries.values() {
            match s {
                PicStatus::Ready(..) => c.0 += 1,
                PicStatus::Pending => c.1 += 1,
                PicStatus::Missing => c.2 += 1,
                PicStatus::Failed(_) => c.3 += 1,
            }
        }
        c
    }
}

/// Names to try for a label: the binomial (first two alphabetic words),
/// then the genus.
pub fn query_names(key: &str) -> Vec<String> {
    let words: Vec<String> = key
        .split(|c: char| c.is_whitespace() || c == '_')
        .filter(|w| !w.is_empty())
        .take(2)
        .take_while(|w| w.chars().all(|c| c.is_alphabetic() || c == '-'))
        .map(|w| w.to_lowercase())
        .collect();
    let mut out = Vec::new();
    if words.len() >= 2 {
        out.push(format!("{} {}", words[0], words[1]));
    }
    if let Some(g) = words.first() {
        if g.len() > 1 {
            out.push(g.clone());
        }
    }
    out
}

fn cache_dir() -> Option<PathBuf> {
    let d = dirs::cache_dir()?.join("Canopy").join("phylopic");
    std::fs::create_dir_all(&d).ok()?;
    Some(d)
}

fn cache_stem(key: &str) -> String {
    key.chars().map(|c| if c.is_ascii_alphanumeric() { c.to_ascii_lowercase() } else { '_' }).collect()
}

fn get(url: &str) -> Result<ureq::Response> {
    ureq::get(url)
        .set("User-Agent", USER_AGENT)
        .set("Accept", "application/vnd.phylopic.v2+json")
        .timeout(Duration::from_secs(25))
        .call()
        .map_err(anyhow::Error::from)
}

fn is_not_found(e: &anyhow::Error) -> bool {
    matches!(e.downcast_ref::<ureq::Error>(), Some(ureq::Error::Status(404, _)))
}

fn current_build() -> Result<u64> {
    static BUILD: OnceLock<u64> = OnceLock::new();
    if let Some(b) = BUILD.get() {
        return Ok(*b);
    }
    let v: serde_json::Value = serde_json::from_str(&get(&format!("{}/", API))?.into_string()?)?;
    let b = v["build"].as_u64().ok_or_else(|| anyhow!("PhyloPic API returned no build number"))?;
    let _ = BUILD.set(b);
    Ok(b)
}

fn urlencode(s: &str) -> String {
    let mut o = String::new();
    for b in s.bytes() {
        match b {
            b'A'..=b'Z' | b'a'..=b'z' | b'0'..=b'9' | b'-' | b'_' | b'.' | b'~' => o.push(b as char),
            _ => o.push_str(&format!("%{:02X}", b)),
        }
    }
    o
}

fn first_name(item: &serde_json::Value) -> String {
    item["names"][0][0]["text"].as_str().unwrap_or("").to_string()
}

/// Search PhyloPic for one name; returns (matched name, primary image JSON).
fn search(name: &str, build: u64) -> Result<Option<(String, serde_json::Value)>> {
    let url = format!(
        "{}/nodes?build={}&filter_name={}&page=0&embed_items=true&embed_primaryImage=true",
        API,
        build,
        urlencode(name)
    );
    // The API answers 404 when nothing matches the name.
    let resp = match get(&url) {
        Ok(r) => r,
        Err(e) if is_not_found(&e) => return Ok(None),
        Err(e) => return Err(e),
    };
    let v: serde_json::Value = serde_json::from_str(&resp.into_string()?)?;
    let Some(items) = v["_embedded"]["items"].as_array() else { return Ok(None) };
    let with_image: Vec<&serde_json::Value> = items.iter().filter(|i| i["_embedded"]["primaryImage"].is_object()).collect();
    let exact = with_image.iter().find(|i| first_name(i).to_lowercase() == name);
    let pick = exact.or_else(|| with_image.first());
    Ok(pick.map(|i| (first_name(i), i["_embedded"]["primaryImage"].clone())))
}

/// Look up one taxon synchronously (uses the disk cache).
pub fn fetch_blocking(key: &str) -> Result<Option<(RgbaImage, PicMeta)>> {
    fetch(key)
}

fn fetch(key: &str) -> Result<Option<(RgbaImage, PicMeta)>> {
    let dir = cache_dir();
    let stem = cache_stem(key);
    if let Some(d) = &dir {
        let png = d.join(format!("{}.png", stem));
        let json = d.join(format!("{}.json", stem));
        let none = d.join(format!("{}.none", stem));
        if png.exists() && json.exists() {
            let img = image::open(&png)?.to_rgba8();
            let meta: PicMeta = serde_json::from_str(&std::fs::read_to_string(&json)?)?;
            return Ok(Some((img, meta)));
        }
        if let Ok(m) = std::fs::metadata(&none) {
            let fresh = m.modified().ok().and_then(|t| t.elapsed().ok()).map(|e| e < Duration::from_secs(30 * 86400)).unwrap_or(true);
            if fresh {
                return Ok(None);
            }
        }
    }
    let build = current_build()?;
    for name in query_names(key) {
        let Some((matched, img)) = search(&name, build)? else { continue };
        let files = img["_links"]["rasterFiles"].as_array().cloned().unwrap_or_default();
        let width = |f: &serde_json::Value| {
            f["sizes"].as_str().and_then(|s| s.split('x').next()).and_then(|w| w.parse::<u32>().ok()).unwrap_or(0)
        };
        // Smallest raster at least 1024 px wide (sharp at 600 dpi for typical sizes), else the largest.
        let mut sorted = files.clone();
        sorted.sort_by_key(width);
        let file = sorted.iter().find(|f| width(f) >= 1024).or_else(|| sorted.last());
        let Some(href) = file.and_then(|f| f["href"].as_str()) else { continue };
        let mut bytes = Vec::new();
        get(href)?.into_reader().take(20_000_000).read_to_end(&mut bytes)?;
        let image = image::load_from_memory(&bytes).context("decoding silhouette")?.to_rgba8();
        let meta = PicMeta {
            query: key.to_string(),
            matched_name: matched,
            attribution: img["attribution"].as_str().map(String::from),
            license: img["_links"]["license"]["href"].as_str().map(String::from),
            image_uuid: img["uuid"].as_str().unwrap_or("").to_string(),
            image_url: href.to_string(),
        };
        if let Some(d) = &dir {
            let _ = image.save(d.join(format!("{}.png", stem)));
            let _ = std::fs::write(d.join(format!("{}.json", stem)), serde_json::to_string_pretty(&meta)?);
        }
        return Ok(Some((image, meta)));
    }
    if let Some(d) = &dir {
        let _ = std::fs::write(d.join(format!("{}.none", stem)), b"");
    }
    Ok(None)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names() {
        assert_eq!(query_names("Homo sapiens"), vec!["homo sapiens", "homo"]);
        assert_eq!(query_names("Canis_lupus_familiaris"), vec!["canis lupus", "canis"]);
        assert_eq!(query_names("Gorilla"), vec!["gorilla"]);
        assert!(query_names("seq_00123").len() == 1);
        assert!(query_names("12345").is_empty());
    }
}
