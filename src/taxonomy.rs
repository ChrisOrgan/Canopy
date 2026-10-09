//! Taxonomic name checking against the Open Tree of Life taxonomy (OTT)
//! through its TNRS service (`/v3/tnrs/match_names`), see
//! <https://github.com/OpenTreeOfLife/germinator/wiki/Open-Tree-of-Life-Web-APIs>.

use anyhow::{Context, Result};
use serde_json::Value;
use std::time::Duration;

const MATCH_NAMES: &str = "https://api.opentreeoflife.org/v3/tnrs/match_names";
const USER_AGENT: &str = concat!("Canopy/", env!("CARGO_PKG_VERSION"), " (phylogenetic tree viewer)");
/// Names sent per request.
const BATCH: usize = 250;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NameStatus {
    /// The name is an accepted OTT name.
    Exact,
    /// The name is a synonym; `accepted` holds the current name.
    Synonym,
    /// No exact match; `accepted` is the closest (fuzzy) match.
    Approximate,
    /// Nothing found.
    Unmatched,
}

impl NameStatus {
    pub fn label(self) -> &'static str {
        match self {
            NameStatus::Exact => "ok",
            NameStatus::Synonym => "synonym",
            NameStatus::Approximate => "misspelled?",
            NameStatus::Unmatched => "not found",
        }
    }
}

/// The result for one name.
#[derive(Clone, Debug)]
pub struct NameCheck {
    /// The name as sent (underscores as spaces).
    pub query: String,
    pub status: NameStatus,
    /// Accepted OTT name of the best match.
    pub accepted: Option<String>,
    /// Name that matched (the synonym or the fuzzy hit).
    pub matched: Option<String>,
    pub score: f64,
    pub ott_id: Option<u64>,
    pub rank: Option<String>,
    /// OTT flags such as "extinct".
    pub flags: Vec<String>,
    /// Other accepted names that matched equally well (homonyms).
    pub alternatives: Vec<String>,
}

impl NameCheck {
    /// The accepted name when it differs from the query.
    pub fn suggestion(&self) -> Option<&str> {
        self.accepted.as_deref().filter(|a| *a != self.query)
    }
}

/// The name to send for a tip label: underscores as spaces, trimmed.
pub fn query_name(label: &str) -> String {
    label.replace('_', " ").split_whitespace().collect::<Vec<_>>().join(" ")
}

/// Check names against the Open Tree taxonomy (blocking; one request per 250
/// names). Results are in the order of `names`.
pub fn check_names(names: &[String]) -> Result<Vec<NameCheck>> {
    let mut out = Vec::with_capacity(names.len());
    for chunk in names.chunks(BATCH) {
        let body = serde_json::json!({ "names": chunk, "do_approximate_matching": true });
        let text = ureq::post(MATCH_NAMES)
            .set("User-Agent", USER_AGENT)
            .set("Content-Type", "application/json")
            .timeout(Duration::from_secs(60))
            .send_string(&body.to_string())
            .context("Open Tree of Life request failed")?
            .into_string()?;
        let resp: Value = serde_json::from_str(&text).context("unexpected reply from Open Tree of Life")?;
        out.extend(parse_response(&resp, chunk));
    }
    Ok(out)
}

/// Turn a `match_names` reply into one `NameCheck` per name in `names`.
pub fn parse_response(resp: &Value, names: &[String]) -> Vec<NameCheck> {
    let results = resp["results"].as_array().cloned().unwrap_or_default();
    names
        .iter()
        .map(|q| {
            let matches = results.iter().find(|r| r["name"].as_str() == Some(q.as_str())).and_then(|r| r["matches"].as_array()).cloned().unwrap_or_default();
            check_from_matches(q, &matches)
        })
        .collect()
}

fn check_from_matches(query: &str, matches: &[Value]) -> NameCheck {
    let mut check = NameCheck {
        query: query.to_string(),
        status: NameStatus::Unmatched,
        accepted: None,
        matched: None,
        score: 0.0,
        ott_id: None,
        rank: None,
        flags: Vec::new(),
        alternatives: Vec::new(),
    };
    // Prefer exact over synonym over fuzzy, then the higher score.
    let rank_of = |m: &Value| {
        let approx = m["is_approximate_match"].as_bool().unwrap_or(false);
        let syn = m["is_synonym"].as_bool().unwrap_or(false);
        (if approx { 0 } else if syn { 1 } else { 2 }, m["score"].as_f64().unwrap_or(0.0))
    };
    let Some(best) = matches.iter().max_by(|a, b| rank_of(a).partial_cmp(&rank_of(b)).unwrap_or(std::cmp::Ordering::Equal)) else { return check };
    let taxon = &best["taxon"];
    check.accepted = taxon["name"].as_str().map(str::to_string);
    check.matched = best["matched_name"].as_str().map(str::to_string);
    check.score = best["score"].as_f64().unwrap_or(0.0);
    check.ott_id = taxon["ott_id"].as_u64();
    check.rank = taxon["rank"].as_str().map(str::to_string);
    check.flags = taxon["flags"].as_array().map(|f| f.iter().filter_map(|x| x.as_str().map(str::to_string)).collect()).unwrap_or_default();
    // OTT reports names of "hidden" taxa as approximate matches even when they
    // match exactly; an identical name is exact.
    let same = check.matched.as_deref().is_some_and(|m| m.eq_ignore_ascii_case(query));
    check.status = match rank_of(best).0 {
        2 => NameStatus::Exact,
        1 => NameStatus::Synonym,
        _ if same && check.accepted.as_deref().is_some_and(|a| a.eq_ignore_ascii_case(query)) => NameStatus::Exact,
        _ => NameStatus::Approximate,
    };
    let top = rank_of(best);
    for m in matches {
        if let Some(n) = m["taxon"]["name"].as_str() {
            if rank_of(m) == top && Some(n) != check.accepted.as_deref() && !check.alternatives.iter().any(|a| a == n) {
                check.alternatives.push(n.to_string());
            }
        }
    }
    check
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_match_names_reply() {
        let resp: Value = serde_json::from_str(
            r#"{"results":[
              {"name":"Homo sapiens","matches":[{"is_approximate_match":false,"is_synonym":false,"matched_name":"Homo sapiens","score":1.0,
                "taxon":{"name":"Homo sapiens","ott_id":770315,"rank":"species","flags":["extinct"]}}]},
              {"name":"Pan troglodites","matches":[{"is_approximate_match":true,"is_synonym":false,"matched_name":"Pan troglodytes","score":0.93,
                "taxon":{"name":"Pan troglodytes","ott_id":417950,"rank":"species","flags":[]}}]},
              {"name":"Hylobates syndactylus","matches":[{"is_approximate_match":false,"is_synonym":true,"matched_name":"Hylobates syndactylus","score":1.0,
                "taxon":{"name":"Symphalangus syndactylus","ott_id":417961,"rank":"species","flags":[]}}]},
              {"name":"Notarealus taxonus","matches":[]},
              {"name":"Australopithecus garhi","matches":[{"is_approximate_match":true,"is_synonym":false,"matched_name":"Australopithecus garhi","score":1.0,
                "taxon":{"name":"Australopithecus garhi","ott_id":3607724,"rank":"species","flags":["hidden","extinct"]}}]}
            ]}"#,
        )
        .unwrap();
        let names: Vec<String> =
            ["Homo sapiens", "Pan troglodites", "Hylobates syndactylus", "Notarealus taxonus", "Australopithecus garhi"].iter().map(|s| s.to_string()).collect();
        let r = parse_response(&resp, &names);
        assert_eq!(r[0].status, NameStatus::Exact);
        assert_eq!(r[0].suggestion(), None);
        assert_eq!(r[0].flags, ["extinct"]);
        assert_eq!(r[1].status, NameStatus::Approximate);
        assert_eq!(r[1].suggestion(), Some("Pan troglodytes"));
        assert_eq!(r[2].status, NameStatus::Synonym);
        assert_eq!(r[2].suggestion(), Some("Symphalangus syndactylus"));
        assert_eq!(r[3].status, NameStatus::Unmatched);
        assert_eq!(r[4].status, NameStatus::Exact, "hidden taxa match exactly");
        assert_eq!(query_name(" Homo_sapiens "), "Homo sapiens");
    }
}
