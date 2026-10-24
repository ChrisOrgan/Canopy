//! Tells BayesTraits VarRates outputs apart: the "Node / Branch" field of the RJ
//! MCMC entries reads "Branch" in variable-rates runs and "LS_Beta" in
//! Fabric-model runs (both also log "Node" scalars).

use anyhow::{bail, Result};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RunKind {
    VariableRates,
    Fabric,
}

impl RunKind {
    pub fn name(self) -> &'static str {
        match self {
            RunKind::VariableRates => "variable rates",
            RunKind::Fabric => "Fabric",
        }
    }
}

/// The kind of run a VarRates file comes from: `Some` if it has Branch or
/// LS_Beta entries, `None` if it only has Node scalars (either kind).
pub fn detect(text: &str) -> Result<Option<RunKind>> {
    let mut lines = text.lines();
    // Entries follow the MCMC header (the line with the "Node ID" column).
    if !lines.by_ref().any(|l| l.split('\t').any(|c| c.trim().eq_ignore_ascii_case("Node ID"))) {
        bail!("no RJ MCMC section (\"Node ID\" column) found: is this a BayesTraits VarRates file?");
    }
    let (mut branch, mut beta) = (false, false);
    for l in lines {
        for c in l.split('\t').map(str::trim) {
            branch |= c.eq_ignore_ascii_case("Branch");
            beta |= c.eq_ignore_ascii_case("LS_Beta");
        }
    }
    match (branch, beta) {
        (true, true) => bail!("the file has both Branch and LS_Beta entries, so it is neither a plain variable-rates nor a Fabric run"),
        (true, false) => Ok(Some(RunKind::VariableRates)),
        (false, true) => Ok(Some(RunKind::Fabric)),
        (false, false) => Ok(None),
    }
}


/// Whether a text file is a BayesTraits VarRates file: it has the RJ MCMC
/// header, starting "It" and with a "Node ID" column.
pub fn is_varrates(text: &str) -> bool {
    text.lines().any(|l| {
        let mut f = l.split('\t').map(str::trim);
        f.next() == Some("It") && f.any(|c| c.eq_ignore_ascii_case("Node ID"))
    })
}

/// Whether a text file is a BayesTraits log: it has a table headed
/// "Iteration", "Lh".
pub fn is_log(text: &str) -> bool {
    text.lines().any(|l| {
        let mut f = l.split('\t').map(str::trim);
        f.next().map(|c| c.eq_ignore_ascii_case("Iteration")).unwrap_or(false) && f.next() == Some("Lh")
    })
}

/// A BayesTraits MCMC log (`*.Log.txt`): the settings block, then a table whose
/// header starts with "Iteration", one row per sample (BayesTraits has already
/// left out the burn-in). Returns each numeric column by name, without the
/// iteration counter; columns with any non-numeric value are left out.
pub fn parse_log(text: &str) -> Result<Vec<(String, Vec<f64>)>> {
    let mut lines = text.lines();
    let Some(header) = lines.by_ref().find(|l| l.split('\t').next().map(|c| c.trim().eq_ignore_ascii_case("Iteration")).unwrap_or(false)) else {
        bail!("no table with an \"Iteration\" column found: is this a BayesTraits log file?")
    };
    let names: Vec<String> = header.split('\t').map(|c| c.trim().to_string()).collect();
    let mut cols: Vec<Vec<f64>> = vec![Vec::new(); names.len()];
    let mut numeric = vec![true; names.len()];
    for l in lines {
        let f: Vec<&str> = l.split('\t').map(str::trim).collect();
        if f.first().map(|s| s.parse::<f64>().is_err()).unwrap_or(true) {
            continue;
        }
        for (i, col) in cols.iter_mut().enumerate() {
            match f.get(i).and_then(|s| s.parse::<f64>().ok()) {
                Some(v) => col.push(v),
                None => numeric[i] = false,
            }
        }
    }
    let out: Vec<(String, Vec<f64>)> = names.into_iter().zip(cols).zip(numeric).skip(1).filter(|((n, c), ok)| *ok && !n.is_empty() && !c.is_empty()).map(|(nc, _)| nc).collect();
    if out.is_empty() {
        bail!("the log has no samples");
    }
    Ok(out)
}
#[cfg(test)]
mod tests {
    use super::*;

    const HEAD: &str = "2\n0\tA\n1\tB\n3\n0\t-1\t2\t0\t1\n1\t1\t1\t0\n2\t1\t1\t1\nIt\tLh\tNode ID\tScaler\tCreat It\tNode / Branch\t\n";


    #[test]
    fn reads_a_log() {
        let text = "Options:\nIterations:\t1010000\nIteration\tLh\tTree No\tAlpha 1\tModel string\t\n10000\t-44.1\t1\t14.15\t'0 1'\t\n11000\t-45.6\t1\t19.58\t'0 0'\t\n";
        let cols = parse_log(text).unwrap();
        let names: Vec<&str> = cols.iter().map(|c| c.0.as_str()).collect();
        // The iteration counter and the non-numeric model string are left out.
        assert_eq!(names, ["Lh", "Tree No", "Alpha 1"]);
        assert_eq!(cols[0].1, vec![-44.1, -45.6]);
        assert!(parse_log("no table here").is_err());
        // The bundled example log.
        if let Ok(t) = std::fs::read_to_string("examples/primates.MCC.VR.Log.txt") {
            let c = parse_log(&t).unwrap();
            assert_eq!(c[0].0, "Lh");
            assert_eq!(c[0].1.len(), 1001);
        }
    }
    #[test]
    fn tells_runs_apart() {
        let vr = format!("{}1\t-1\t1\t2.0\t5\tNode\t2\t1.5\t6\tBranch\t\n", HEAD);
        let fab = format!("{}1\t-1\t1\t2.0\t5\tNode\t2\t-0.4\t6\tLS_Beta\t\n", HEAD);
        let nodes_only = format!("{}1\t-1\t1\t2.0\t5\tNode\t\n", HEAD);
        let both = format!("{}1\t-1\t1\t1.5\t6\tBranch\t2\t-0.4\t6\tLS_Beta\t\n", HEAD);
        assert_eq!(detect(&vr).unwrap(), Some(RunKind::VariableRates));
        assert_eq!(detect(&fab).unwrap(), Some(RunKind::Fabric));
        assert_eq!(detect(&nodes_only).unwrap(), None);
        assert!(detect(&both).is_err());
        assert!(detect("not a varrates file").is_err());
    }
}
