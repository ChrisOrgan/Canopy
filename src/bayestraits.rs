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

#[cfg(test)]
mod tests {
    use super::*;

    const HEAD: &str = "2\n0\tA\n1\tB\n3\n0\t-1\t2\t0\t1\n1\t1\t1\t0\n2\t1\t1\t1\nIt\tLh\tNode ID\tScaler\tCreat It\tNode / Branch\t\n";

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
