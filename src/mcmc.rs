//! Convergence diagnostics for several MCMC runs of the same analysis
//! (Vehtari, Gelman, Simpson, Carpenter & Bürkner 2021): rank-normalized split
//! R-hat (the larger of the bulk and tail versions) and bulk effective sample size,
//! as in Stan and ArviZ. Chains are cut to the shortest run's length.

/// Inverse of the standard normal CDF (Acklam's rational approximation,
/// relative error below 1.2e-9).
pub fn normal_quantile(p: f64) -> f64 {
    const A: [f64; 6] = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00];
    const B: [f64; 5] = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01];
    const C: [f64; 6] = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00];
    const D: [f64; 4] = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00];
    let p = p.clamp(1e-300, 1.0 - 1e-16);
    let (lo, hi) = (0.02425, 1.0 - 0.02425);
    if p < lo {
        let q = (-2.0 * p.ln()).sqrt();
        (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5]) / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    } else if p <= hi {
        let q = p - 0.5;
        let r = q * q;
        (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0)
    } else {
        let q = (-2.0 * (1.0 - p).ln()).sqrt();
        -(((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5]) / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    }
}

/// The chains cut to the shortest length, then each split in half (split
/// chains detect trends within a run as well as differences between runs).
fn split(chains: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = chains.iter().map(Vec::len).min().unwrap_or(0);
    let half = n / 2;
    if half < 2 {
        return Vec::new();
    }
    chains.iter().flat_map(|c| [c[..half].to_vec(), c[n - half..n].to_vec()]).collect()
}

/// Replace values by normal scores of their ranks over all chains (ties get
/// their average rank): z = Φ⁻¹((r − 3/8) / (S + 1/4)).
fn rank_normalize(chains: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let mut all: Vec<(f64, usize, usize)> = chains.iter().enumerate().flat_map(|(c, v)| v.iter().enumerate().map(move |(i, &x)| (x, c, i))).collect();
    all.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let s = all.len() as f64;
    let mut out: Vec<Vec<f64>> = chains.iter().map(|c| vec![0.0; c.len()]).collect();
    let mut i = 0;
    while i < all.len() {
        let mut j = i;
        while j + 1 < all.len() && all[j + 1].0 == all[i].0 {
            j += 1;
        }
        let rank = (i + j) as f64 / 2.0 + 1.0;
        let z = normal_quantile((rank - 0.375) / (s + 0.25));
        for e in &all[i..=j] {
            out[e.1][e.2] = z;
        }
        i = j + 1;
    }
    out
}

fn mean(v: &[f64]) -> f64 {
    v.iter().sum::<f64>() / v.len().max(1) as f64
}

fn var(v: &[f64]) -> f64 {
    let m = mean(v);
    v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (v.len() as f64 - 1.0).max(1.0)
}

/// Classic R-hat of equal-length chains; NaN when the chains do not vary.
fn basic_rhat(chains: &[Vec<f64>]) -> f64 {
    let n = chains[0].len() as f64;
    let w = mean(&chains.iter().map(|c| var(c)).collect::<Vec<_>>());
    let b = n * var(&chains.iter().map(|c| mean(c)).collect::<Vec<_>>());
    if w <= 0.0 {
        return f64::NAN;
    }
    (((n - 1.0) / n * w + b / n) / w).sqrt()
}

/// Rank-normalized split R-hat: the larger of the bulk (ranks of the values) and
/// tail (ranks of distances from the median) versions. Values near 1 mean the
/// runs agree; above 1.1, inadequate convergence. NaN without variation
/// (e.g. a branch never rescaled) or with too few iterations.
pub fn rhat(chains: &[Vec<f64>]) -> f64 {
    let s = split(chains);
    if s.len() < 2 || !s.iter().flatten().any(|x| *x != s[0][0]) {
        return f64::NAN;
    }
    let bulk = basic_rhat(&rank_normalize(&s));
    let mut pooled: Vec<f64> = s.iter().flatten().copied().collect();
    pooled.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let med = crate::io::data::quantile(&pooled, 0.5);
    let folded: Vec<Vec<f64>> = s.iter().map(|c| c.iter().map(|x| (x - med).abs()).collect()).collect();
    let tail = basic_rhat(&rank_normalize(&folded));
    match (bulk.is_nan(), tail.is_nan()) {
        (true, true) => f64::NAN,
        (true, false) => tail,
        (false, true) => bulk,
        _ => bulk.max(tail),
    }
}

/// Effective sample size of equal-length chains, as `ess_basic` in the R package
/// posterior (Stan's reference): the multi-chain autocorrelation estimate
/// ρ̂_t = 1 − (W − mean autocovariance_t) / var⁺, Geyer's initial positive
/// sequence (pairs summed while positive), made monotone by averaging, and
/// τ = −1 + 2 Σ ρ̂_t (+ the next even lag), at least 1/log10(draws).
fn basic_ess(chains: &[Vec<f64>]) -> f64 {
    let m = chains.len();
    let n = chains[0].len();
    if n < 4 {
        return f64::NAN;
    }
    let means: Vec<f64> = chains.iter().map(|c| mean(c)).collect();
    // Autocovariances of each chain (biased, divided by n), averaged over chains.
    let acov = |t: usize| -> f64 { chains.iter().zip(&means).map(|(c, &mu)| (0..n - t).map(|i| (c[i] - mu) * (c[i + t] - mu)).sum::<f64>() / n as f64).sum::<f64>() / m as f64 };
    let mean_var = acov(0) * n as f64 / (n as f64 - 1.0);
    let mut var_plus = mean_var * (n as f64 - 1.0) / n as f64;
    if m > 1 {
        var_plus += var(&means);
    }
    if !(var_plus > 0.0) {
        return f64::NAN;
    }
    let rho_at = |t: usize| 1.0 - (mean_var - acov(t)) / var_plus;
    let mut rho = vec![0.0; n];
    rho[0] = 1.0;
    let mut even = 1.0;
    let mut odd = rho_at(1);
    rho[1] = odd;
    let mut t = 0;
    while t + 5 < n && even + odd > 0.0 {
        t += 2;
        even = rho_at(t);
        odd = rho_at(t + 1);
        if even + odd >= 0.0 {
            rho[t] = even;
            rho[t + 1] = odd;
        }
    }
    let max_t = t;
    if even > 0.0 {
        rho[max_t] = even;
    }
    // Initial monotone sequence.
    let mut t = 0;
    while t + 4 <= max_t {
        t += 2;
        if rho[t] + rho[t + 1] > rho[t - 2] + rho[t - 1] {
            rho[t] = (rho[t - 2] + rho[t - 1]) / 2.0;
            rho[t + 1] = rho[t];
        }
    }
    let total = (m * n) as f64;
    let tau = -1.0 + 2.0 * rho[..max_t].iter().sum::<f64>() + rho[max_t];
    total / tau.max(1.0 / total.log10())
}

/// Bulk effective sample size (rank-normalized split chains): about how many
/// independent draws the runs together are worth for estimating the
/// posterior's centre. Above 400 is the usual target. NaN without variation.
pub fn ess(chains: &[Vec<f64>]) -> f64 {
    let s = split(chains);
    if s.is_empty() || !s.iter().flatten().any(|x| *x != s[0][0]) {
        return f64::NAN;
    }
    basic_ess(&rank_normalize(&s))
}

/// Sample standard deviation across runs (divided by k − 1), as MrBayes computes
/// the standard deviation of split frequencies for its ASDSF.
pub fn sd(v: &[f64]) -> f64 {
    if v.len() < 2 {
        return 0.0;
    }
    var(v).sqrt()
}

/// Across-run comparison of one branch: posterior probabilities of a scalar
/// (or β) in each run, and R-hat / ESS of its value trace.
#[derive(Clone, Debug)]
pub struct BranchCheck {
    pub id: usize,
    pub n_taxa: usize,
    /// Per run: posterior probability of the branch's scalar (variable rates) or β (Fabric).
    pub pp: Vec<f64>,
    /// Per run: posterior probability of a node scalar at this node.
    pub node_pp: Vec<f64>,
    pub rhat: f64,
    pub ess: f64,
}

/// Convergence report for runs combined into one posterior.
#[derive(Clone, Debug)]
pub struct Report {
    /// Per run: name, iterations, mean log-likelihood.
    pub runs: Vec<(String, usize, f64)>,
    pub lh_rhat: f64,
    pub lh_ess: f64,
    pub branches: Vec<BranchCheck>,
    /// Name of the per-branch quantity ("rate scalar" or "LS_Beta").
    pub quantity: &'static str,
}

impl Report {
    /// Average standard deviation of the posterior probabilities across runs
    /// (branch scalars or β, and node scalars), over those reaching 0.1 in some
    /// run, as MrBayes' ASDSF does for clades; and the largest difference.
    pub fn asd_pp(&self) -> (f64, f64) {
        let mut sds = Vec::new();
        let mut max_diff: f64 = 0.0;
        for b in &self.branches {
            for v in [&b.pp, &b.node_pp] {
                let hi = v.iter().cloned().fold(0.0, f64::max);
                let lo = v.iter().cloned().fold(1.0, f64::min);
                max_diff = max_diff.max(hi - lo);
                if hi >= 0.1 {
                    sds.push(sd(v));
                }
            }
        }
        (if sds.is_empty() { 0.0 } else { mean(&sds) }, max_diff)
    }

    /// Plain-language findings, each with whether it passes the usual targets
    /// (log-likelihood R-hat at most 1.1; average SD of posterior
    /// probabilities < 0.01, MrBayes' guide for ASDSF, up to 0.05 acceptable).
    pub fn findings(&self) -> Vec<(bool, String)> {
        let f = |x: f64, d: usize| if x.is_finite() { crate::tree::format_num(x, d) } else { "—".into() };
        let total: usize = self.runs.iter().map(|r| r.1).sum();
        let mut out = vec![(true, format!("{} iterations pooled from {} runs", total, self.runs.len()))];
        let lens: Vec<usize> = self.runs.iter().map(|r| r.1).collect();
        if lens.iter().any(|&l| l != lens[0]) {
            out.push((true, format!("Runs differ in length; R-hat and ESS use the first {} iterations of each", lens.iter().min().unwrap_or(&0))));
        }
        out.push((self.lh_rhat.is_nan() || self.lh_rhat <= 1.1, format!("Log-likelihood: R-hat = {}, ESS = {}", f(self.lh_rhat, 3), f(self.lh_ess, 0))));
        let (asd, diff) = self.asd_pp();
        out.push((asd < 0.05, format!("Posterior probabilities across runs: average SD {} (largest difference {})", f(asd, 3), f(diff, 3))));
        out
    }

    /// Per-branch comparison table: posterior probabilities in each run with
    /// their SD across runs (R-hat and ESS of each branch's value are in the
    /// diagnostics' parameter table).
    pub fn table(&self) -> (Vec<String>, Vec<Vec<String>>) {
        let k = self.runs.len();
        let what = if self.quantity == "LS_Beta" { "β" } else { "Branch scalar" };
        let mut header = vec!["Node".to_string(), "Taxa".to_string()];
        header.extend((1..=k).map(|i| format!("{} PP run {}", what, i)));
        header.push(format!("{} PP SD", what));
        header.extend((1..=k).map(|i| format!("Node scalar PP run {}", i)));
        header.push("Node scalar PP SD".to_string());
        let f = |x: f64, d: usize| if x.is_finite() { crate::tree::format_num(x, d) } else { "—".into() };
        let rows = self
            .branches
            .iter()
            .map(|b| {
                let mut r = vec![b.id.to_string(), b.n_taxa.to_string()];
                r.extend(b.pp.iter().map(|x| f(*x, 3)));
                r.push(f(sd(&b.pp), 3));
                r.extend(b.node_pp.iter().map(|x| f(*x, 3)));
                r.push(f(sd(&b.node_pp), 3));
                r
            })
            .collect();
        (header, rows)
    }

    /// The largest R-hat and the smallest ESS over the log-likelihood and branches.
    pub fn worst(&self) -> (f64, f64) {
        let rh = std::iter::once(self.lh_rhat).chain(self.branches.iter().map(|b| b.rhat)).filter(|x| x.is_finite()).fold(f64::NAN, f64::max);
        let es = std::iter::once(self.lh_ess).chain(self.branches.iter().map(|b| b.ess)).filter(|x| x.is_finite()).fold(f64::NAN, f64::min);
        (rh, es)
    }
}


/// One parameter's samples in each run.
#[derive(Clone, Debug)]
pub struct Param {
    pub name: String,
    pub chains: Vec<Vec<f64>>,
}

impl Param {
    /// All runs' samples together.
    pub fn pooled(&self) -> Vec<f64> {
        self.chains.iter().flatten().copied().collect()
    }
}

/// Everything the MCMC diagnostics window shows: the runs, each parameter's
/// samples per run, and for variable-rates or Fabric runs the comparison of
/// posterior probabilities across runs.
#[derive(Clone, Debug)]
pub struct Diagnostics {
    pub title: String,
    /// Per run: name and number of samples.
    pub runs: Vec<(String, usize)>,
    pub params: Vec<Param>,
    /// Posterior probabilities of scalars per run (two or more VarRates runs).
    pub report: Option<Report>,
}

impl Diagnostics {
    /// From BayesTraits logs (one per run: its name and its columns). Columns
    /// are matched by name; those missing from a run are left out.
    pub fn from_logs(title: &str, logs: &[(String, Vec<(String, Vec<f64>)>)]) -> Diagnostics {
        let names: Vec<String> = logs
            .first()
            .map(|l| l.1.iter().map(|c| c.0.clone()).filter(|n| logs.iter().all(|r| r.1.iter().any(|c| c.0 == *n))).collect())
            .unwrap_or_default();
        let params = names
            .into_iter()
            .map(|n| Param { chains: logs.iter().map(|r| r.1.iter().find(|c| c.0 == n).map(|c| c.1.clone()).unwrap_or_default()).collect(), name: n })
            .collect();
        Diagnostics { title: title.to_string(), runs: logs.iter().map(|l| (l.0.clone(), l.1.first().map(|c| c.1.len()).unwrap_or(0))).collect(), params, report: None }
    }

    /// Leave out run `i` (its samples of every parameter); the title, if it
    /// named the runs joined by " + ", names the runs that are left.
    pub fn remove_run(&mut self, i: usize) {
        if i >= self.runs.len() {
            return;
        }
        let joined = |runs: &[(String, usize)]| runs.iter().map(|r| r.0.as_str()).collect::<Vec<_>>().join(" + ");
        let renamed = self.title == joined(&self.runs);
        self.runs.remove(i);
        for p in self.params.iter_mut() {
            if i < p.chains.len() {
                p.chains.remove(i);
            }
        }
        self.report = None;
        if renamed {
            self.title = joined(&self.runs);
        }
    }

    /// Add parameters not already present (e.g. a run's log next to its VarRates).
    pub fn add_params(&mut self, more: Vec<Param>) {
        for p in more {
            if !self.params.iter().any(|q| q.name == p.name) {
                self.params.push(p);
            }
        }
    }

    /// One row per parameter: mean, median and 95% HPD of all runs pooled, then
    /// R-hat and ESS across runs (split within a single run).
    pub fn table(&self) -> (Vec<String>, Vec<Vec<String>>) {
        let header: Vec<String> = ["Parameter", "Mean", "Median", "95% HPD lower", "95% HPD upper", "R-hat", "ESS"].iter().map(|s| s.to_string()).collect();
        let f = |x: f64, d: usize| if x.is_finite() { crate::tree::format_num(x, d) } else { "—".into() };
        let rows = self
            .params
            .iter()
            .map(|p| {
                let mut v = p.pooled();
                v.retain(|x| x.is_finite());
                v.sort_by(|a, b| a.partial_cmp(b).unwrap());
                if v.is_empty() {
                    return vec![p.name.clone(), "—".into(), "—".into(), "—".into(), "—".into(), "—".into(), "—".into()];
                }
                let (lo, hi) = crate::consensus::hpd(&v, 0.95).unwrap_or((f64::NAN, f64::NAN));
                vec![p.name.clone(), f(mean(&v), 3), f(crate::io::data::quantile(&v, 0.5), 3), f(lo, 3), f(hi, 3), f(rhat(&p.chains), 3), f(ess(&p.chains), 0)]
            })
            .collect();
        (header, rows)
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic standard normal draws (xorshift + Box–Muller).
    fn normals(seed: u64, n: usize) -> Vec<f64> {
        let mut s = seed.max(1);
        let mut u = || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            (s >> 11) as f64 / (1u64 << 53) as f64
        };
        (0..n).map(|_| (-2.0 * u().max(1e-300).ln()).sqrt() * (2.0 * std::f64::consts::PI * u()).cos()).collect()
    }

    #[test]
    fn quantile_matches_known_values() {
        assert!(normal_quantile(0.5).abs() < 1e-9);
        assert!((normal_quantile(0.975) - 1.959963985).abs() < 1e-6);
        assert!((normal_quantile(0.001) + 3.090232306).abs() < 1e-6);
    }

    #[test]
    fn rhat_and_ess_behave() {
        // Four well-mixed runs of independent draws: R-hat ≈ 1, ESS ≈ all draws.
        let good: Vec<Vec<f64>> = (1..=4).map(|k| normals(k * 7919, 1000)).collect();
        let r = rhat(&good);
        assert!(r < 1.01, "{}", r);
        let e = ess(&good);
        assert!(e > 3500.0 && e < 4600.0, "iid ESS {} (theory 4000)", e);
        // One run stuck elsewhere: R-hat well above 1.
        let mut bad = good.clone();
        bad[3].iter_mut().for_each(|x| *x += 3.0);
        assert!(rhat(&bad) > 1.1);
        // Strong autocorrelation (AR(1), φ = 0.9): ESS far below the draws (≈ n·(1−φ)/(1+φ)).
        let ar: Vec<Vec<f64>> = (1..=4)
            .map(|k| {
                let z = normals(k * 104729, 1000);
                let mut x = 0.0;
                z.iter().map(|e| {
                    x = 0.9 * x + e;
                    x
                }).collect()
            })
            .collect();
        let e = ess(&ar);
        assert!(e > 160.0 && e < 270.0, "AR(1) ESS {} (theory 211)", e);
        // R-hat of the AR(1) runs is still close to 1 (they mix, just slowly).
        assert!(rhat(&ar) < 1.05, "{}", rhat(&ar));
        // No variation: undefined.
        assert!(rhat(&[vec![1.0; 50], vec![1.0; 50]]).is_nan());
    }

    #[test]
    fn removes_a_run() {
        let col = |v: f64, n: usize| vec![("Lh".to_string(), vec![v; n])];
        let logs = vec![("a".to_string(), col(1.0, 3)), ("b".to_string(), col(2.0, 4)), ("c".to_string(), col(3.0, 5))];
        let mut d = Diagnostics::from_logs("a + b + c", &logs);
        d.remove_run(1);
        assert_eq!(d.title, "a + c");
        assert_eq!(d.runs, vec![("a".to_string(), 3), ("c".to_string(), 5)]);
        assert_eq!(d.params[0].chains, vec![vec![1.0; 3], vec![3.0; 5]]);
    }
}
