//! The delta test for the node-density artifact (Webster, Payne & Pagel 2003;
//! assessed by Venditti, Meade & Pagel 2006, Syst. Biol. 55: 637–643).
//!
//! For each species: x = root-to-tip path length, n = internal nodes along the
//! path (the root counts as zero; the tip is not counted). The artifact makes
//! longer paths go with more nodes along a curve n = βx^δ with δ > 1.
//!
//! 1. Association: phylogenetic GLS of x on n (random walk vs directional
//!    model, slope β*), likelihood-ratio test against χ²₁.
//! 2. Curvature: δ by maximum likelihood from n = βx^δ (and, for comparison,
//!    from x = β* n^(1/δ)).
//!
//! Residuals are phylogenetically correlated: the covariance of two species is
//! their shared root path, in branch-length units when path length is the
//! response and in node counts (unit branches) when nodes are the response.

use crate::tree::{NodeId, Tree};
use anyhow::{bail, Result};

#[derive(Clone)]
pub struct DeltaTest {
    pub labels: Vec<String>,
    /// Root-to-tip path length per species.
    pub path_lengths: Vec<f64>,
    /// Internal nodes along each path (root = 0, tip not counted).
    pub nodes: Vec<f64>,
    /// Slope of path length on nodes (directional model).
    pub beta_star: f64,
    /// Likelihood ratio, directional vs random-walk model, and its χ²₁ p-value.
    pub lr: f64,
    pub p_value: f64,
    /// δ and β from n = βx^δ (the preferred estimate).
    pub delta: f64,
    pub beta: f64,
    /// δ and β* from x = β* n^(1/δ) (for comparison; better when δ ≫ 1).
    pub delta_alt: f64,
    pub beta_star_alt: f64,
}

impl DeltaTest {
    /// Significant positive association together with δ > 1.
    pub fn artifact(&self) -> bool {
        self.beta_star > 0.0 && self.p_value < 0.05 && self.delta > 1.0
    }

    /// δ for display: values at the edge of the search range are not real estimates
    /// (e.g. all path lengths equal makes the curve vertical).
    pub fn delta_text(d: f64) -> String {
        if d >= 49.5 {
            "> 50 (not determinable)".into()
        } else if d <= 0.051 {
            "< 0.05 (not determinable)".into()
        } else {
            crate::tree::format_num(d, 3)
        }
    }

    /// Plain-text summary plus the per-species data (tab-separated).
    pub fn report(&self) -> String {
        use crate::tree::format_num as f;
        let mut s = String::from("Node-density artifact: delta test\n");
        s += &format!("Species\t{}\n", self.labels.len());
        s += &format!("Slope of path length on nodes (beta*)\t{}\n", f(self.beta_star, 4));
        s += &format!("Likelihood ratio (directional vs random walk)\t{}\n", f(self.lr, 3));
        s += &format!("p (chi-square, 1 df)\t{}\n", f(self.p_value, 4));
        s += &format!("delta (from n = beta x^delta)\t{}\n", Self::delta_text(self.delta));
        s += &format!("delta (from x = beta* n^(1/delta))\t{}\n", Self::delta_text(self.delta_alt));
        s += &format!("Result\t{}\n", self.verdict());
        s += "\nspecies\tpath_length\tnodes\n";
        for i in 0..self.labels.len() {
            s += &format!("{}\t{}\t{}\n", self.labels[i], f(self.path_lengths[i], 6), self.nodes[i]);
        }
        s
    }

    pub fn verdict(&self) -> &'static str {
        if self.artifact() {
            "Node-density artifact indicated (significant positive association and δ > 1)"
        } else if self.beta_star > 0.0 && self.p_value < 0.05 {
            "Significant association, but δ ≤ 1: not the pattern of the node-density artifact"
        } else {
            "No significant association between path length and nodes: no evidence of the artifact"
        }
    }
}

/// Lower Cholesky factor of a symmetric positive-definite matrix (row-major).
fn cholesky(a: &[f64], n: usize) -> Option<Vec<f64>> {
    let mut l = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..=i {
            let mut s = a[i * n + j];
            for k in 0..j {
                s -= l[i * n + k] * l[j * n + k];
            }
            if i == j {
                if s <= 0.0 {
                    return None;
                }
                l[i * n + i] = s.sqrt();
            } else {
                l[i * n + j] = s / l[j * n + j];
            }
        }
    }
    Some(l)
}

/// Solve L y = b (forward substitution): "whitens" a vector.
fn whiten(l: &[f64], n: usize, b: &[f64]) -> Vec<f64> {
    let mut y = vec![0.0; n];
    for i in 0..n {
        let mut s = b[i];
        for k in 0..i {
            s -= l[i * n + k] * y[k];
        }
        y[i] = s / l[i * n + i];
    }
    y
}

/// A phylogenetic covariance structure, factored once.
struct Gls {
    l: Vec<f64>,
    n: usize,
    logdet: f64,
}

impl Gls {
    fn new(v: Vec<f64>, n: usize) -> Result<Self> {
        // A tiny ridge keeps zero-length terminal branches from making V singular.
        let mean_diag = (0..n).map(|i| v[i * n + i]).sum::<f64>() / n as f64;
        let mut v = v;
        for i in 0..n {
            v[i * n + i] += 1e-9 * mean_diag.max(1e-12);
        }
        let Some(l) = cholesky(&v, n) else { bail!("the phylogenetic covariance matrix is not positive definite") };
        let logdet = 2.0 * (0..n).map(|i| l[i * n + i].ln()).sum::<f64>();
        Ok(Gls { l, n, logdet })
    }

    /// ML fit of y on the columns of `x` (each a regressor). Returns (coefficients, log-likelihood).
    fn fit(&self, y: &[f64], x: &[Vec<f64>]) -> (Vec<f64>, f64) {
        let n = self.n;
        let ys = whiten(&self.l, n, y);
        let xs: Vec<Vec<f64>> = x.iter().map(|c| whiten(&self.l, n, c)).collect();
        let k = xs.len();
        // Normal equations (k is 1 or 2 here).
        let mut a = vec![0.0; k * k];
        let mut b = vec![0.0; k];
        for i in 0..k {
            b[i] = xs[i].iter().zip(&ys).map(|(p, q)| p * q).sum();
            for j in 0..k {
                a[i * k + j] = xs[i].iter().zip(&xs[j]).map(|(p, q)| p * q).sum();
            }
        }
        let coef = match k {
            1 => vec![if a[0] > 0.0 { b[0] / a[0] } else { 0.0 }],
            _ => {
                let det = a[0] * a[3] - a[1] * a[2];
                if det.abs() < 1e-300 {
                    vec![0.0, 0.0]
                } else {
                    vec![(a[3] * b[0] - a[1] * b[1]) / det, (a[0] * b[1] - a[2] * b[0]) / det]
                }
            }
        };
        let rss: f64 = (0..n)
            .map(|i| {
                let fit: f64 = (0..k).map(|j| coef[j] * xs[j][i]).sum();
                (ys[i] - fit).powi(2)
            })
            .sum();
        let s2 = (rss / n as f64).max(1e-300);
        let ll = -0.5 * n as f64 * ((2.0 * std::f64::consts::PI * s2).ln() + 1.0) - 0.5 * self.logdet;
        (coef, ll)
    }
}

/// Maximize `ll(δ)` over δ in [0.05, 50]: log-spaced grid, then golden-section refinement.
fn maximize(ll: impl Fn(f64) -> f64) -> f64 {
    let (lo, hi) = (0.05f64.ln(), 50f64.ln());
    let steps = 120;
    let grid: Vec<f64> = (0..=steps).map(|i| lo + (hi - lo) * i as f64 / steps as f64).collect();
    let vals: Vec<f64> = grid.iter().map(|&g| ll(g.exp())).collect();
    let best = (0..vals.len()).max_by(|&a, &b| vals[a].partial_cmp(&vals[b]).unwrap_or(std::cmp::Ordering::Equal)).unwrap_or(0);
    let (mut a, mut b) = (grid[best.saturating_sub(1)], grid[(best + 1).min(steps)]);
    let r = (5f64.sqrt() - 1.0) / 2.0;
    for _ in 0..60 {
        let c = b - r * (b - a);
        let d = a + r * (b - a);
        if ll(c.exp()) > ll(d.exp()) {
            b = d;
        } else {
            a = c;
        }
    }
    ((a + b) / 2.0).exp()
}

/// Run the delta test on the whole tree (it should be rooted and have branch lengths).
pub fn delta_test(tree: &Tree) -> Result<DeltaTest> {
    if !tree.has_lengths() {
        bail!("the delta test needs branch lengths");
    }
    let tips: Vec<NodeId> = tree.tips();
    let m = tips.len();
    if m < 4 {
        bail!("the delta test needs at least 4 species");
    }
    let depth = tree.depths(true);
    let mut edges = vec![0usize; tree.nodes.len()];
    for u in tree.preorder() {
        if let Some(p) = tree.nodes[u].parent {
            edges[u] = edges[p] + 1;
        }
    }
    let x: Vec<f64> = tips.iter().map(|&t| depth[t]).collect();
    // Internal nodes strictly between the root and the tip.
    let nn: Vec<f64> = tips.iter().map(|&t| edges[t].saturating_sub(1) as f64).collect();
    if nn.iter().all(|&v| v == nn[0]) {
        bail!("every species has the same number of nodes to the root (a perfectly balanced tree cannot show the artifact)");
    }
    // Shared root paths: branch lengths (for x) and edge counts (for n).
    let mut vx = vec![0.0; m * m];
    let mut vn = vec![0.0; m * m];
    for i in 0..m {
        for j in i..m {
            let a = if i == j { tips[i] } else { tree.mrca(&[tips[i], tips[j]]).unwrap_or(tree.root) };
            vx[i * m + j] = depth[a];
            vx[j * m + i] = depth[a];
            vn[i * m + j] = edges[a] as f64;
            vn[j * m + i] = edges[a] as f64;
        }
    }
    let gx = Gls::new(vx, m)?;
    let gn = Gls::new(vn, m)?;
    let ones = vec![1.0; m];

    // 1. Random walk (x = α) vs directional (x = α + β* n).
    let (_, ll0) = gx.fit(&x, &[ones.clone()]);
    let (c1, ll1) = gx.fit(&x, &[ones.clone(), nn.clone()]);
    let lr = (2.0 * (ll1 - ll0)).max(0.0);
    let p_value = 2.0 * (1.0 - crate::tree::normal_cdf(lr.sqrt()));

    // 2. δ from n = β x^δ (nodes as response).
    let ll_n = |d: f64| gn.fit(&nn, &[x.iter().map(|v| v.max(0.0).powf(d)).collect()]).1;
    let delta = maximize(ll_n);
    let beta = gn.fit(&nn, &[x.iter().map(|v| v.max(0.0).powf(delta)).collect()]).0[0];

    // δ from x = β* n^(1/δ) (path length as response).
    let ll_x = |d: f64| gx.fit(&x, &[nn.iter().map(|v| v.powf(1.0 / d)).collect()]).1;
    let delta_alt = maximize(ll_x);
    let beta_star_alt = gx.fit(&x, &[nn.iter().map(|v| v.powf(1.0 / delta_alt)).collect()]).0[0];

    Ok(DeltaTest {
        labels: tips.iter().map(|&t| tree.label(t).to_string()).collect(),
        path_lengths: x,
        nodes: nn,
        beta_star: c1[1],
        lr,
        p_value,
        delta,
        beta,
        delta_alt,
        beta_star_alt,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::newick::parse_newick;

    #[test]
    fn ultrametric_tree_shows_no_artifact() {
        // Ultrametric caterpillar: path lengths are identical, so no association.
        let t = parse_newick("(((((A:1,B:1):1,C:2):1,D:3):1,E:4):1,F:5);").unwrap();
        let r = delta_test(&t).unwrap();
        assert_eq!(r.nodes, vec![4.0, 4.0, 3.0, 2.0, 1.0, 0.0]);
        assert!(r.lr < 1e-6, "LR {}", r.lr);
        assert!(!r.artifact());
    }

    #[test]
    fn curved_path_lengths_show_the_artifact() {
        // A ladder: the species on spine node k has k nodes to the root, and path
        // lengths that rise with k but level off (x ∝ (k+1)^(1/3)): the curvilinear
        // signature of the artifact, so δ > 1.
        let kmax = 14;
        fn spine(k: usize, kmax: usize) -> String {
            let x = |k: usize| 2.0 * ((k + 1) as f64).powf(1.0 / 3.0) + if k % 2 == 0 { 0.02 } else { -0.02 };
            let d = 0.1 * k as f64;
            if k == kmax {
                format!("(T{}a:{:.4},T{}b:{:.4})", k, x(k) - d, k, x(k) - d + 0.01)
            } else {
                format!("(T{}:{:.4},{}:0.1)", k, x(k) - d, spine(k + 1, kmax))
            }
        }
        let s = spine(0, kmax);        let t = parse_newick(&format!("{};", s)).unwrap();
        let r = delta_test(&t).unwrap();
        assert!(r.beta_star > 0.0 && r.p_value < 0.05, "β* {} p {}", r.beta_star, r.p_value);
        assert!(r.delta > 1.0, "δ {}", r.delta);
    }

    #[test]
    fn needs_lengths_and_unequal_nodes() {
        assert!(delta_test(&parse_newick("((A,B),(C,D));").unwrap()).is_err());
        assert!(delta_test(&parse_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap()).is_err());
    }
}
