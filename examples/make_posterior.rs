//! Generate a small BEAST-style posterior sample for trying the consensus /
//! MCC tools: `cargo run --example make_posterior`
//! Writes examples/primates_posterior.trees (NEXUS with a TRANSLATE table).

use canopy::io::newick::{parse_newick, WriteOptions};
use canopy::io::nexus::write_nexus;

/// Small deterministic PRNG (xorshift) so the example needs no dependencies.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
    /// Approximately normal jitter around `h` (sd = 8%).
    fn jitter(&mut self, h: f64) -> f64 {
        let z: f64 = (0..6).map(|_| self.next()).sum::<f64>() - 3.0;
        (h * (1.0 + 0.08 * z * 1.41)).max(0.1)
    }
}

fn sample(r: &mut Rng) -> String {
    let swap_hpg = r.next() < 0.18;
    let swap_colobine = r.next() < 0.3;
    let hp = r.jitter(6.4);
    let hpg = r.jitter(8.8).max(hp + 0.2);
    let hominini = if swap_hpg {
        format!("(Homo_sapiens:{hp:.3},Gorilla_gorilla:{hp:.3}):{:.3},Pan_troglodytes:{hpg:.3}", hpg - hp)
    } else {
        format!("(Homo_sapiens:{hp:.3},Pan_troglodytes:{hp:.3}):{:.3},Gorilla_gorilla:{hpg:.3}", hpg - hp)
    };
    let hpgp = r.jitter(15.7).max(hpg + 0.2);
    let hom = r.jitter(20.2).max(hpgp + 0.2);
    let apes = format!("((({hominini}):{:.3},Pongo_abelii:{hpgp:.3}):{:.3},Hylobates_lar:{hom:.3})", hpgp - hpg, hom - hpgp);
    let mp = r.jitter(10.6);
    let co = r.jitter(11.1);
    let owm = mp.max(co) + r.jitter(7.5);
    let (a, b, c, d) = if swap_colobine {
        ("Macaca_mulatta", "Colobus_guereza", "Papio_anubis", "Nasalis_larvatus")
    } else {
        ("Macaca_mulatta", "Papio_anubis", "Colobus_guereza", "Nasalis_larvatus")
    };
    let cerco = format!("(({a}:{mp:.3},{b}:{mp:.3}):{:.3},({c}:{co:.3},{d}:{co:.3}):{:.3})", owm - mp, owm - co);
    let root = hom.max(owm) + r.jitter(9.6);
    format!("({apes}:{:.3},{cerco}:{:.3});", root - hom, root - owm)
}

fn main() {
    let mut r = Rng(0x9E3779B97F4A7C15);
    let mut trees = Vec::new();
    for i in 0..200 {
        let mut t = parse_newick(&sample(&mut r)).expect("valid newick");
        t.name = Some(format!("STATE_{}", i * 1000));
        t.rooted = Some(true);
        trees.push(t);
    }
    let refs: Vec<_> = trees.iter().collect();
    let text = write_nexus(&refs, &WriteOptions::default());
    std::fs::write("examples/primates_posterior.trees", text).expect("write file");
    println!("wrote examples/primates_posterior.trees ({} trees)", trees.len());
}
