//! Statistical plots drawn as renderer-independent scenes, so the same drawing
//! appears on screen and in saved PNG / SVG / PDF files.

use crate::io::data::{histogram, kernel_density, quantile};
use crate::scene::{Prim, Scene, P};
use crate::style::Color;
use crate::tree::format_num;

/// A round step for about `target` axis ticks over `range`.
pub fn tick_step(range: f64, target: f64) -> f64 {
    let raw = (range / target).max(1e-12);
    let mag = 10f64.powf(raw.log10().floor());
    let f = raw / mag;
    mag * if f < 1.5 {
        1.0
    } else if f < 3.5 {
        2.0
    } else if f < 7.5 {
        5.0
    } else {
        10.0
    }
}

fn text(s: &mut Scene, pos: P, t: String, size: f32, color: Color, halign: f32, valign: f32) {
    s.prims.push(Prim::Text { pos, text: t, size, color, angle: 0.0, halign, valign, italic: false, bold: false });
}

fn line(s: &mut Scene, a: P, b: P, color: Color, width: f32, dashed: bool) {
    s.prims.push(Prim::Path { pts: vec![a, b], color, width, dashed });
}

fn rect(s: &mut Scene, x0: f32, y0: f32, x1: f32, y1: f32, fill: Color) {
    s.prims.push(Prim::Poly { pts: vec![[x0, y0], [x1, y0], [x1, y1], [x0, y1]], fill, stroke: None });
}

/// Distribution of one column, `w` × `h` pixels with sizes in points times `pt`:
/// for numbers, a histogram scaled to a density with a kernel density curve and
/// the central 95% interval dashed; for categories (`numbers` None), a bar chart
/// of the states (at most 30). `title` goes at the top left.
pub fn distribution(numbers: Option<&[f64]>, states: &[(String, usize)], title: &str, w: f32, h: f32, pt: f32) -> Scene {
    let mut s = Scene { width: w, height: h, ..Default::default() };
    let (left, right, top, bottom) = (52.0 * pt, w - 14.0 * pt, 26.0 * pt, h - 34.0 * pt);
    let ink = Color::rgb(50, 50, 50);
    let grid = Color::rgb(232, 232, 232);
    let fill = Color::rgb(86, 140, 200).with_alpha(130);
    let red = Color::rgb(198, 40, 40);
    let size = 9.0 * pt;
    text(&mut s, [left, 6.0 * pt], title.to_string(), 11.0 * pt, ink, 0.0, 0.0);
    let empty = |s: &mut Scene, msg: &str| text(s, [(left + right) / 2.0, (top + bottom) / 2.0], msg.to_string(), 11.0 * pt, Color::GREY, 0.5, 0.5);
    match numbers {
        Some(xs) if xs.len() >= 2 => {
            let (edges, heights) = histogram(xs);
            let curve = kernel_density(xs, 200);
            let mut sorted = xs.to_vec();
            sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let (lo, hi) = (curve[0].0.min(edges[0]), curve[curve.len() - 1].0.max(edges[edges.len() - 1]));
            let ymax = heights.iter().cloned().chain(curve.iter().map(|p| p.1)).fold(0.0, f64::max) * 1.08;
            let fx = |x: f64| left + ((x - lo) / (hi - lo)) as f32 * (right - left);
            let fy = |y: f64| bottom - (y / ymax) as f32 * (bottom - top);
            let step = tick_step(hi - lo, 6.0);
            let mut t = (lo / step).ceil() * step;
            while t <= hi + step * 1e-9 {
                let x = fx(t);
                line(&mut s, [x, top], [x, bottom], grid, 0.8 * pt, false);
                text(&mut s, [x, bottom + 4.0 * pt], format_num(t, 4), size, ink, 0.5, 0.0);
                t += step;
            }
            let ystep = tick_step(ymax, 4.0);
            let mut y = 0.0;
            while y <= ymax {
                text(&mut s, [left - 4.0 * pt, fy(y)], format_num(y, 4), size, ink, 1.0, 0.5);
                y += ystep;
            }
            for (i, hgt) in heights.iter().enumerate() {
                rect(&mut s, fx(edges[i]) + 0.5 * pt, fy(*hgt), fx(edges[i + 1]) - 0.5 * pt, bottom, fill);
            }
            s.prims.push(Prim::Path { pts: curve.iter().map(|p| [fx(p.0), fy(p.1)]).collect(), color: Color::rgb(20, 60, 120), width: 1.6 * pt, dashed: false });
            for p in [0.025, 0.975] {
                let x = fx(quantile(&sorted, p));
                line(&mut s, [x, top], [x, bottom], red, 1.2 * pt, true);
            }
            text(&mut s, [right, 6.0 * pt], "dashed: 95% interval".into(), size, red, 1.0, 0.0);
            text(&mut s, [4.0 * pt, top - 14.0 * pt], "Density".into(), size, ink, 0.0, 0.0);
        }
        Some(_) => empty(&mut s, "Needs at least two values"),
        None if states.is_empty() => empty(&mut s, "No values"),
        None => {
            let shown = &states[..states.len().min(30)];
            let ymax = shown.iter().map(|x| x.1).max().unwrap_or(1) as f64 * 1.08;
            let bw = (right - left) / shown.len() as f32;
            let fy = |y: f64| bottom - (y / ymax) as f32 * (bottom - top);
            let step = tick_step(ymax, 5.0).max(1.0);
            let mut t = 0.0;
            while t <= ymax {
                let y = fy(t);
                line(&mut s, [left, y], [right, y], grid, 0.8 * pt, false);
                text(&mut s, [left - 4.0 * pt, y], format_num(t, 0), size, ink, 1.0, 0.5);
                t += step;
            }
            for (i, (name, n)) in shown.iter().enumerate() {
                let x0 = left + i as f32 * bw;
                rect(&mut s, x0 + bw * 0.15, fy(*n as f64), x0 + bw * 0.85, bottom, fill);
                text(&mut s, [x0 + bw / 2.0, bottom + 4.0 * pt], name.clone(), size, ink, 0.5, 0.0);
            }
            text(&mut s, [4.0 * pt, top - 14.0 * pt], "Count".into(), size, ink, 0.0, 0.0);
        }
    }
    line(&mut s, [left, bottom], [right, bottom], ink, 1.0 * pt, false);
    line(&mut s, [left, bottom], [left, top], ink, 1.0 * pt, false);
    s
}



/// A smooth curve through `pts` (x increasing): monotone cubic Hermite
/// interpolation (Fritsch & Carlson 1980), so it never overshoots between points,
/// sampled `per` times per segment.
pub fn smooth_curve(pts: &[P], per: usize) -> Vec<P> {
    let n = pts.len();
    if n < 3 {
        return pts.to_vec();
    }
    let h: Vec<f32> = (0..n - 1).map(|i| (pts[i + 1][0] - pts[i][0]).max(1e-6)).collect();
    let s: Vec<f32> = (0..n - 1).map(|i| (pts[i + 1][1] - pts[i][1]) / h[i]).collect();
    let mut m = vec![0.0f32; n];
    m[0] = s[0];
    m[n - 1] = s[n - 2];
    for i in 1..n - 1 {
        m[i] = if s[i - 1] * s[i] <= 0.0 { 0.0 } else { 3.0 * (h[i - 1] + h[i]) / ((2.0 * h[i] + h[i - 1]) / s[i - 1] + (h[i] + 2.0 * h[i - 1]) / s[i]) };
    }
    let mut out = Vec::with_capacity((n - 1) * per + 1);
    for i in 0..n - 1 {
        for k in 0..per {
            let t = k as f32 / per as f32;
            let (t2, t3) = (t * t, t * t * t);
            let y = (2.0 * t3 - 3.0 * t2 + 1.0) * pts[i][1] + (t3 - 2.0 * t2 + t) * h[i] * m[i] + (-2.0 * t3 + 3.0 * t2) * pts[i + 1][1] + (t3 - t2) * h[i] * m[i + 1];
            out.push([pts[i][0] + t * h[i], y]);
        }
    }
    out.push(pts[n - 1]);
    out
}

/// Settings of a lineages-through-time plot.
#[derive(Clone, Copy, Default)]
pub struct LttPlot<'a> {
    /// Lineages on a log scale.
    pub log: bool,
    /// Smooth curves through the steps instead of a step plot.
    pub smooth: bool,
    /// Geologic timescale rows under the plot (finest first; ages in Ma), with
    /// faint boundary lines of the finest across it.
    pub geo: &'a [crate::geotime::Level],
    /// Title at the top left.
    pub title: &'a str,
    /// The tree of the solid curve, drawn above the plot on the same time axis
    /// with dashed lines from each node down to the curve; with the age of its
    /// youngest tip and the time per branch-length unit (as for the curve).
    pub tree: Option<(&'a crate::tree::Tree, f64, f64)>,
}

/// Lineages-through-time plot, `w` × `h` pixels with sizes in points times `pt`:
/// `main` (time before present, lineages) solid over `others` drawn faint.
pub fn ltt(main: &[(f64, usize)], others: &[Vec<(f64, usize)>], opt: &LttPlot, w: f32, h: f32, pt: f32) -> Scene {
    let (log, smooth, geo, title) = (opt.log, opt.smooth, opt.geo, opt.title);
    let mut s = Scene { width: w, height: h, ..Default::default() };
    let rh = 13.0 * pt;
    let bands = geo.len() as f32 * rh;
    let (left, right, top0, bottom) = (52.0 * pt, w - 14.0 * pt, if title.is_empty() { 12.0 } else { 26.0 } * pt, h - 34.0 * pt - bands);
    // With the tree: it takes the upper part, the plot the rest.
    let tree_h = if opt.tree.is_some() { ((bottom - top0) * 0.45).max(0.0) } else { 0.0 };
    let top = if opt.tree.is_some() { top0 + tree_h + 16.0 * pt } else { top0 };
    let ink = Color::rgb(50, 50, 50);
    let grid = Color::rgb(232, 232, 232);
    let size = 9.0 * pt;
    if !title.is_empty() {
        text(&mut s, [left, 6.0 * pt], title.to_string(), 11.0 * pt, ink, 0.0, 0.0);
    }
    let all = || main.iter().chain(others.iter().flatten());
    let (old, young) = all().fold((f64::NEG_INFINITY, f64::INFINITY), |(o, y), p| (o.max(p.0), y.min(p.0)));
    let max_n = all().map(|p| p.1).max().unwrap_or(1).max(2) as f64;
    if !old.is_finite() || old <= young {
        text(&mut s, [(left + right) / 2.0, (top + bottom) / 2.0], "The tree needs branch lengths".into(), 11.0 * pt, Color::GREY, 0.5, 0.5);
        return s;
    }
    let fx = |age: f64| left + ((old - age) / (old - young)) as f32 * (right - left);
    let fy = |n: f64| {
        let t = if log { n.max(1.0).ln() / max_n.ln() } else { n / max_n };
        bottom - t as f32 * (bottom - top)
    };
    // Geologic timescale: rows under the plot, boundaries of the finest across it.
    if let Some(&finest) = geo.first() {
        use crate::scene::SceneEnv;
        for iv in crate::geotime::overlapping(finest, young, old) {
            if iv.start < old {
                let x = fx(iv.start);
                line(&mut s, [x, top], [x, bottom], Color::rgb(150, 150, 150).with_alpha(110), 0.6 * pt, false);
            }
        }
        for (r, &level) in geo.iter().enumerate() {
            let y0 = bottom + r as f32 * rh;
            for iv in crate::geotime::overlapping(level, young, old) {
                let (x0, x1) = (fx(iv.start.min(old)), fx(iv.end.max(young)));
                s.prims.push(Prim::Poly { pts: vec![[x0, y0], [x1, y0], [x1, y0 + rh], [x0, y0 + rh]], fill: Color::hex(iv.color), stroke: Some((Color::WHITE, 0.5 * pt)) });
                let ls = 7.0 * pt;
                if let Some(t) = [iv.name, iv.abbr].into_iter().find(|t| crate::scene::ApproxEnv.text_width(t, ls, false) <= (x1 - x0).abs() - 2.0 * pt) {
                    text(&mut s, [(x0 + x1) / 2.0, y0 + rh / 2.0], t.to_string(), ls, Color::BLACK, 0.5, 0.5);
                }
            }
        }
    }
    let ticks_y = bottom + bands;
    // Time ticks (ages before present).
    let step = tick_step(old - young, 6.0);
    let mut a = (young / step).ceil() * step;
    while a <= old + step * 1e-9 {
        let x = fx(a);
        line(&mut s, [x, top], [x, bottom], grid, 0.8 * pt, false);
        text(&mut s, [x, ticks_y + 4.0 * pt], format_num(a, 4), size, ink, 0.5, 0.0);
        a += step;
    }
    text(&mut s, [(left + right) / 2.0, h - 4.0 * pt], if geo.is_empty() { "Time before present" } else { "Time before present (Ma)" }.into(), size, ink, 0.5, 1.0);
    // Lineage ticks: 1, 2, 5, 10… on a log scale, round steps otherwise.
    let ticks: Vec<f64> = if log {
        let mut v = Vec::new();
        let mut m = 1.0;
        while m <= max_n * 1.0001 {
            v.extend([1.0, 2.0, 5.0].iter().map(|k| k * m).filter(|x| *x <= max_n * 1.0001));
            m *= 10.0;
        }
        v
    } else {
        let st = tick_step(max_n, 5.0).max(1.0);
        (0..).map(|i| i as f64 * st).take_while(|v| *v <= max_n * 1.0001).collect()
    };
    for n in ticks {
        let y = fy(n);
        line(&mut s, [left, y], [right, y], grid, 0.8 * pt, false);
        text(&mut s, [left - 4.0 * pt, y], format_num(n, 0), size, ink, 1.0, 0.5);
    }
    text(&mut s, [4.0 * pt, top - 14.0 * pt], "Lineages".into(), size, ink, 0.0, 0.0);
    let curve = |c: &[(f64, usize)]| -> Vec<P> {
        if smooth {
            // Through each step's new count, then a shape-preserving cubic.
            return smooth_curve(&c.iter().map(|p| [fx(p.0), fy(p.1 as f64)]).collect::<Vec<P>>(), 8);
        }
        let mut pts = Vec::with_capacity(c.len() * 2);
        for (i, p) in c.iter().enumerate() {
            if i > 0 {
                pts.push([fx(p.0), fy(c[i - 1].1 as f64)]);
            }
            pts.push([fx(p.0), fy(p.1 as f64)]);
        }
        pts
    };
    // Fainter the more curves there are, so overlaps build up into a band.
    let alpha = (3000 / others.len().max(1)).clamp(10, 160) as u8;
    for c in others {
        s.prims.push(Prim::Path { pts: curve(c), color: Color::rgb(70, 110, 170).with_alpha(alpha), width: 1.2 * pt, dashed: false });
    }
    s.prims.push(Prim::Path { pts: curve(main), color: Color::BLACK, width: 1.8 * pt, dashed: false });
    // The tree above, its nodes joined to the curve by dashed lines at their ages.
    if let Some((t, young_tip, mpu)) = opt.tree {
        let d = t.depths(true);
        let deepest = t.tips().into_iter().map(|n| d[n]).fold(0.0, f64::max);
        let age = |n: crate::tree::NodeId| young_tip + (deepest - d[n]) * mpu;
        let tips = t.tips();
        let mut y = vec![0.0f32; t.nodes.len()];
        for (i, &n) in tips.iter().enumerate() {
            y[n] = top0 + (i as f32 + 0.5) / tips.len().max(1) as f32 * tree_h;
        }
        let order = t.postorder();
        for &n in &order {
            let ch = &t.nodes[n].children;
            if let (Some(&a), Some(&b)) = (ch.first(), ch.last()) {
                y[n] = (y[a] + y[b]) / 2.0;
            }
        }
        let red = Color::rgb(213, 60, 40);
        for &n in &order {
            let ch = &t.nodes[n].children;
            if ch.is_empty() {
                continue;
            }
            let x = fx(age(n));
            // The count just after this node: the curve's step at its age.
            if let Some(p) = main.iter().min_by(|a, b| (a.0 - age(n)).abs().partial_cmp(&(b.0 - age(n)).abs()).unwrap()) {
                s.prims.push(Prim::Path { pts: vec![[x, y[n]], [x, fy(p.1 as f64)]], color: red, width: 0.9 * pt, dashed: true });
            }
        }
        let w_tree = 1.4 * pt;
        for &n in &order {
            let x = fx(age(n));
            match t.nodes[n].parent {
                Some(p) => line(&mut s, [fx(age(p)), y[n]], [x, y[n]], Color::BLACK, w_tree, false),
                None => line(&mut s, [x - 8.0 * pt, y[n]], [x, y[n]], Color::BLACK, w_tree, false),
            }
            let ch = &t.nodes[n].children;
            if let (Some(&a), Some(&b)) = (ch.first(), ch.last()) {
                line(&mut s, [x, y[a]], [x, y[b]], Color::BLACK, w_tree, false);
            }
        }
    }
    line(&mut s, [left, bottom], [right, bottom], ink, 1.0 * pt, false);
    line(&mut s, [left, bottom], [left, top], ink, 1.0 * pt, false);
    s
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ltt_scene_and_smooth_curve() {
        let main = vec![(10.0, 2), (6.0, 3), (2.0, 5), (0.0, 5)];
        let others = vec![vec![(11.0, 2), (5.0, 4), (0.0, 5)]];
        let s = ltt(&main, &others, &LttPlot { smooth: true, title: "LTT", ..Default::default() }, 400.0, 300.0, 1.0);
        // With the timescale: a band per interval of the rows asked for.
        let g = ltt(&[(30.0, 2), (0.0, 4)], &[], &LttPlot { geo: &[crate::geotime::Level::Epoch, crate::geotime::Level::Period], ..Default::default() }, 400.0, 300.0, 1.0);
        // With the tree: one dashed line per internal node, from the node to the curve.
        let t = crate::io::newick::parse_newick("((A:1,B:1):1,(C:1.5,D:1.5):0.5);").unwrap();
        let curve: Vec<(f64, usize)> = crate::ltt::ltt(&t).into_iter().map(|(x, k)| (2.0 - x, k)).collect();
        let tv = ltt(&curve, &[], &LttPlot { tree: Some((&t, 0.0, 1.0)), ..Default::default() }, 400.0, 300.0, 1.0);
        assert_eq!(tv.prims.iter().filter(|p| matches!(p, Prim::Path { dashed: true, .. })).count(), 3);
        let names: Vec<String> = g.prims.iter().filter_map(|p| if let Prim::Text { text, .. } = p { Some(text.clone()) } else { None }).collect();
        assert!(names.iter().any(|t| t == "Neogene" || t == "Ng") && names.iter().any(|t| t == "Miocene" || t == "Mio"), "{:?}", names);
        let paths: Vec<&Prim> = s.prims.iter().filter(|p| matches!(p, Prim::Path { color, .. } if *color == Color::BLACK || color.b == 170)).collect();
        assert_eq!(paths.len(), 2);
        // The smooth curve passes through every step and never overshoots.
        let pts: Vec<P> = [(0.0, 300.0), (40.0, 250.0), (45.0, 200.0), (120.0, 150.0), (200.0, 150.0), (210.0, 100.0)].iter().map(|&(x, y)| [x, y]).collect();
        let c = smooth_curve(&pts, 8);
        for p in &pts {
            assert!(c.iter().any(|q| (q[0] - p[0]).abs() < 1e-3 && (q[1] - p[1]).abs() < 1e-3));
        }
        for w in c.windows(2) {
            assert!(w[1][0] >= w[0][0] && w[1][1] <= w[0][1] + 1e-3, "{:?}", w);
        }
    }

    #[test]
    fn distribution_scenes() {
        let xs: Vec<f64> = (1..=50).map(|i| i as f64).collect();
        let s = distribution(Some(&xs), &[], "mass", 600.0, 300.0, 1.0);
        // Bars, the density curve, two dashed interval lines.
        assert!(s.prims.iter().filter(|p| matches!(p, Prim::Poly { .. })).count() >= 5);
        assert_eq!(s.prims.iter().filter(|p| matches!(p, Prim::Path { dashed: true, .. })).count(), 2);
        let c = distribution(None, &[("omni".into(), 3), ("frugi".into(), 1)], "diet", 600.0, 300.0, 1.0);
        assert_eq!(c.prims.iter().filter(|p| matches!(p, Prim::Poly { .. })).count(), 2);
        assert!(c.prims.iter().any(|p| matches!(p, Prim::Text { text, .. } if text == "omni")));
    }
}
