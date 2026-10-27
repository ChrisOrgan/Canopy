//! Statistical plots drawn as renderer-independent scenes, so the same drawing
//! appears on screen and in saved PNG / SVG / PDF files.
//!
//! The style follows Tufte (The Visual Display of Quantitative Information):
//! as much of the ink as possible shows data. No grid lines or boxes; axes are
//! range frames (a thin line spanning only the data's range, set off from it,
//! with short outward ticks); series are told apart by colour, named in their
//! own colour rather than in a legend box; fills are translucent.

use crate::io::data::{histogram, kernel_density, quantile};
use crate::scene::{Prim, Scene, P};
use crate::style::Color;
use crate::tree::format_num;

/// Ink for titles; axis lines, ticks and tick labels in 25% grey; axis titles
/// in 50% grey.
const INK: Color = Color { r: 70, g: 70, b: 70, a: 255 };
const FRAME: Color = Color { r: 191, g: 191, b: 191, a: 255 };
const LABEL: Color = Color { r: 128, g: 128, b: 128, a: 255 };
/// The one accent colour (intervals).
const ACCENT: Color = Color { r: 190, g: 50, b: 40, a: 255 };
/// Gap between the data area and its range frames.
const OFFSET: f32 = 6.0;

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

/// Round tick values within `lo..=hi`, about `target` of them.
fn ticks(lo: f64, hi: f64, target: f64) -> Vec<f64> {
    let step = tick_step(hi - lo, target);
    let mut v = Vec::new();
    let mut t = (lo / step - 1e-9).ceil() * step;
    while t <= hi + step * 1e-9 {
        v.push(if t.abs() < step * 1e-9 { 0.0 } else { t });
        t += step;
    }
    v
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

/// A horizontal range frame at `y` from `x0` to `x1`, with ticks (position,
/// label) pointing down and their labels under them.
fn x_frame(s: &mut Scene, y: f32, x0: f32, x1: f32, marks: &[(f32, String)], pt: f32) {
    line(s, [x0, y], [x1, y], FRAME, 0.7 * pt, false);
    for (x, l) in marks {
        line(s, [*x, y], [*x, y + 3.0 * pt], FRAME, 0.7 * pt, false);
        text(s, [*x, y + 5.0 * pt], l.clone(), 9.0 * pt, FRAME, 0.5, 0.0);
    }
}

/// A vertical range frame at `x` from `y0` to `y1`, ticks pointing left.
fn y_frame(s: &mut Scene, x: f32, y0: f32, y1: f32, marks: &[(f32, String)], pt: f32) {
    line(s, [x, y0], [x, y1], FRAME, 0.7 * pt, false);
    for (y, l) in marks {
        line(s, [x, *y], [x - 3.0 * pt, *y], FRAME, 0.7 * pt, false);
        text(s, [x - 5.0 * pt, *y], l.clone(), 9.0 * pt, FRAME, 1.0, 0.5);
    }
}

/// Names of the runs in their colours along `y`, from `x` (no swatches).
fn key(s: &mut Scene, names: &[String], col: &dyn Fn(usize) -> Color, x: f32, y: f32, pt: f32) {
    let size = 9.0 * pt;
    let mut x = x;
    for (i, r) in names.iter().enumerate() {
        text(s, [x, y], r.clone(), size, col(i), 0.0, 0.0);
        x += (r.chars().count() as f32 * 0.55 + 2.5) * size;
    }
}

/// The value of a sampled curve (x increasing) at `x`, linearly interpolated.
fn curve_at(c: &[(f64, f64)], x: f64) -> f64 {
    match c.iter().position(|p| p.0 >= x) {
        Some(0) => c[0].1,
        Some(i) => {
            let (a, b) = (c[i - 1], c[i]);
            a.1 + (b.1 - a.1) * ((x - a.0) / (b.0 - a.0).max(1e-300))
        }
        None => c.last().map(|p| p.1).unwrap_or(0.0),
    }
}

/// Width in points of one group in a box plot or bar chart: room for its label
/// (at least `min`), the same whatever the size of the window.
fn group_slot(names: &[String], min: f32) -> f32 {
    names.iter().map(|n| n.chars().count() as f32 * 0.55 * 9.0 + 16.0).fold(min, f32::max)
}

/// The groups of an MCMC box plot: each run, then all runs pooled.
fn box_groups(runs: &[String]) -> Vec<String> {
    let mut g = runs.to_vec();
    if runs.len() > 1 {
        g.push("All runs".into());
    }
    g
}

/// Natural width in points of a `mcmc_param` plot without a value axis across
/// (box plots), which keeps its size rather than stretching with the window;
/// None for plots that fill the width.
pub fn mcmc_width(kind: McmcPlot, runs: &[String]) -> Option<f32> {
    let g = box_groups(runs);
    (kind == McmcPlot::Box).then(|| (70.0 + g.len() as f32 * group_slot(&g, 70.0)).max(400.0))
}

/// Natural width in points of a `distribution` bar chart (categories), None
/// for histograms, which fill the width.
pub fn distribution_width(numbers: bool, states: &[(String, usize)]) -> Option<f32> {
    let names: Vec<String> = states.iter().take(30).map(|s| s.0.clone()).collect();
    (!numbers && !names.is_empty()).then(|| (66.0 + names.len() as f32 * group_slot(&names, 40.0)).max(240.0))
}

/// Distribution of one column, `w` × `h` pixels with sizes in points times `pt`:
/// for numbers, a histogram scaled to a density with a kernel density curve and
/// the central 95% interval dashed; for categories (`numbers` None), a bar chart
/// of the states (at most 30). Bars are 75% grey with white lines at the tick
/// values in place of grid lines. `density` draws the kernel density curve,
/// with the area under it shaded at opacity `shade` (0 = none). `title` goes
/// at the top left.
pub fn distribution(numbers: Option<&[f64]>, states: &[(String, usize)], density: bool, shade: f32, title: &str, w: f32, h: f32, pt: f32) -> Scene {
    let mut s = Scene { width: w, height: h, ..Default::default() };
    let (left, right, top, bottom) = (52.0 * pt, w - 14.0 * pt, 28.0 * pt, h - 34.0 * pt);
    let bar = Color::rgb(64, 64, 64);
    let curve_col = Color::rgb(86, 180, 233);
    let size = 9.0 * pt;
    text(&mut s, [left, 6.0 * pt], title.to_string(), 11.0 * pt, INK, 0.0, 0.0);
    let empty = |s: &mut Scene, msg: &str| text(s, [(left + right) / 2.0, (top + bottom) / 2.0], msg.to_string(), 11.0 * pt, Color::GREY, 0.5, 0.5);
    // White lines across the bars at the tick values.
    let erase = |s: &mut Scene, ys: &[f32], x0: f32, x1: f32| {
        for &y in ys {
            line(s, [x0, y], [x1, y], Color::WHITE, 1.0 * pt, false);
        }
    };
    match numbers {
        Some(xs) if xs.len() >= 2 => {
            let (edges, heights) = histogram(xs);
            let curve = kernel_density(xs, 200);
            let mut sorted = xs.to_vec();
            sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let (e0, e1) = (edges[0], edges[edges.len() - 1]);
            let (lo, hi) = if density { (curve[0].0.min(e0), curve[curve.len() - 1].0.max(e1)) } else { (e0, e1) };
            let top_y = heights.iter().cloned().chain(curve.iter().filter(|_| density).map(|p| p.1)).fold(0.0, f64::max);
            // Height of the bar holding `x`.
            let bar_at = |x: f64| edges.windows(2).position(|e| x >= e[0] && x <= e[1]).map(|i| heights[i]).unwrap_or(0.0);
            let ymax = top_y * 1.05;
            let fx = |x: f64| left + ((x - lo) / (hi - lo)) as f32 * (right - left);
            let fy = |y: f64| bottom - (y / ymax) as f32 * (bottom - top);
            for (i, hgt) in heights.iter().enumerate() {
                rect(&mut s, fx(edges[i]) + 0.4 * pt, fy(*hgt), fx(edges[i + 1]) - 0.4 * pt, bottom, bar);
            }
            let yt = ticks(0.0, top_y, 4.0);
            erase(&mut s, &yt.iter().skip(1).map(|&y| fy(y)).collect::<Vec<_>>(), left, right);
            if density {
                let pts: Vec<P> = curve.iter().map(|p| [fx(p.0), fy(p.1)]).collect();
                if shade > 0.0 {
                    let base: Vec<P> = pts.iter().map(|p| [p[0], bottom]).collect();
                    s.prims.push(Prim::Strip { a: pts.clone(), b: base, fill: curve_col.with_alpha((shade.clamp(0.0, 1.0) * 255.0).round() as u8), stroke: None });
                }
                s.prims.push(Prim::Path { pts, color: curve_col, width: 1.4 * pt, dashed: false });
            }
            // The central 95%: dashed from the axis up to the curve (or the bar).
            for p in [0.025, 0.975] {
                let q = quantile(&sorted, p);
                let y = if density { curve_at(&curve, q) } else { bar_at(q) };
                line(&mut s, [fx(q), bottom], [fx(q), fy(y)], ACCENT, 1.0 * pt, true);
            }
            text(&mut s, [right, 6.0 * pt], "dashed: central 95%".into(), size, ACCENT, 1.0, 0.0);
            // Range frames: the data's extent across, zero to the highest bar up.
            let (dlo, dhi) = (sorted[0], sorted[sorted.len() - 1]);
            let xt: Vec<(f32, String)> = ticks(dlo, dhi, 5.0).into_iter().map(|t| (fx(t), format_num(t, 4))).collect();
            x_frame(&mut s, bottom + OFFSET * pt, fx(dlo), fx(dhi), &xt, pt);
            let ym: Vec<(f32, String)> = yt.iter().map(|&y| (fy(y), format_num(y, 4))).collect();
            y_frame(&mut s, left - OFFSET * pt, fy(0.0), fy(yt.last().copied().unwrap_or(0.0)), &ym, pt);
            text(&mut s, [left - OFFSET * pt, top - 17.0 * pt], "Density".into(), size, LABEL, 1.0, 0.0);
        }
        Some(_) => empty(&mut s, "Needs at least two values"),
        None if states.is_empty() => empty(&mut s, "No values"),
        None => {
            let shown = &states[..states.len().min(30)];
            let most = shown.iter().map(|x| x.1).max().unwrap_or(1) as f64;
            let ymax = most * 1.05;
            let names: Vec<String> = shown.iter().map(|x| x.0.clone()).collect();
            let bw = (group_slot(&names, 40.0) * pt).min((right - left) / shown.len() as f32);
            let fy = |y: f64| bottom - (y / ymax) as f32 * (bottom - top);
            for (i, (name, n)) in shown.iter().enumerate() {
                let x0 = left + i as f32 * bw;
                rect(&mut s, x0 + bw * 0.2, fy(*n as f64), x0 + bw * 0.8, bottom, bar);
                text(&mut s, [x0 + bw / 2.0, bottom + 5.0 * pt], name.clone(), size, FRAME, 0.5, 0.0);
            }
            let step = tick_step(most, 4.0).max(1.0);
            let yt: Vec<f64> = (0..).map(|i| i as f64 * step).take_while(|v| *v <= most + 1e-9).collect();
            erase(&mut s, &yt.iter().skip(1).map(|&y| fy(y)).collect::<Vec<_>>(), left, right);
            let ym: Vec<(f32, String)> = yt.iter().map(|&y| (fy(y), format_num(y, 0))).collect();
            y_frame(&mut s, left - OFFSET * pt, fy(0.0), fy(yt.last().copied().unwrap_or(0.0)), &ym, pt);
            text(&mut s, [left - OFFSET * pt, top - 17.0 * pt], "Count".into(), size, LABEL, 1.0, 0.0);
        }
    }
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
    let gap = if geo.is_empty() { 0.0 } else { 3.0 * pt };
    let (left, right, top0, bottom) = (52.0 * pt, w - 14.0 * pt, if title.is_empty() { 14.0 } else { 28.0 } * pt, h - 38.0 * pt - bands - gap);
    // With the tree: it takes the upper part, the plot the rest.
    let tree_h = if opt.tree.is_some() { ((bottom - top0) * 0.45).max(0.0) } else { 0.0 };
    let top = if opt.tree.is_some() { top0 + tree_h + 16.0 * pt } else { top0 };
    let size = 9.0 * pt;
    if !title.is_empty() {
        text(&mut s, [left, 6.0 * pt], title.to_string(), 11.0 * pt, INK, 0.0, 0.0);
    }
    let all = || main.iter().chain(others.iter().flatten());
    let (old, young) = all().fold((f64::NEG_INFINITY, f64::INFINITY), |(o, y), p| (o.max(p.0), y.min(p.0)));
    let max_n = all().map(|p| p.1).max().unwrap_or(1).max(2) as f64;
    let min_n = all().map(|p| p.1).min().unwrap_or(1).max(1) as f64;
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
                line(&mut s, [x, top], [x, bottom], Color::rgb(150, 150, 150).with_alpha(90), 0.5 * pt, false);
            }
        }
        for (r, &level) in geo.iter().enumerate() {
            let y0 = bottom + gap + r as f32 * rh;
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
    // Time: a range frame from the oldest to the youngest age, under the timescale.
    let xt: Vec<(f32, String)> = ticks(young, old, 6.0).into_iter().map(|a| (fx(a), format_num(a, 4))).collect();
    x_frame(&mut s, bottom + gap + bands + OFFSET * pt, fx(old), fx(young), &xt, pt);
    text(&mut s, [(left + right) / 2.0, h - 4.0 * pt], if geo.is_empty() { "Time before present" } else { "Time before present (Ma)" }.into(), size, LABEL, 0.5, 1.0);
    // Lineages: 1, 2, 5, 10… on a log scale, round steps otherwise; the frame
    // spans the fewest to the most lineages.
    let yt: Vec<f64> = if log {
        let mut v = Vec::new();
        let mut m = 1.0;
        while m <= max_n * 1.0001 {
            v.extend([1.0, 2.0, 5.0].iter().map(|k| k * m).filter(|x| *x >= min_n * 0.9999 && *x <= max_n * 1.0001));
            m *= 10.0;
        }
        v
    } else {
        // Round steps inside the range, plus the fewest and the most lineages.
        let st = tick_step(max_n - min_n, 5.0).max(1.0);
        let mut v: Vec<f64> = (0..).map(|i| i as f64 * st).take_while(|v| *v <= max_n).filter(|v| *v > min_n + st * 0.4 && *v < max_n - st * 0.4).collect();
        v.insert(0, min_n);
        v.push(max_n);
        v
    };
    let ym: Vec<(f32, String)> = yt.iter().map(|&n| (fy(n), format_num(n, 0))).collect();
    y_frame(&mut s, left - OFFSET * pt, fy(min_n), fy(max_n), &ym, pt);
    text(&mut s, [left - OFFSET * pt, top - 17.0 * pt], "Lineages".into(), size, LABEL, 1.0, 0.0);
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
    let alpha = (3000 / others.len().max(1)).clamp(10, 140) as u8;
    for c in others {
        s.prims.push(Prim::Path { pts: curve(c), color: Color::rgb(70, 110, 170).with_alpha(alpha), width: 0.9 * pt, dashed: false });
    }
    s.prims.push(Prim::Path { pts: curve(main), color: Color::BLACK, width: 1.4 * pt, dashed: false });
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
        for &n in &order {
            let ch = &t.nodes[n].children;
            if ch.is_empty() {
                continue;
            }
            let x = fx(age(n));
            // The count just after this node: the curve's step at its age.
            if let Some(p) = main.iter().min_by(|a, b| (a.0 - age(n)).abs().partial_cmp(&(b.0 - age(n)).abs()).unwrap()) {
                s.prims.push(Prim::Path { pts: vec![[x, y[n]], [x, fy(p.1 as f64)]], color: ACCENT.with_alpha(200), width: 0.7 * pt, dashed: true });
            }
        }
        let w_tree = 1.1 * pt;
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
    s
}

/// The node-density (delta) test plot, `w` × `h` pixels with sizes in points
/// times `pt`: nodes along each root-to-tip path against its length, one dot
/// per species, with the fitted curve n = βx^δ.
pub fn node_density(r: &crate::nodedensity::DeltaTest, w: f32, h: f32, pt: f32) -> Scene {
    let mut s = Scene { width: w, height: h, ..Default::default() };
    let (left, right, top, bottom) = (48.0 * pt, w - 14.0 * pt, 22.0 * pt, h - 38.0 * pt);
    let size = 9.0 * pt;
    let fold = |v: &[f64]| v.iter().filter(|x| x.is_finite()).fold((f64::INFINITY, f64::NEG_INFINITY), |a, &x| (a.0.min(x), a.1.max(x)));
    let ((xlo, xhi), (ylo, yhi)) = (fold(&r.path_lengths), fold(&r.nodes));
    if !xlo.is_finite() || !ylo.is_finite() {
        text(&mut s, [(left + right) / 2.0, (top + bottom) / 2.0], "No species".into(), 11.0 * pt, Color::GREY, 0.5, 0.5);
        return s;
    }
    // Equal path lengths (ultrametric trees; differences below 0.1% are rounding)
    // and equal node counts are widened, so the dots stack in the middle.
    let widen = |lo: f64, hi: f64| {
        let mid = (lo + hi) / 2.0;
        if hi - lo > 1e-3 * mid.abs().max(1e-9) { (lo, hi) } else { (mid - mid.abs().max(1.0) * 0.1, mid + mid.abs().max(1.0) * 0.1) }
    };
    let ((xlo, xhi), (ylo, yhi)) = (widen(xlo, xhi), if yhi > ylo { (ylo, yhi) } else { (ylo - 1.0, yhi + 1.0) });
    let fx = |x: f64| left + ((x - xlo) / (xhi - xlo)) as f32 * (right - left);
    let fy = |y: f64| bottom - ((y - ylo) / (yhi - ylo)) as f32 * (bottom - top);
    // The fitted curve over the data's range, clipped to it.
    let curve: Vec<P> = (0..=100)
        .map(|i| xlo + (xhi - xlo) * i as f64 / 100.0)
        .map(|x| (x, r.beta * x.powf(r.delta)))
        .filter(|(_, y)| y.is_finite() && *y >= ylo && *y <= yhi)
        .map(|(x, y)| [fx(x), fy(y)])
        .collect();
    if curve.len() > 1 {
        s.prims.push(Prim::Path { pts: curve, color: Color::rgb(214, 96, 0), width: 1.3 * pt, dashed: false });
    }
    for (x, y) in r.path_lengths.iter().zip(&r.nodes) {
        if x.is_finite() && y.is_finite() {
            s.prims.push(Prim::Circle { c: [fx(*x), fy(*y)], r: 2.0 * pt, fill: Color::rgb(30, 95, 191).with_alpha(190), stroke: None });
        }
    }
    let xt: Vec<(f32, String)> = ticks(xlo, xhi, 5.0).into_iter().map(|t| (fx(t), format_num(t, 3))).collect();
    x_frame(&mut s, bottom + OFFSET * pt, fx(xlo), fx(xhi), &xt, pt);
    let yt: Vec<(f32, String)> = ticks(ylo, yhi, 5.0).into_iter().filter(|t| t.fract() == 0.0).map(|t| (fy(t), format_num(t, 0))).collect();
    y_frame(&mut s, left - OFFSET * pt, fy(ylo), fy(yhi), &yt, pt);
    text(&mut s, [(left + right) / 2.0, h - 4.0 * pt], "Root-to-tip path length".into(), size, LABEL, 0.5, 1.0);
    text(&mut s, [left - OFFSET * pt, top - 17.0 * pt], "Nodes".into(), size, LABEL, 1.0, 0.0);
    s
}

/// How an MCMC parameter is drawn by `mcmc_param`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum McmcPlot {
    /// A kernel density curve per run.
    Density,
    /// A box and whiskers per run, and one for all runs pooled.
    Box,
    /// Each run's samples in order (a trace), runs overlaid.
    Trace,
}

/// Whether samples are discrete (counts, states): every value a whole number,
/// and at most 50 distinct values. Plotted as bars instead of a density curve.
pub fn is_discrete(chains: &[Vec<f64>]) -> bool {
    let mut seen: Vec<f64> = Vec::new();
    for &x in chains.iter().flatten().filter(|x| x.is_finite()) {
        if x.fract() != 0.0 {
            return false;
        }
        if !seen.contains(&x) {
            if seen.len() == 50 {
                return false;
            }
            seen.push(x);
        }
    }
    !seen.is_empty()
}

/// Default opacity of the shading under density curves.
pub const SHADE_OPACITY: f32 = 0.25;

/// One MCMC parameter's samples, run by run (`chains`, named `runs`), as density
/// curves, box plots or traces, `w` × `h` pixels with sizes in points times
/// `pt`. Densities: the area under each curve is shaded in its run's colour at
/// opacity `shade` (0 = none); no density axis, since only the shape matters.
/// Box plots in Tufte's form: whiskers (to the most extreme samples within
/// 1.5 IQR) as lines, a gap for the quartiles, a dot at the median, and the 95%
/// HPD as a line beside each. `title` at the top left.
pub fn mcmc_param(chains: &[Vec<f64>], runs: &[String], kind: McmcPlot, shade: f32, title: &str, w: f32, h: f32, pt: f32) -> Scene {
    let mut s = Scene { width: w, height: h, ..Default::default() };
    let size = 9.0 * pt;
    let (left, right, top, bottom) = (56.0 * pt, w - 14.0 * pt, 40.0 * pt, h - 36.0 * pt);
    text(&mut s, [left, 6.0 * pt], title.to_string(), 11.0 * pt, INK, 0.0, 0.0);
    let col = |i: usize| crate::style::Palette::OkabeIto.discrete(i, runs.len().max(2));
    let clean: Vec<Vec<f64>> = chains.iter().map(|c| c.iter().copied().filter(|x| x.is_finite()).collect()).collect();
    let all: Vec<f64> = clean.iter().flatten().copied().collect();
    let (lo, hi) = all.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |a, &x| (a.0.min(x), a.1.max(x)));
    if !lo.is_finite() {
        text(&mut s, [(left + right) / 2.0, (top + bottom) / 2.0], "No samples".into(), 11.0 * pt, Color::GREY, 0.5, 0.5);
        return s;
    }
    if hi <= lo {
        text(&mut s, [(left + right) / 2.0, (top + bottom) / 2.0], format!("Every sample is {}", format_num(lo, 4)), 11.0 * pt, Color::GREY, 0.5, 0.5);
        return s;
    }
    // The runs, named in their colours (the box plot names its groups below).
    if kind != McmcPlot::Box && runs.len() > 1 {
        key(&mut s, runs, &col, left, 22.0 * pt, pt);
    }
    // A value axis (range frame over the samples) on the left.
    let value_frame = |s: &mut Scene, fy: &dyn Fn(f64) -> f32| {
        let m: Vec<(f32, String)> = ticks(lo, hi, 5.0).into_iter().map(|t| (fy(t), format_num(t, 4))).collect();
        y_frame(s, left - OFFSET * pt, fy(lo), fy(hi), &m, pt);
    };
    match kind {
        McmcPlot::Density if is_discrete(&clean) => {
            // Whole-number values (counts): a bar per value and run, the
            // proportion of that run's samples, with white lines in place of a grid.
            let values: Vec<f64> = {
                let mut v: Vec<f64> = all.clone();
                v.sort_by(|a, b| a.partial_cmp(b).unwrap());
                v.dedup();
                v
            };
            let props: Vec<Vec<f64>> = clean.iter().map(|c| values.iter().map(|&v| c.iter().filter(|&&x| x == v).count() as f64 / c.len().max(1) as f64).collect()).collect();
            let pmax = props.iter().flatten().fold(0.0f64, |a, &b| a.max(b));
            let (xlo, xhi) = (lo - 0.5, hi + 0.5);
            let fx = |v: f64| left + ((v - xlo) / (xhi - xlo)) as f32 * (right - left);
            let fy = |p: f64| bottom - (p / (pmax * 1.05).max(1e-300)) as f32 * (bottom - top);
            let k = clean.len().max(1) as f32;
            let unit = (fx(1.0) - fx(0.0)).abs();
            let bw = (unit * 0.8 / k).min(40.0 * pt);
            for (j, &v) in values.iter().enumerate() {
                for (i, p) in props.iter().enumerate() {
                    let x0 = fx(v) - bw * k / 2.0 + bw * i as f32;
                    s.prims.push(Prim::Poly { pts: vec![[x0, fy(p[j])], [x0 + bw * 0.9, fy(p[j])], [x0 + bw * 0.9, bottom], [x0, bottom]], fill: col(i), stroke: None });
                }
            }
            let yt = ticks(0.0, pmax, 4.0);
            for &y in yt.iter().skip(1) {
                line(&mut s, [left, fy(y)], [right, fy(y)], Color::WHITE, 1.0 * pt, false);
            }
            let ym: Vec<(f32, String)> = yt.iter().map(|&y| (fy(y), format_num(y, 3))).collect();
            y_frame(&mut s, left - OFFSET * pt, fy(0.0), fy(yt.last().copied().unwrap_or(0.0)), &ym, pt);
            // Above the axis, clear of the run names along the top.
            text(&mut s, [left - OFFSET * pt, fy(yt.last().copied().unwrap_or(0.0)) - 14.0 * pt], "Proportion".into(), size, LABEL, 1.0, 0.5);
            // Every value labelled when there are few, round steps otherwise.
            let xt: Vec<f64> = if hi - lo <= 20.0 { values.clone() } else { ticks(lo, hi, 6.0).into_iter().filter(|t| t.fract() == 0.0).collect() };
            let m: Vec<(f32, String)> = xt.into_iter().map(|t| (fx(t), format_num(t, 0))).collect();
            x_frame(&mut s, bottom + OFFSET * pt, fx(lo), fx(hi), &m, pt);
        }
        McmcPlot::Density => {
            let curves: Vec<Vec<(f64, f64)>> = clean.iter().map(|c| if c.len() > 1 { crate::io::data::kernel_density(c, 200) } else { Vec::new() }).collect();
            let (xlo, xhi) = curves.iter().flatten().fold((lo, hi), |a, p| (a.0.min(p.0), a.1.max(p.0)));
            let ymax = curves.iter().flatten().map(|p| p.1).fold(0.0, f64::max) * 1.05;
            let fx = |v: f64| left + ((v - xlo) / (xhi - xlo)) as f32 * (right - left);
            let fy = |v: f64| bottom - (v / ymax.max(1e-300)) as f32 * (bottom - top);
            // Translucent fills first, so every curve stays visible where runs overlap.
            if shade > 0.0 {
                let a = (shade.clamp(0.0, 1.0) * 255.0).round() as u8;
                for (i, c) in curves.iter().enumerate().filter(|(_, c)| c.len() > 1) {
                    let top_pts: Vec<P> = c.iter().map(|p| [fx(p.0), fy(p.1)]).collect();
                    let base: Vec<P> = top_pts.iter().map(|p| [p[0], bottom]).collect();
                    s.prims.push(Prim::Strip { a: top_pts, b: base, fill: col(i).with_alpha(a), stroke: None });
                }
            }
            for (i, c) in curves.iter().enumerate() {
                if c.len() > 1 {
                    s.prims.push(Prim::Path { pts: c.iter().map(|p| [fx(p.0), fy(p.1)]).collect(), color: col(i), width: 1.2 * pt, dashed: false });
                }
            }
            let m: Vec<(f32, String)> = ticks(lo, hi, 6.0).into_iter().map(|t| (fx(t), format_num(t, 4))).collect();
            x_frame(&mut s, bottom + OFFSET * pt, fx(lo), fx(hi), &m, pt);
        }
        McmcPlot::Trace => {
            // Sample number across, value up; well-mixed runs look like overlapping noise.
            let longest = clean.iter().map(Vec::len).max().unwrap_or(1).max(2);
            let fx = |i: usize| left + (i as f64 / (longest - 1) as f64) as f32 * (right - left);
            let fy = |v: f64| bottom - ((v - lo) / (hi - lo)) as f32 * (bottom - top);
            for (i, c) in clean.iter().enumerate() {
                if c.len() > 1 {
                    s.prims.push(Prim::Path { pts: c.iter().enumerate().map(|(j, &v)| [fx(j), fy(v)]).collect(), color: col(i).with_alpha(170), width: 0.6 * pt, dashed: false });
                }
            }
            value_frame(&mut s, &fy);
            let m: Vec<(f32, String)> = ticks(0.0, (longest - 1) as f64, 6.0).into_iter().map(|t| (fx(t as usize), format_num(t, 0))).collect();
            x_frame(&mut s, bottom + OFFSET * pt, fx(0), fx(longest - 1), &m, pt);
            text(&mut s, [(left + right) / 2.0, h - 4.0 * pt], "Sample".into(), size, LABEL, 0.5, 1.0);
        }
        McmcPlot::Box => {
            // One per run, then all runs pooled.
            let mut groups: Vec<(String, Vec<f64>, Color)> = clean.iter().enumerate().map(|(i, c)| (runs.get(i).cloned().unwrap_or_default(), c.clone(), col(i))).collect();
            if clean.len() > 1 {
                groups.push(("All runs".into(), all.clone(), Color::rgb(90, 90, 90)));
            }
            let fy = |v: f64| bottom - ((v - lo) / (hi - lo)) as f32 * (bottom - top);
            let names: Vec<String> = groups.iter().map(|g| g.0.clone()).collect();
            let slot = (group_slot(&names, 70.0) * pt).min((right - left) / groups.len() as f32);
            for (g, (name, v, c)) in groups.iter().enumerate() {
                if v.is_empty() {
                    continue;
                }
                let mut v = v.clone();
                v.sort_by(|a, b| a.partial_cmp(b).unwrap());
                let (q1, med, q3) = (quantile(&v, 0.25), quantile(&v, 0.5), quantile(&v, 0.75));
                let iqr = q3 - q1;
                let wlo = v.iter().copied().find(|&x| x >= q1 - 1.5 * iqr).unwrap_or(v[0]);
                let whi = v.iter().rev().copied().find(|&x| x <= q3 + 1.5 * iqr).unwrap_or(v[v.len() - 1]);
                let cx = left + slot * (g as f32 + 0.5);
                // Whiskers, the gap between the quartiles, and the median.
                line(&mut s, [cx, fy(wlo)], [cx, fy(q1)], *c, 1.2 * pt, false);
                line(&mut s, [cx, fy(q3)], [cx, fy(whi)], *c, 1.2 * pt, false);
                s.prims.push(Prim::Circle { c: [cx, fy(med)], r: 2.6 * pt, fill: *c, stroke: None });
                // Outliers (at most 200, spread evenly) as small faint dots.
                let out: Vec<f64> = v.iter().copied().filter(|&x| x < wlo || x > whi).collect();
                let every = out.len().div_ceil(200).max(1);
                for &o in out.iter().step_by(every) {
                    s.prims.push(Prim::Circle { c: [cx, fy(o)], r: 1.0 * pt, fill: c.with_alpha(140), stroke: None });
                }
                if let Some((a, b)) = crate::consensus::hpd(&v, 0.95) {
                    line(&mut s, [cx + 7.0 * pt, fy(a)], [cx + 7.0 * pt, fy(b)], ACCENT, 1.2 * pt, false);
                }
                text(&mut s, [cx, bottom + OFFSET * pt], name.clone(), size, *c, 0.5, 0.0);
            }
            value_frame(&mut s, &fy);
            text(&mut s, [right, 22.0 * pt], "lines: whiskers to 1.5 IQR · gap: quartiles · dot: median · red: 95% HPD".into(), size, LABEL, 1.0, 0.0);
        }
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Grid lines are gone: no full-width or full-height light lines.
    fn no_grid(s: &Scene) -> bool {
        !s.prims.iter().any(|p| matches!(p, Prim::Path { color, .. } if *color == Color::rgb(232, 232, 232)))
    }

    #[test]
    fn mcmc_parameter_plots() {
        let runs = vec!["a".to_string(), "b".to_string()];
        let chains = vec![vec![1.0, 2.0, 3.0, 2.5, 1.5], vec![1.2, 2.2, 2.8, 2.4, 9.0]];
        // Density: one curve per run, each shaded at the opacity asked for.
        let d = mcmc_param(&chains, &runs, McmcPlot::Density, 0.3, "x", 400.0, 300.0, 1.0);
        assert_eq!(d.prims.iter().filter(|p| matches!(p, Prim::Path { pts, .. } if pts.len() > 50)).count(), 2);
        let fills: Vec<u8> = d.prims.iter().filter_map(|p| if let Prim::Strip { fill, .. } = p { Some(fill.a) } else { None }).collect();
        assert_eq!(fills, vec![77, 77]);
        let plain = mcmc_param(&chains, &runs, McmcPlot::Density, 0.0, "x", 400.0, 300.0, 1.0);
        assert!(!plain.prims.iter().any(|p| matches!(p, Prim::Strip { .. })) && no_grid(&d));
        // Box plots: a median dot per run plus all runs; run b's 9.0 is an outlier.
        let b = mcmc_param(&chains, &runs, McmcPlot::Box, 0.3, "x", 400.0, 300.0, 1.0);
        assert_eq!(b.prims.iter().filter(|p| matches!(p, Prim::Circle { r, .. } if *r > 2.0)).count(), 3);
        assert!(b.prims.iter().any(|p| matches!(p, Prim::Circle { r, .. } if *r < 2.0)));
        assert!(!b.prims.iter().any(|p| matches!(p, Prim::Poly { .. })) && no_grid(&b));
        let labels: Vec<String> = b.prims.iter().filter_map(|p| if let Prim::Text { text, .. } = p { Some(text.clone()) } else { None }).collect();
        assert!(labels.contains(&"All runs".to_string()) && !labels.iter().any(|t| t == "-0"));
    }

    #[test]
    fn plots_without_an_x_axis_keep_their_width() {
        let runs = vec!["a".to_string(), "b".to_string()];
        let chains = vec![vec![1.0, 2.0, 3.0], vec![1.5, 2.5, 3.5]];
        let medians = |w: f32| -> Vec<f32> { mcmc_param(&chains, &runs, McmcPlot::Box, 0.0, "x", w, 300.0, 1.0).prims.iter().filter_map(|p| if let Prim::Circle { c, r, .. } = p { (*r > 2.0).then_some(c[0]) } else { None }).collect() };
        let natural = mcmc_width(McmcPlot::Box, &runs).unwrap();
        assert_eq!(medians(natural), medians(natural * 3.0));
        assert!(mcmc_width(McmcPlot::Density, &runs).is_none());
        let states = vec![("x".to_string(), 2), ("y".to_string(), 1)];
        let bars = |w: f32| -> Vec<f32> { distribution(None, &states, true, 0.25, "t", w, 300.0, 1.0).prims.iter().filter_map(|p| if let Prim::Poly { pts, .. } = p { Some(pts[0][0]) } else { None }).collect() };
        let natural = distribution_width(false, &states).unwrap();
        assert_eq!(bars(natural), bars(natural * 3.0));
        assert!(distribution_width(true, &states).is_none());
    }

    #[test]
    fn discrete_parameters_get_bars() {
        let runs = vec!["a".to_string(), "b".to_string()];
        let counts = vec![vec![0.0, 1.0, 1.0, 2.0], vec![1.0, 1.0, 2.0, 3.0]];
        assert!(is_discrete(&counts) && !is_discrete(&[vec![0.5, 1.0]]));
        let s = mcmc_param(&counts, &runs, McmcPlot::Density, 0.25, "x", 400.0, 300.0, 1.0);
        // Four values (0 to 3), a bar for each run: proportions, no density curve.
        assert_eq!(s.prims.iter().filter(|p| matches!(p, Prim::Poly { .. })).count(), 8);
        assert!(!s.prims.iter().any(|p| matches!(p, Prim::Strip { .. })) && no_grid(&s));
    }

    #[test]
    fn range_frames_span_the_data() {
        // The value axis of a trace runs from the smallest to the largest sample.
        let chains = vec![vec![2.0, 4.0, 3.0, 5.0]];
        let s = mcmc_param(&chains, &["a".to_string()], McmcPlot::Trace, 0.0, "x", 400.0, 300.0, 1.0);
        let (top, bottom) = (40.0, 300.0 - 36.0);
        let x = 56.0 - OFFSET;
        let frame = s.prims.iter().find_map(|p| match p { Prim::Path { pts, color, .. } if *color == FRAME && pts[0][0] == x && pts[1][0] == x => Some((pts[0][1], pts[1][1])), _ => None }).unwrap();
        assert!((frame.0 - bottom).abs() < 0.01 && (frame.1 - top).abs() < 0.01, "{:?}", frame);
    }

    #[test]
    fn node_density_plot_has_offset_grey_axes() {
        let t = crate::io::newick::parse_newick("(((A:1,B:2):1,(C:1,D:3):2):1,((E:2,F:1):1,G:4):1);").unwrap();
        let r = crate::nodedensity::delta_test(&t).unwrap();
        let s = node_density(&r, 400.0, 300.0, 1.0);
        assert_eq!(s.prims.iter().filter(|p| matches!(p, Prim::Circle { .. })).count(), 7);
        // Range frames set off from the data, in 25% grey.
        assert_eq!(FRAME, Color::rgb(191, 191, 191));
        let x_axis = 48.0 - OFFSET;
        assert!(s.prims.iter().any(|p| matches!(p, Prim::Path { pts, color, .. } if *color == FRAME && pts[0][0] == x_axis && pts[1][0] == x_axis)));
        assert!(no_grid(&s));
    }

    #[test]
    fn ltt_scene_and_smooth_curve() {
        let main = vec![(10.0, 2), (6.0, 3), (2.0, 5), (0.0, 5)];
        let others = vec![vec![(11.0, 2), (5.0, 4), (0.0, 5)]];
        let s = ltt(&main, &others, &LttPlot { smooth: true, title: "LTT", ..Default::default() }, 400.0, 300.0, 1.0);
        assert!(no_grid(&s));
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
        let s = distribution(Some(&xs), &[], true, 0.25, "mass", 600.0, 300.0, 1.0);
        // Bars, the density curve, two dashed interval lines, white lines in place of a grid.
        assert!(s.prims.iter().filter(|p| matches!(p, Prim::Poly { .. })).count() >= 5);
        assert_eq!(s.prims.iter().filter(|p| matches!(p, Prim::Path { dashed: true, .. })).count(), 2);
        assert!(s.prims.iter().any(|p| matches!(p, Prim::Path { color, .. } if *color == Color::WHITE)) && no_grid(&s));
        // Bars in 75% grey; the curve shaded at the opacity asked for, or left out.
        assert!(s.prims.iter().any(|p| matches!(p, Prim::Poly { fill, .. } if *fill == Color::rgb(64, 64, 64))));
        assert!(s.prims.iter().any(|p| matches!(p, Prim::Strip { fill, .. } if fill.a == 64)));
        let bare = distribution(Some(&xs), &[], false, 0.25, "mass", 600.0, 300.0, 1.0);
        assert!(!bare.prims.iter().any(|p| matches!(p, Prim::Strip { .. })));
        assert!(!bare.prims.iter().any(|p| matches!(p, Prim::Path { pts, .. } if pts.len() > 50)));
        assert_eq!(bare.prims.iter().filter(|p| matches!(p, Prim::Path { dashed: true, .. })).count(), 2);
        let c = distribution(None, &[("omni".into(), 3), ("frugi".into(), 1)], true, 0.25, "diet", 600.0, 300.0, 1.0);
        assert_eq!(c.prims.iter().filter(|p| matches!(p, Prim::Poly { .. })).count(), 2);
        assert!(c.prims.iter().any(|p| matches!(p, Prim::Text { text, .. } if text == "omni")));
    }
}


