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


#[cfg(test)]
mod tests {
    use super::*;

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
