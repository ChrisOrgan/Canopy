//! The geologic time scale (International Chronostratigraphic Chart, v2024/12):
//! eras, periods, epochs and ages (stages) of the Phanerozoic, with their standard
//! CGMW colors, for timescale bars under time-calibrated trees (like deeptime's
//! `coord_geo`). Ages in Ma. Names, boundaries and colors from Macrostrat
//! (macrostrat.org, "international" timescales, CC-BY 4.0).

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Level {
    Era,
    Period,
    Epoch,
    /// Ages (stages), the finest level.
    Age,
}

pub struct Interval {
    pub name: &'static str,
    pub abbr: &'static str,
    /// Older boundary (Ma).
    pub start: f64,
    /// Younger boundary (Ma).
    pub end: f64,
    pub color: &'static str,
    pub level: Level,
}

const fn iv(name: &'static str, abbr: &'static str, start: f64, end: f64, color: &'static str, level: Level) -> Interval {
    Interval { name, abbr, start, end, color, level }
}

use Level::*;

pub const INTERVALS: &[Interval] = &[
    // Eras
    iv("Cenozoic", "Cz", 66.0, 0.0, "#F2F91D", Era),
    iv("Mesozoic", "Mz", 251.902, 66.0, "#67C5CA", Era),
    iv("Paleozoic", "Pz", 538.8, 251.902, "#99C08D", Era),
    // Periods
    iv("Quaternary", "Q", 2.58, 0.0, "#F9F97F", Period),
    iv("Neogene", "Ng", 23.04, 2.58, "#FFE619", Period),
    iv("Paleogene", "Pg", 66.0, 23.04, "#FD9A52", Period),
    iv("Cretaceous", "K", 143.1, 66.0, "#7FC64E", Period),
    iv("Jurassic", "J", 201.4, 143.1, "#34B2C9", Period),
    iv("Triassic", "Tr", 251.902, 201.4, "#812B92", Period),
    iv("Permian", "P", 298.9, 251.902, "#F04028", Period),
    iv("Carboniferous", "C", 358.86, 298.9, "#67A599", Period),
    iv("Devonian", "D", 419.62, 358.86, "#CB8C37", Period),
    iv("Silurian", "S", 443.1, 419.62, "#B3E1B6", Period),
    iv("Ordovician", "O", 486.85, 443.1, "#009270", Period),
    iv("Cambrian", "Cm", 538.8, 486.85, "#7FA056", Period),
    // Epochs
    iv("Holocene", "Hol", 0.0117, 0.0, "#FEF2E0", Epoch),
    iv("Pleistocene", "Ple", 2.58, 0.0117, "#FFF2AE", Epoch),
    iv("Pliocene", "Pli", 5.333, 2.58, "#FFFF99", Epoch),
    iv("Miocene", "Mio", 23.04, 5.333, "#FFFF00", Epoch),
    iv("Oligocene", "Ol", 33.9, 23.04, "#FDC07A", Epoch),
    iv("Eocene", "Eo", 56.0, 33.9, "#FDB46C", Epoch),
    iv("Paleocene", "Pal", 66.0, 56.0, "#FDA75F", Epoch),
    iv("Late Cretaceous", "LK", 100.5, 66.0, "#A6D84A", Epoch),
    iv("Early Cretaceous", "EK", 143.1, 100.5, "#8CCD57", Epoch),
    iv("Late Jurassic", "LJ", 161.5, 143.1, "#B3E3EE", Epoch),
    iv("Middle Jurassic", "MJ", 174.7, 161.5, "#80CFD8", Epoch),
    iv("Early Jurassic", "EJ", 201.4, 174.7, "#42AED0", Epoch),
    iv("Late Triassic", "LTr", 237.0, 201.4, "#BD8CC3", Epoch),
    iv("Middle Triassic", "MTr", 246.7, 237.0, "#B168B1", Epoch),
    iv("Early Triassic", "ETr", 251.902, 246.7, "#983999", Epoch),
    iv("Lopingian", "Lop", 259.51, 251.902, "#FBA794", Epoch),
    iv("Guadalupian", "Gua", 274.4, 259.51, "#FB745C", Epoch),
    iv("Cisuralian", "Cis", 298.9, 274.4, "#EF5845", Epoch),
    iv("Pennsylvanian", "Penn", 323.4, 298.9, "#99C2B5", Epoch),
    iv("Mississippian", "Miss", 358.86, 323.4, "#678F66", Epoch),
    iv("Late Devonian", "LD", 382.31, 358.86, "#F1E19D", Epoch),
    iv("Middle Devonian", "MD", 393.47, 382.31, "#F1C868", Epoch),
    iv("Early Devonian", "ED", 419.62, 393.47, "#E5AC4D", Epoch),
    iv("Pridoli", "Pri", 422.7, 419.62, "#E6F5E1", Epoch),
    iv("Ludlow", "Lud", 426.7, 422.7, "#BFE6CF", Epoch),
    iv("Wenlock", "Wen", 432.9, 426.7, "#B3E1C2", Epoch),
    iv("Llandovery", "Lla", 443.1, 432.9, "#99D7B3", Epoch),
    iv("Late Ordovician", "LO", 458.2, 443.1, "#7FCA93", Epoch),
    iv("Middle Ordovician", "MO", 471.3, 458.2, "#4DB47E", Epoch),
    iv("Early Ordovician", "EO", 486.85, 471.3, "#1A9D6F", Epoch),
    iv("Furongian", "Fur", 497.0, 486.85, "#B3E095", Epoch),
    iv("Miaolingian", "Mia", 506.5, 497.0, "#A6CF86", Epoch),
    iv("Series 2", "S2", 521.0, 506.5, "#99C078", Epoch),
    iv("Terreneuvian", "Ter", 538.8, 521.0, "#8CB06C", Epoch),
    // Ages
    iv("Meghalayan", "Meg", 0.0042, 0.0, "#FEF2E0", Age),
    iv("Northgrippian", "Nor", 0.0082, 0.0042, "#FEF2E0", Age),
    iv("Greenlandian", "Gre", 0.0117, 0.0082, "#FEF2E0", Age),
    iv("Late Pleistocene", "LP", 0.129, 0.0117, "#FFF2C7", Age),
    iv("Chibanian", "Chi", 0.774, 0.129, "#FFF2C7", Age),
    iv("Calabrian", "Cal", 1.8, 0.774, "#FFF2C7", Age),
    iv("Gelasian", "Gel", 2.58, 1.8, "#FFEDB3", Age),
    iv("Piacenzian", "Pia", 3.6, 2.58, "#FFFFBF", Age),
    iv("Zanclean", "Zan", 5.333, 3.6, "#FFFFB3", Age),
    iv("Messinian", "Mes", 7.246, 5.333, "#FFFF73", Age),
    iv("Tortonian", "Tor", 11.63, 7.246, "#FFFF66", Age),
    iv("Serravallian", "Ser", 13.82, 11.63, "#FFFF59", Age),
    iv("Langhian", "Lan", 15.98, 13.82, "#FFFF4D", Age),
    iv("Burdigalian", "Bur", 20.45, 15.98, "#FFFF41", Age),
    iv("Aquitanian", "Aqu", 23.04, 20.45, "#FFFF33", Age),
    iv("Chattian", "Cha", 27.3, 23.04, "#FEE6AA", Age),
    iv("Rupelian", "Rup", 33.9, 27.3, "#FED99A", Age),
    iv("Priabonian", "Pri", 37.71, 33.9, "#FDCDA1", Age),
    iv("Bartonian", "Bar", 41.03, 37.71, "#FDC091", Age),
    iv("Lutetian", "Lut", 48.07, 41.03, "#FCB482", Age),
    iv("Ypresian", "Ypr", 56.0, 48.07, "#FCA773", Age),
    iv("Thanetian", "Tha", 59.24, 56.0, "#FDBF6F", Age),
    iv("Selandian", "Sel", 61.66, 59.24, "#FEBF65", Age),
    iv("Danian", "Dan", 66.0, 61.66, "#FDB462", Age),
    iv("Maastrichtian", "Maa", 72.2, 66.0, "#F2FA8C", Age),
    iv("Campanian", "Cam", 83.6, 72.2, "#E6F47F", Age),
    iv("Santonian", "San", 85.7, 83.6, "#D9EF74", Age),
    iv("Coniacian", "Con", 89.8, 85.7, "#CCE968", Age),
    iv("Turonian", "Tur", 93.9, 89.8, "#BFE35D", Age),
    iv("Cenomanian", "Cen", 100.5, 93.9, "#B3DE53", Age),
    iv("Albian", "Alb", 113.2, 100.5, "#CCEA97", Age),
    iv("Aptian", "Apt", 121.4, 113.2, "#BFE48A", Age),
    iv("Barremian", "Bar", 125.77, 121.4, "#B3DF7F", Age),
    iv("Hauterivian", "Hau", 132.6, 125.77, "#A6D975", Age),
    iv("Valanginian", "Val", 137.05, 132.6, "#99D36A", Age),
    iv("Berriasian", "Ber", 143.1, 137.05, "#8CCD60", Age),
    iv("Tithonian", "Tit", 149.2, 143.1, "#D9F1F7", Age),
    iv("Kimmeridgian", "Kim", 154.8, 149.2, "#CCECF4", Age),
    iv("Oxfordian", "Oxf", 161.5, 154.8, "#BFE7F1", Age),
    iv("Callovian", "Cal", 165.3, 161.5, "#BFE7E5", Age),
    iv("Bathonian", "Bat", 168.2, 165.3, "#B3E2E3", Age),
    iv("Bajocian", "Baj", 170.9, 168.2, "#A6DDE0", Age),
    iv("Aalenian", "Aal", 174.7, 170.9, "#9AD9DD", Age),
    iv("Toarcian", "Toa", 184.2, 174.7, "#99CEE3", Age),
    iv("Pliensbachian", "Pli", 192.9, 184.2, "#80C5DD", Age),
    iv("Sinemurian", "Sin", 199.5, 192.9, "#67BCD8", Age),
    iv("Hettangian", "Het", 201.4, 199.5, "#4EB3D3", Age),
    iv("Rhaetian", "Rha", 205.7, 201.4, "#E3B9DB", Age),
    iv("Norian", "Nor", 227.3, 205.7, "#D6AAD3", Age),
    iv("Carnian", "Car", 237.0, 227.3, "#C99BCB", Age),
    iv("Ladinian", "Lad", 241.464, 237.0, "#C983BF", Age),
    iv("Anisian", "Ani", 246.7, 241.464, "#BC75B7", Age),
    iv("Olenekian", "Ole", 249.9, 246.7, "#B051A5", Age),
    iv("Induan", "Ind", 251.902, 249.9, "#A4469F", Age),
    iv("Changhsingian", "Cha", 254.14, 251.902, "#FCC0B2", Age),
    iv("Wuchiapingian", "Wuc", 259.51, 254.14, "#FCB4A2", Age),
    iv("Capitanian", "Cap", 264.28, 259.51, "#FB9A85", Age),
    iv("Wordian", "Wor", 266.9, 264.28, "#FB8D76", Age),
    iv("Roadian", "Roa", 274.4, 266.9, "#FB8069", Age),
    iv("Kungurian", "Kun", 283.3, 274.4, "#E38776", Age),
    iv("Artinskian", "Art", 290.1, 283.3, "#E37B68", Age),
    iv("Sakmarian", "Sak", 293.52, 290.1, "#E36F5C", Age),
    iv("Asselian", "Ass", 298.9, 293.52, "#E36350", Age),
    iv("Gzhelian", "Gzh", 303.7, 298.9, "#CCD4C7", Age),
    iv("Kasimovian", "Kas", 307.0, 303.7, "#BFD0C5", Age),
    iv("Moscovian", "Mos", 315.2, 307.0, "#C7CBB9", Age),
    iv("Bashkirian", "Bas", 323.4, 315.2, "#99C2B5", Age),
    iv("Serpukhovian", "Ser", 330.3, 323.4, "#BFC26B", Age),
    iv("Visean", "Vis", 346.7, 330.3, "#A6B96C", Age),
    iv("Tournaisian", "Tou", 358.86, 346.7, "#8CB06C", Age),
    iv("Famennian", "Fam", 372.15, 358.86, "#F2EDC5", Age),
    iv("Frasnian", "Fra", 382.31, 372.15, "#F2EDAD", Age),
    iv("Givetian", "Giv", 387.95, 382.31, "#F1E185", Age),
    iv("Eifelian", "Eif", 393.47, 387.95, "#F1D576", Age),
    iv("Emsian", "Ems", 410.62, 393.47, "#E5D075", Age),
    iv("Pragian", "Pra", 413.02, 410.62, "#E5C468", Age),
    iv("Lochkovian", "Loc", 419.62, 413.02, "#E5B75A", Age),
    iv("Pridoli", "Pri", 422.7, 419.62, "#E6F5E1", Age),
    iv("Ludfordian", "Lud", 425.0, 422.7, "#D9F0DF", Age),
    iv("Gorstian", "Gor", 426.7, 425.0, "#CCECDD", Age),
    iv("Homerian", "Hom", 430.6, 426.7, "#CCEBD1", Age),
    iv("Sheinwoodian", "She", 432.9, 430.6, "#BFE6C3", Age),
    iv("Telychian", "Tel", 438.6, 432.9, "#BFE6CF", Age),
    iv("Aeronian", "Aer", 440.5, 438.6, "#B3E1C2", Age),
    iv("Rhuddanian", "Rhu", 443.1, 440.5, "#A6DCB5", Age),
    iv("Hirnantian", "Hir", 445.2, 443.1, "#A6DBAB", Age),
    iv("Katian", "Kat", 452.8, 445.2, "#99D69F", Age),
    iv("Sandbian", "San", 458.2, 452.8, "#8CD094", Age),
    iv("Darriwilian", "Dar", 469.4, 458.2, "#74C69C", Age),
    iv("Dapingian", "Dap", 471.3, 469.4, "#66C092", Age),
    iv("Floian", "Flo", 477.1, 471.3, "#41B087", Age),
    iv("Tremadocian", "Tre", 486.85, 477.1, "#33A97E", Age),
    iv("Stage 10", "S10", 491.0, 486.85, "#E6F5C9", Age),
    iv("Jiangshanian", "Jia", 494.2, 491.0, "#D9F0BB", Age),
    iv("Paibian", "Pai", 497.0, 494.2, "#CCEBAE", Age),
    iv("Guzhangian", "Guz", 500.5, 497.0, "#CCDFAA", Age),
    iv("Drumian", "Dru", 504.5, 500.5, "#BFD99D", Age),
    iv("Wuliuan", "Wul", 506.5, 504.5, "#B3D492", Age),
    iv("Stage 4", "S4", 514.5, 506.5, "#B3CA8E", Age),
    iv("Stage 3", "S3", 521.0, 514.5, "#A6C583", Age),
    iv("Stage 2", "S2", 529.0, 521.0, "#A6BA80", Age),
    iv("Fortunian", "For", 538.8, 529.0, "#99B575", Age),
];

/// Intervals of a level overlapping the age range [young, old].
pub fn overlapping(level: Level, young: f64, old: f64) -> impl Iterator<Item = &'static Interval> {
    INTERVALS.iter().filter(move |i| i.level == level && i.start > young && i.end < old)
}

/// Boundary ages (Ma) between the intervals of a level, youngest first, starting
/// with the present (0).
pub fn boundaries(level: Level) -> Vec<f64> {
    let mut b: Vec<f64> = std::iter::once(0.0).chain(INTERVALS.iter().filter(|i| i.level == level).map(|i| i.start)).collect();
    b.sort_by(|a, c| a.partial_cmp(c).unwrap());
    b.dedup();
    b
}

/// Snap a dragged node to the nearest boundary age within `tol` (branch-length
/// units). `len` is its branch length, `parent_depth` the root-to-parent
/// distance, and the present lies at depth `present_depth`, with `youngest` Ma
/// there and `ma_per_unit` Ma per branch-length unit. Returns the snapped length
/// and the boundary's age. The present wins when it is within reach (the Holocene
/// starts too close to it to tell apart at most scales).
pub fn snap_length(len: f64, parent_depth: f64, present_depth: f64, youngest: f64, ma_per_unit: f64, ages: &[f64], tol: f64) -> Option<(f64, f64)> {
    let mpu = if ma_per_unit > 0.0 { ma_per_unit } else { 1.0 };
    let near: Vec<(f64, f64)> = ages
        .iter()
        .filter(|&&a| a >= youngest - 1e-9)
        .map(|&a| (present_depth - (a - youngest) / mpu - parent_depth, a))
        .filter(|&(l, _)| l >= 0.0 && (l - len).abs() <= tol)
        .collect();
    near.iter().find(|x| x.1 == 0.0).or_else(|| near.iter().min_by(|x, y| (x.0 - len).abs().partial_cmp(&(y.0 - len).abs()).unwrap())).copied()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn snaps_to_boundaries_and_the_present() {
        let b = boundaries(Level::Epoch);
        assert_eq!(b[0], 0.0);
        assert!(b.contains(&33.9) && b.contains(&66.0));
        // Present at depth 100 (1 Ma per unit), parent at depth 50: the Eocene–Oligocene
        // boundary (33.9 Ma) is at length 16.1, the present at 50.
        assert_eq!(snap_length(16.5, 50.0, 100.0, 0.0, 1.0, &b, 1.0), Some((100.0 - 33.9 - 50.0, 33.9)));
        assert_eq!(snap_length(49.2, 50.0, 100.0, 0.0, 1.0, &b, 1.0), Some((50.0, 0.0)));
        assert_eq!(snap_length(20.0, 50.0, 100.0, 0.0, 1.0, &b, 1.0), None);
        // With 10 Ma per unit, the boundary at 33.9 Ma is 3.39 units below the present.
        let (l, a) = snap_length(46.5, 50.0, 100.0, 0.0, 10.0, &b, 0.5).unwrap();
        assert!((l - 46.61).abs() < 1e-9 && a == 33.9);
    }

    #[test]
    fn levels_are_contiguous() {
        for level in [Era, Period, Epoch, Age] {
            let mut v: Vec<&Interval> = INTERVALS.iter().filter(|i| i.level == level).collect();
            v.sort_by(|a, b| a.end.partial_cmp(&b.end).unwrap());
            for w in v.windows(2) {
                assert!((w[0].start - w[1].end).abs() < 1e-9, "gap between {} and {}", w[0].name, w[1].name);
            }
        }
        let names: Vec<&str> = overlapping(Period, 0.0, 70.0).map(|i| i.name).collect();
        assert_eq!(names, ["Quaternary", "Neogene", "Paleogene", "Cretaceous"]);
    }

    #[test]
    fn each_level_nests_in_the_next() {
        // Every boundary of a coarser level is also a boundary of each finer level.
        let pairs = [(Era, Period), (Period, Epoch), (Epoch, Age)];
        for (coarse, fine) in pairs {
            let fine_b = boundaries(fine);
            for b in boundaries(coarse) {
                assert!(fine_b.iter().any(|f| (f - b).abs() < 1e-9), "{:?} boundary {} is not a {:?} boundary", coarse, b, fine);
            }
        }
        assert_eq!(INTERVALS.iter().filter(|i| i.level == Age).count(), 102);
    }
}
