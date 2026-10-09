//! The geologic time scale (International Chronostratigraphic Chart,
//! v2023/09) with the standard CGMW colors, for timescale bars under
//! time-calibrated trees (like deeptime's `coord_geo`). Ages in Ma.

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Level {
    Era,
    Period,
    Epoch,
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
    iv("Neogene", "N", 23.03, 2.58, "#FFE619", Period),
    iv("Paleogene", "Pg", 66.0, 23.03, "#FD9A52", Period),
    iv("Cretaceous", "K", 145.0, 66.0, "#7FC64E", Period),
    iv("Jurassic", "J", 201.4, 145.0, "#34B2C9", Period),
    iv("Triassic", "Tr", 251.902, 201.4, "#812B92", Period),
    iv("Permian", "P", 298.9, 251.902, "#F04028", Period),
    iv("Carboniferous", "C", 358.9, 298.9, "#67A599", Period),
    iv("Devonian", "D", 419.2, 358.9, "#CB8C37", Period),
    iv("Silurian", "S", 443.8, 419.2, "#B3E1B6", Period),
    iv("Ordovician", "O", 485.4, 443.8, "#009270", Period),
    iv("Cambrian", "Cm", 538.8, 485.4, "#7FA056", Period),
    // Epochs
    iv("Holocene", "Hol", 0.0117, 0.0, "#FEF2E0", Epoch),
    iv("Pleistocene", "Ple", 2.58, 0.0117, "#FFF2AE", Epoch),
    iv("Pliocene", "Pli", 5.333, 2.58, "#FFFF99", Epoch),
    iv("Miocene", "Mio", 23.03, 5.333, "#FFFF00", Epoch),
    iv("Oligocene", "Oli", 33.9, 23.03, "#FDC07A", Epoch),
    iv("Eocene", "Eoc", 56.0, 33.9, "#FDB46C", Epoch),
    iv("Paleocene", "Pal", 66.0, 56.0, "#FDA75F", Epoch),
    iv("Late Cretaceous", "LK", 100.5, 66.0, "#A6D84A", Epoch),
    iv("Early Cretaceous", "EK", 145.0, 100.5, "#8CCD57", Epoch),
    iv("Late Jurassic", "LJ", 161.5, 145.0, "#B3E3EE", Epoch),
    iv("Middle Jurassic", "MJ", 174.7, 161.5, "#80CFD8", Epoch),
    iv("Early Jurassic", "EJ", 201.4, 174.7, "#42AED0", Epoch),
    iv("Late Triassic", "LTr", 237.0, 201.4, "#BD8CC3", Epoch),
    iv("Middle Triassic", "MTr", 247.2, 237.0, "#B168B1", Epoch),
    iv("Early Triassic", "ETr", 251.902, 247.2, "#A4469F", Epoch),
    iv("Lopingian", "Lop", 259.51, 251.902, "#FBA794", Epoch),
    iv("Guadalupian", "Gua", 274.4, 259.51, "#FB745C", Epoch),
    iv("Cisuralian", "Cis", 298.9, 274.4, "#EF5845", Epoch),
    iv("Pennsylvanian", "Penn", 323.2, 298.9, "#99C2B5", Epoch),
    iv("Mississippian", "Miss", 358.9, 323.2, "#678F66", Epoch),
    iv("Late Devonian", "LD", 382.7, 358.9, "#F1E19D", Epoch),
    iv("Middle Devonian", "MD", 393.3, 382.7, "#F1C868", Epoch),
    iv("Early Devonian", "ED", 419.2, 393.3, "#E5AC4D", Epoch),
    iv("Pridoli", "Pri", 423.0, 419.2, "#E6F5E1", Epoch),
    iv("Ludlow", "Lud", 427.4, 423.0, "#BFE6CF", Epoch),
    iv("Wenlock", "Wen", 433.4, 427.4, "#B3E1C2", Epoch),
    iv("Llandovery", "Lla", 443.8, 433.4, "#99D7B3", Epoch),
    iv("Late Ordovician", "LO", 458.4, 443.8, "#7FCA93", Epoch),
    iv("Middle Ordovician", "MO", 470.0, 458.4, "#4DB47E", Epoch),
    iv("Early Ordovician", "EO", 485.4, 470.0, "#1A9D6F", Epoch),
    iv("Furongian", "Fur", 497.0, 485.4, "#B3E095", Epoch),
    iv("Miaolingian", "Mia", 509.0, 497.0, "#A6CF86", Epoch),
    iv("Series 2", "S2", 521.0, 509.0, "#99C078", Epoch),
    iv("Terreneuvian", "Ter", 538.8, 521.0, "#8CB06C", Epoch),
];

/// Intervals of a level overlapping the age range [young, old].
pub fn overlapping(level: Level, young: f64, old: f64) -> impl Iterator<Item = &'static Interval> {
    INTERVALS.iter().filter(move |i| i.level == level && i.start > young && i.end < old)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn levels_are_contiguous() {
        for level in [Era, Period, Epoch] {
            let mut v: Vec<&Interval> = INTERVALS.iter().filter(|i| i.level == level).collect();
            v.sort_by(|a, b| a.end.partial_cmp(&b.end).unwrap());
            for w in v.windows(2) {
                assert!((w[0].start - w[1].end).abs() < 1e-9, "gap between {} and {}", w[0].name, w[1].name);
            }
        }
        let names: Vec<&str> = overlapping(Period, 0.0, 70.0).map(|i| i.name).collect();
        assert_eq!(names, ["Quaternary", "Neogene", "Paleogene", "Cretaceous"]);
    }
}
