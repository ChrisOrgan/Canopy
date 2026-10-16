//! Click-to-sort column headers for the app's tables. A table's sort (column,
//! ascending) lives in egui memory under the table's id; clicking a header
//! sorts ascending, then descending, then back to the original order.

use egui::{Id, RichText, Ui};
use std::cmp::Ordering;

type SortState = (usize, bool);

fn state(ui: &Ui, id: Id) -> Option<SortState> {
    ui.data(|d| d.get_temp::<SortState>(id))
}

/// A column header that sorts its table when clicked, `width` wide if given.
pub fn header_sized(ui: &mut Ui, id: Id, col: usize, text: &str, width: Option<f32>) {
    let st = state(ui, id);
    let arrow = match st {
        Some((c, true)) if c == col => " ▲",
        Some((c, false)) if c == col => " ▼",
        _ => "",
    };
    // Full width so the column fits its header; only fixed-width columns (the
    // data editor) shorten it, with the whole name on hover.
    let label = egui::Label::new(RichText::new(format!("{}{}", text, arrow)).strong()).sense(egui::Sense::click()).extend();
    let resp = match width {
        Some(w) => ui.add_sized([w, 18.0], label.truncate()).on_hover_text(text),
        None => ui.add(label),
    };
    if resp.on_hover_cursor(egui::CursorIcon::PointingHand).on_hover_text("Click to sort").clicked() {
        let next = match st {
            Some((c, true)) if c == col => Some((col, false)),
            Some((c, false)) if c == col => None,
            _ => Some((col, true)),
        };
        ui.data_mut(|d| match next {
            Some(s) => d.insert_temp(id, s),
            None => d.remove::<SortState>(id),
        });
    }
}

/// A header row of sortable columns.
pub fn headers<S: AsRef<str>>(ui: &mut Ui, id: Id, names: &[S]) {
    for (i, h) in names.iter().enumerate() {
        header_sized(ui, id, i, h.as_ref(), None);
    }
}

/// Row indices in display order for the table's current sort.
pub fn order(ui: &Ui, id: Id, rows: &[Vec<String>]) -> Vec<usize> {
    sorted(rows, state(ui, id))
}

/// Rows sorted by a column (numbers numerically, then text without regard to
/// case; empty cells last either way), else in their original order.
pub fn sorted(rows: &[Vec<String>], sort: Option<SortState>) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..rows.len()).collect();
    if let Some((c, asc)) = sort {
        let cell = |i: usize| rows[i].get(c).map(|s| s.trim()).unwrap_or("");
        idx.sort_by(|&a, &b| match (cell(a).is_empty(), cell(b).is_empty()) {
            (true, true) => Ordering::Equal,
            (true, false) => Ordering::Greater,
            (false, true) => Ordering::Less,
            _ => {
                let o = compare(cell(a), cell(b));
                if asc { o } else { o.reverse() }
            }
        });
    }
    idx
}

/// The number a cell starts with: "0.5", "1e-3", "[0.1, 0.9]" (interval), "12 (3%)".
fn number(s: &str) -> Option<f64> {
    let s = s.trim_start_matches(['[', '{', '(']);
    let end = s.find([',', ' ', ']', '}', ')']).unwrap_or(s.len());
    s[..end].parse::<f64>().ok().filter(|x| !x.is_nan())
}

/// Numbers before text; numbers by value, text alphabetically (ignoring case).
pub fn compare(a: &str, b: &str) -> Ordering {
    match (number(a), number(b)) {
        (Some(x), Some(y)) => x.partial_cmp(&y).unwrap_or(Ordering::Equal),
        (Some(_), None) => Ordering::Less,
        (None, Some(_)) => Ordering::Greater,
        _ => a.to_lowercase().cmp(&b.to_lowercase()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sorts_numbers_text_and_blanks() {
        let rows: Vec<Vec<String>> = [["b", "10"], ["A", "9"], ["c", ""], ["a2", "[0.5, 1]"], ["d", "-1e-3"]].iter().map(|r| r.iter().map(|s| s.to_string()).collect()).collect();
        assert_eq!(sorted(&rows, None), [0, 1, 2, 3, 4]);
        assert_eq!(sorted(&rows, Some((1, true))), [4, 3, 1, 0, 2]);
        // Descending keeps blanks last.
        assert_eq!(sorted(&rows, Some((1, false))), [0, 1, 3, 4, 2]);
        assert_eq!(sorted(&rows, Some((0, true))), [1, 3, 0, 2, 4]);
    }
}
