//! The legacy `.props` and `.nodes` tables, regenerated from a report head.

use super::{ReportHead, VertexStatus};
use std::collections::HashMap;
use std::fmt::Write as _;

/// `.props`: `label\tπ_ℓ` per final-tree leaf with mass, sorted by label, no header.
pub fn props_tsv(head: &ReportHead) -> String {
    let labels = head.tree.labels();
    let mut rows: Vec<(&str, f64)> = head
        .props
        .iter()
        .map(|(node, p)| (labels.get(node).copied().unwrap_or(""), *p))
        .collect();
    rows.sort_by(|a, b| a.0.cmp(b.0));
    let mut out = String::new();
    for (label, p) in rows {
        // Writing into a String cannot fail.
        let _ = writeln!(out, "{}\t{:.10e}", label, p);
    }
    out
}

/// `.nodes`: every final-tree vertex in pre-order with the columns `query` used to
/// write (`TaxID Parent Rank Contigs Count Proportion Conditional`).
pub fn nodes_tsv(head: &ReportHead) -> String {
    let label_of: HashMap<usize, &str> = head
        .tree
        .vertices
        .iter()
        .map(|v| (v.id, v.label.as_str()))
        .collect();
    let mut out = String::from("TaxID\tParent\tRank\tContigs\tCount\tProportion\tConditional\n");
    for v in &head.tree.vertices {
        if !matches!(v.status, VertexStatus::Active | VertexStatus::Collapsed) {
            continue;
        }
        let parent = v
            .parent
            .and_then(|p| label_of.get(&p).copied())
            .unwrap_or("-");
        let _ = writeln!(
            out,
            "{}\t{}\t{}\t{}\t{:.3}\t{:.10e}\t{:.10e}",
            v.label,
            parent,
            v.rank,
            v.contigs,
            v.count.unwrap_or(0.0),
            v.proportion,
            v.conditional.unwrap_or(0.0),
        );
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::report::read::tests::sample_head;

    #[test]
    fn props_are_sorted_by_label_without_header() {
        let head = sample_head();
        assert_eq!(
            props_tsv(&head),
            "s2|x\t6.0000000000e-1\nsp1\t4.0000000000e-1\n"
        );
    }

    #[test]
    fn nodes_lists_only_final_vertices_in_preorder() {
        let head = sample_head();
        let out = nodes_tsv(&head);
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(
            lines[0],
            "TaxID\tParent\tRank\tContigs\tCount\tProportion\tConditional"
        );
        assert_eq!(
            lines[1],
            "root\t-\troot\t0\t55.000\t1.0000000000e0\t1.0000000000e0"
        );
        assert_eq!(
            lines[2].split('\t').take(4).collect::<Vec<_>>(),
            ["A", "root", "genus", "0"]
        );
        assert_eq!(
            lines[3].split('\t').take(4).collect::<Vec<_>>(),
            ["sp1", "A", "species", "2"]
        );
        assert_eq!(
            lines[4].split('\t').take(2).collect::<Vec<_>>(),
            ["sp2", "A"]
        );
        assert_eq!(
            lines[5].split('\t').take(2).collect::<Vec<_>>(),
            ["s2|x", "sp2"]
        );
        assert_eq!(lines.len(), 6);
    }
}
