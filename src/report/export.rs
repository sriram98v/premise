//! Tables regenerated from a report head: the NEXUS `.props` tree, and the flat leaf
//! proportions the web UI charts.

use super::nexus::quote;
use super::write::newick;
use super::{ReportHead, ReportTree};
use std::collections::HashMap;
use std::fmt::Write as _;

/// Version of the `.props` layout, in its `[PREMISE proportions vN]` header comment.
pub const PROPS_VERSION: u32 = 1;

/// The final-tree vertices with `π_v > 0`, in pre-order. `π` never grows down the tree,
/// so the parent of every kept vertex is kept and the result is a tree under the root.
fn weighted_tree(head: &ReportHead) -> ReportTree {
    ReportTree {
        vertices: head
            .tree
            .vertices
            .iter()
            .filter(|v| v.status.is_final() && v.proportion > 0.0)
            .cloned()
            .collect(),
    }
}

/// `.props`: the taxonomy as NEXUS (TAXA, TREES, PROPORTIONS) with the end-of-EM `π_v`,
/// the mass of the subtree under `v`, on every final-tree vertex that has any.
pub fn props_nexus(head: &ReportHead) -> String {
    let tree = weighted_tree(head);
    let label_of: HashMap<usize, &str> = tree
        .vertices
        .iter()
        .map(|v| (v.id, v.label.as_str()))
        .collect();
    let mut out = String::new();
    // Writing into a String cannot fail.
    let _ = writeln!(
        out,
        "#NEXUS\n[PREMISE proportions v{}]\n[{}]\n",
        PROPS_VERSION,
        head.meta.em_summary()
    );
    let _ = writeln!(out, "BEGIN TAXA;");
    let _ = writeln!(out, "\tDIMENSIONS NTAX={};", tree.vertices.len());
    let _ = writeln!(out, "\tTAXLABELS");
    for v in &tree.vertices {
        let _ = writeln!(out, "\t\t{}", quote(&v.label));
    }
    let _ = writeln!(out, "\t;\nEND;\n");
    let _ = writeln!(out, "BEGIN TREES;");
    let _ = writeln!(out, "\tTREE proportions = [&R] {};", newick(&tree));
    let _ = writeln!(out, "END;\n");
    let _ = writeln!(out, "BEGIN PROPORTIONS;");
    let _ = writeln!(out, "\tDIMENSIONS N={};", tree.vertices.len());
    let _ = writeln!(
        out,
        "\t[label parent rank proportion; row i is TAXLABELS i, proportion the mass of the subtree under it]"
    );
    let _ = writeln!(out, "\tMATRIX");
    for v in &tree.vertices {
        let parent = v.parent.and_then(|p| label_of.get(&p).copied());
        let _ = writeln!(
            out,
            "\t{} {} {} {:e}",
            quote(&v.label),
            parent.map_or_else(|| "-".into(), quote),
            quote(&v.rank),
            v.proportion
        );
    }
    let _ = writeln!(out, "\t;\nEND;");
    out
}

/// `label\tπ_ℓ` per final-tree leaf with mass, sorted by label, no header: what the web
/// UI charts.
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
        let _ = writeln!(out, "{}\t{:.10e}", label, p);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::report::nexus::tokens;
    use crate::report::read::tests::sample_head;

    #[test]
    fn props_tsv_is_sorted_by_label_without_header() {
        let head = sample_head();
        assert_eq!(
            props_tsv(&head),
            "s2|x\t6.0000000000e-1\nsp1\t4.0000000000e-1\n"
        );
    }

    /// The token rows of a block's MATRIX.
    fn matrix(text: &str, block: &str) -> Vec<Vec<String>> {
        let start = text.find(&format!("BEGIN {block};")).unwrap();
        let body = &text[start..];
        let m = body.find("MATRIX").unwrap() + "MATRIX".len();
        let end = body[m..].find("\t;").unwrap() + m;
        body[m..end]
            .lines()
            .map(|l| tokens(l).unwrap())
            .filter(|t| !t.is_empty())
            .collect()
    }

    #[test]
    fn props_nexus_holds_every_final_vertex_with_mass() {
        let head = sample_head();
        let out = props_nexus(&head);
        assert!(
            out.starts_with(
                "#NEXUS\n[PREMISE proportions v1]\n[PREMISE 0.6.0 EM: init=split iter=100 rho=1.5e2 omega=1e-10 em_threshold=1e-6 em_stop=loglik penalty=1; 1 iterations, final penalized log-likelihood -8e0]\n"
            ),
            "{out}"
        );
        // A refit names the report it started from; brackets cannot end the comment.
        let mut refit = sample_head();
        refit.meta.source = Some("base [1].report".into());
        let summary = refit.meta.em_summary();
        assert!(
            summary.ends_with("-8e0, refit of base _1_.report"),
            "{summary}"
        );
        let rows = matrix(&out, "PROPORTIONS");
        let flat: Vec<Vec<&str>> = rows
            .iter()
            .map(|r| r.iter().map(String::as_str).collect())
            .collect();
        assert_eq!(
            flat,
            vec![
                vec!["root", "-", "root", "1e0"],
                vec!["A", "root", "genus", "1e0"],
                vec!["sp1", "A", "species", "4e-1"],
                vec!["sp2", "A", "species", "6e-1"],
                vec!["s2|x", "sp2", "sequence", "6e-1"],
            ]
        );
        // pruned genus B and the pooled sequences under sp1 are absent
        assert!(!out.contains("\tB ") && !out.contains("s0"), "{out}");
        assert!(out.contains("\tDIMENSIONS NTAX=5;") && out.contains("\tDIMENSIONS N=5;"));
        // TAXA labels in row order, quoted where NEXUS needs it
        let taxa = &out[out.find("TAXLABELS").unwrap()..out.find("BEGIN TREES").unwrap()];
        let labels: Vec<String> = taxa
            .lines()
            .skip(1)
            .flat_map(|l| tokens(l).unwrap())
            .collect();
        assert_eq!(labels, ["root", "A", "sp1", "sp2", "s2|x", ";", "END", ";"]);
        assert!(taxa.contains("'s2|x'"));
        assert!(
            out.contains("\tTREE proportions = [&R] ((3,(5)4)2)1;"),
            "{out}"
        );
    }
}
