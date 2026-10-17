//! `.matches`, written from `(EM leaf, posterior, likelihood, position)` candidates: the
//! vertex each read is assigned to at `gamma`. `query` feeds it straight from the fit and
//! `analyze` from a report, so both produce the same bytes for the same `gamma`.

use crate::tree::{arena_len, children, is_synthetic_root, label, SubtreeIndex, TaxTree, NO_NODE};
use crate::utils::EMProb;
use anyhow::Result;
use phylo::node::NodeID;
use phylo::prelude::*;
use std::io::Write;

/// One sequence a read aligned to: the EM leaf owning it, the read's posterior on it,
/// the alignment log-likelihood and position.
#[derive(Clone, Copy, Debug)]
pub struct Candidate {
    pub leaf: NodeID,
    pub posterior: EMProb,
    pub loglik: EMProb,
    pub pos: (usize, usize),
}

/// Where a classified read landed: the vertex, its subtree posterior `h`, and the
/// position of the read's best-scoring sequence beneath it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Assignment {
    pub vertex: NodeID,
    pub h: EMProb,
    pub pos: (usize, usize),
}

/// Assigns reads on the final EM tree: from the root, descend to the child with the
/// largest subtree posterior while that posterior exceeds `gamma` (on a tie, the later
/// child, as [`crate::tree::assign_read`]). Without a taxonomy the tree is the root over
/// the sequences, so a read reaches a sequence only when its posterior on it exceeds
/// `gamma`.
///
/// A read costs O(hits × depth): the subtree posteriors live in dense scratch arrays
/// reset after each read, and the descent looks only at the vertices the read touched,
/// never at every child of a vertex (a species of a per-contig reference has thousands).
/// The sums are added in the order `assign_read` adds them, so the result is bitwise its.
///
/// The classifier is immutable and shared between threads; each thread classifies with
/// its own [`Scratch`].
pub struct Classifier {
    tree: TaxTree,
    subtree: SubtreeIndex,
    gamma: EMProb,
    root: NodeID,
    /// Parent of every vertex by id ([`NO_NODE`] at the root and for unused ids).
    parent: Vec<NodeID>,
    /// Position of every vertex among its parent's children.
    order: Vec<usize>,
}

/// Per-thread working memory of a [`Classifier`], all zero between reads.
pub struct Scratch {
    /// Posterior per EM leaf of the current read.
    leaf_mass: Vec<EMProb>,
    /// The current read's EM leaves with mass, in first-appearance order.
    hits: Vec<NodeID>,
    /// Subtree posterior per vertex of the current read.
    h: Vec<EMProb>,
    /// The vertices with `h > 0`, in first-touch order.
    touched: Vec<NodeID>,
}

impl Scratch {
    pub fn new(classifier: &Classifier) -> Self {
        let n = classifier.parent.len();
        Self {
            leaf_mass: vec![0.0; n],
            hits: Vec::new(),
            h: vec![0.0; n],
            touched: Vec::new(),
        }
    }

    /// Zero what the last read used.
    fn reset(&mut self) {
        for &v in &self.touched {
            self.h[v] = 0.0;
        }
        for &l in &self.hits {
            self.leaf_mass[l] = 0.0;
        }
        self.touched.clear();
        self.hits.clear();
    }
}

impl Classifier {
    /// `tree` is the final EM tree (collapsed vertices are leaves).
    pub fn new(tree: TaxTree, gamma: EMProb) -> Self {
        let subtree = SubtreeIndex::new(&tree);
        let n = arena_len(&tree);
        let mut parent = vec![NO_NODE; n];
        let mut order = vec![0; n];
        for v in tree.get_node_ids() {
            for (i, c) in children(&tree, v).into_iter().enumerate() {
                parent[c] = v;
                order[c] = i;
            }
        }
        Self {
            root: tree.get_root_id(),
            tree,
            subtree,
            gamma,
            parent,
            order,
        }
    }

    /// The deepest vertex the read loaded in `s` descends to, and its subtree posterior.
    fn descend(&self, s: &mut Scratch) -> (NodeID, EMProb) {
        for &leaf in &s.hits {
            let p = s.leaf_mass[leaf];
            let mut v = leaf;
            while v != NO_NODE {
                if s.h[v] == 0.0 {
                    s.touched.push(v);
                }
                s.h[v] += p;
                v = self.parent[v];
            }
        }
        let (mut v, mut hv) = (self.root, s.h[self.root]);
        loop {
            // The touched child with the largest h; ties to the later child.
            let best = s
                .touched
                .iter()
                .filter(|&&c| self.parent[c] == v)
                .max_by(|&&a, &&b| {
                    s.h[a]
                        .total_cmp(&s.h[b])
                        .then(self.order[a].cmp(&self.order[b]))
                });
            match best {
                Some(&c) if s.h[c] > 0.0 && s.h[c] > self.gamma => {
                    v = c;
                    hv = s.h[c];
                }
                _ => break,
            }
        }
        (v, hv)
    }

    /// Classify one read from its candidates (in [`crate::SeqId`] order): `None` when
    /// the read is unclassified (no mass, or stopped at the synthetic root).
    pub fn classify(&self, s: &mut Scratch, candidates: &[Candidate]) -> Option<Assignment> {
        for c in candidates {
            if c.posterior == 0.0 {
                continue;
            }
            if s.leaf_mass[c.leaf] == 0.0 {
                s.hits.push(c.leaf);
            }
            s.leaf_mass[c.leaf] += c.posterior;
        }
        if s.hits.is_empty() {
            return None;
        }
        let (vertex, h) = self.descend(s);
        s.reset();
        if is_synthetic_root(&self.tree, vertex) {
            return None;
        }
        let best = candidates
            .iter()
            .filter(|c| self.subtree.contains(vertex, c.leaf))
            .max_by(|a, b| EMProb::total_cmp(&a.loglik, &b.loglik))?;
        Some(Assignment {
            vertex,
            h,
            pos: best.pos,
        })
    }

    pub fn label(&self, v: NodeID) -> &str {
        label(&self.tree, v)
    }

    /// The `.matches` header line.
    /// The `.matches` head: a `#` comment line with the EM that produced the proportions
    /// ([`crate::report::ReportMeta::em_summary`]) and the gamma, then the column header.
    pub fn write_header<W: Write>(
        &self,
        out: &mut W,
        em_summary: &str,
        label_header: &str,
    ) -> Result<()> {
        writeln!(out, "# {}; gamma={}", em_summary, self.gamma)?;
        writeln!(
            out,
            "ReadID\t{}\tPosterior\tForward Position\tReverse Position",
            label_header
        )?;
        Ok(())
    }

    /// Classify one read from its candidates (in [`crate::SeqId`] order) and write its
    /// `.matches` row.
    pub fn write_row<W: Write>(
        &self,
        s: &mut Scratch,
        out: &mut W,
        read_id: &str,
        candidates: &[Candidate],
    ) -> Result<()> {
        match self.classify(s, candidates) {
            Some(a) => writeln!(
                out,
                "{}\t{}\t{:.5e}\t{}\t{}",
                read_id,
                self.label(a.vertex),
                a.h,
                a.pos.0,
                a.pos.1
            )?,
            None => writeln!(out, "{}\tunclassified\t-\t-\t-", read_id)?,
        }
        Ok(())
    }
}

/// The `.matches` table, written read by read on one thread.
pub struct MatchesWriter<W: Write> {
    classifier: Classifier,
    scratch: Scratch,
    out: W,
}

impl<W: Write> MatchesWriter<W> {
    /// Writes the head ([`Classifier::write_header`]). `label_header` is `RefID` or
    /// `TaxID`; `tree` is the final EM tree (collapsed vertices are leaves).
    pub fn new(
        mut out: W,
        em_summary: &str,
        label_header: &str,
        tree: TaxTree,
        gamma: EMProb,
    ) -> Result<Self> {
        let classifier = Classifier::new(tree, gamma);
        classifier.write_header(&mut out, em_summary, label_header)?;
        let scratch = Scratch::new(&classifier);
        Ok(Self {
            classifier,
            scratch,
            out,
        })
    }

    /// Classify one read from its candidates (in [`crate::SeqId`] order).
    pub fn push(&mut self, read_id: &str, candidates: &[Candidate]) -> Result<()> {
        self.classifier
            .write_row(&mut self.scratch, &mut self.out, read_id, candidates)
    }

    /// Flush the table and hand back its sink.
    pub fn finish(mut self) -> Result<W> {
        self.out.flush()?;
        Ok(self.out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::report::fixtures::fitted;
    use crate::report::ReportTree;
    use crate::tree::{assign_read, new_node, new_tree};

    fn c(leaf: NodeID, posterior: EMProb, loglik: EMProb, pos: (usize, usize)) -> Candidate {
        Candidate {
            leaf,
            posterior,
            loglik,
            pos,
        }
    }

    /// root → genus A (1) → sp1 (3, collapsed over two sequences), sp2 (4) → s2|x (8).
    fn final_tree() -> TaxTree {
        let (stored, tips, fit_tree, part, fit) = fitted();
        ReportTree::from_fit(&stored, &tips, &fit_tree, &part, &fit).to_tax_tree()
    }

    /// The tree of a run without a taxonomy: the root over sequences s0 (1) and s1 (2).
    fn flat_tree() -> TaxTree {
        let mut t = new_tree();
        t.add_child(0, new_node(1, "s0", "sequence"));
        t.add_child(0, new_node(2, "s1", "sequence"));
        t
    }

    fn write_on(tree: TaxTree, gamma: EMProb, reads: &[(&str, Vec<Candidate>)]) -> String {
        let mut w = MatchesWriter::new(Vec::new(), "EM test", "TaxID", tree, gamma).unwrap();
        for (id, cands) in reads {
            w.push(id, cands).unwrap();
        }
        String::from_utf8(w.finish().unwrap()).unwrap()
    }

    fn write(gamma: EMProb, reads: &[(&str, Vec<Candidate>)]) -> String {
        write_on(final_tree(), gamma, reads)
    }

    /// The data row of a single read.
    fn row(tree: TaxTree, gamma: EMProb, cands: Vec<Candidate>) -> String {
        write_on(tree, gamma, &[("r", cands)])
            .lines()
            .nth(2)
            .unwrap()
            .to_string()
    }

    #[test]
    fn header_and_unclassified_rows() {
        let out = write(
            0.6,
            &[("r0", vec![]), ("r1", vec![c(3, 0.0, 1e-3, (1, 2))])],
        );
        assert_eq!(
            out,
            "# EM test; gamma=0.6\nReadID\tTaxID\tPosterior\tForward Position\tReverse Position\nr0\tunclassified\t-\t-\t-\nr1\tunclassified\t-\t-\t-\n"
        );
    }

    #[test]
    fn descends_while_gamma_allows_and_sums_per_leaf() {
        // Two sequences pooled under collapsed species 3 (0.35 + 0.35) vs sequence 8 (0.3).
        let cands = vec![
            c(3, 0.35, 1e-9, (5, 6)),
            c(3, 0.35, 1e-8, (7, 8)),
            c(8, 0.3, 1e-7, (9, 10)),
        ];
        // gamma 0.6: genus 1 holds 1.0, species 3 holds 0.7 → stop at 3, best position under 3.
        assert_eq!(
            row(final_tree(), 0.6, cands.clone()),
            "r\tsp1\t7.00000e-1\t7\t8"
        );
        // gamma 0.8: 3 has 0.7 < 0.8 → stop at genus A with 1.0; best position over all.
        assert_eq!(
            row(final_tree(), 0.8, cands.clone()),
            "r\tA\t1.00000e0\t9\t10"
        );
        // gamma 0.7: 3 has exactly 0.7, which does not exceed it → the genus.
        assert_eq!(row(final_tree(), 0.7, cands), "r\tA\t1.00000e0\t9\t10");
    }

    #[test]
    fn without_a_taxonomy_a_sequence_needs_a_posterior_above_gamma() {
        let cands = vec![c(1, 0.8, 1e-9, (5, 6)), c(2, 0.2, 1e-7, (9, 10))];
        assert_eq!(
            row(flat_tree(), 0.6, cands.clone()),
            "r\ts0\t8.00000e-1\t5\t6"
        );
        assert_eq!(
            row(flat_tree(), 0.8, cands.clone()),
            "r\tunclassified\t-\t-\t-"
        );
        assert_eq!(row(flat_tree(), 0.9, cands), "r\tunclassified\t-\t-\t-");
        // twins split every read 0.5 / 0.5: unclassified at the default gamma
        let twins = vec![c(1, 0.5, 1e-9, (5, 6)), c(2, 0.5, 1e-9, (9, 10))];
        assert_eq!(row(flat_tree(), 0.6, twins), "r\tunclassified\t-\t-\t-");
    }

    #[test]
    fn zero_posterior_rows_are_position_candidates_only() {
        // Only sequence 8 carries mass, but a zero-weight row under 8's subtree with a
        // higher likelihood still supplies the position (as `write_matches` did).
        let cands = vec![c(8, 0.0, 1e-2, (1, 1)), c(8, 1.0, 1e-5, (2, 2))];
        assert_eq!(row(final_tree(), 0.6, cands), "r\ts2|x\t1.00000e0\t1\t1");
    }

    #[test]
    fn loglik_ties_go_to_the_later_row() {
        let cands = vec![c(8, 0.5, 1e-2, (1, 1)), c(8, 0.5, 1e-2, (2, 2))];
        assert_eq!(row(final_tree(), 0.6, cands), "r\ts2|x\t1.00000e0\t2\t2");
    }

    /// root → genus 1 → species 2, 3, each over `per_species` sequence leaves: the shape
    /// of a per-contig reference, where a species holds thousands of tips.
    fn wide_tree(per_species: usize) -> TaxTree {
        let mut t = new_tree();
        t.add_child(0, new_node(1, "g", "genus"));
        t.add_child(1, new_node(2, "s1", "species"));
        t.add_child(1, new_node(3, "s2", "species"));
        for i in 0..2 * per_species {
            let id = 4 + i;
            t.add_child(2 + i % 2, new_node(id, &format!("c{id}"), "sequence"));
        }
        t
    }

    /// Deterministic pseudo-random candidates of `n_reads` reads over the leaves 4.., with
    /// some exact posterior ties.
    fn wide_reads(n_reads: usize, per_read: usize, n_leaves: usize) -> Vec<Vec<Candidate>> {
        let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        (0..n_reads)
            .map(|_| {
                let mut leaves: Vec<NodeID> = (0..per_read)
                    .map(|_| 4 + (next() as usize) % n_leaves)
                    .collect();
                leaves.sort_unstable();
                let w: Vec<EMProb> = leaves.iter().map(|_| (next() % 4 + 1) as EMProb).collect();
                let total: EMProb = w.iter().sum();
                leaves
                    .iter()
                    .zip(&w)
                    .enumerate()
                    .map(|(i, (&leaf, &w))| c(leaf, w / total, -((next() % 3) as EMProb), (i, i)))
                    .collect()
            })
            .collect()
    }

    #[test]
    fn classifier_matches_the_reference_descent_bit_for_bit() {
        let tree = wide_tree(50);
        for gamma in [0.5, 0.6, 0.9, 0.99] {
            let cls = Classifier::new(tree.clone(), gamma);
            let mut s = Scratch::new(&cls);
            for (r, cands) in wide_reads(2000, 6, 100)
                .into_iter()
                .chain([vec![], vec![c(4, 0.0, 0.0, (0, 0))]])
                .enumerate()
            {
                let mut hits: Vec<(NodeID, EMProb)> = Vec::new();
                for x in cands.iter().filter(|x| x.posterior != 0.0) {
                    match hits.iter_mut().find(|(l, _)| *l == x.leaf) {
                        Some((_, p)) => *p += x.posterior,
                        None => hits.push((x.leaf, x.posterior)),
                    }
                }
                let want = (!hits.is_empty())
                    .then(|| assign_read(&tree, &hits, gamma))
                    .filter(|&(v, _)| !is_synthetic_root(&tree, v));
                let got = cls.classify(&mut s, &cands);
                assert_eq!(
                    got.map(|a| (a.vertex, a.h.to_bits())),
                    want.map(|(v, h)| (v, h.to_bits())),
                    "gamma {gamma} read {r}: {cands:?}"
                );
            }
        }
    }

    /// `cargo test --release -- --ignored --nocapture classify_wide`: per-read cost on a
    /// per-contig-shaped tree.
    #[test]
    #[ignore]
    fn classify_wide_tree_timing() {
        let per_species = 5000;
        for per_read in [200, 1000] {
            let reads = wide_reads(20_000, per_read, 2 * per_species);
            let cls = Classifier::new(wide_tree(per_species), 0.6);
            let mut s = Scratch::new(&cls);
            let t = std::time::Instant::now();
            let n = reads.iter().filter_map(|r| cls.classify(&mut s, r)).count();
            let s = t.elapsed().as_secs_f64();
            eprintln!(
                "{} reads x {per_read} hits: {:.3} s ({:.1} us/read), {n} classified",
                reads.len(),
                s,
                s * 1e6 / reads.len() as f64
            );
            // The previous classifier: per-leaf sums by linear search, then assign_read.
            let tree = wide_tree(per_species);
            let t = std::time::Instant::now();
            let mut hits: Vec<(NodeID, EMProb)> = Vec::new();
            for r in &reads {
                hits.clear();
                for x in r.iter().filter(|x| x.posterior != 0.0) {
                    match hits.iter_mut().find(|(l, _)| *l == x.leaf) {
                        Some((_, p)) => *p += x.posterior,
                        None => hits.push((x.leaf, x.posterior)),
                    }
                }
                std::hint::black_box(assign_read(&tree, &hits, 0.6));
            }
            let s = t.elapsed().as_secs_f64();
            eprintln!(
                "  previous: {:.3} s ({:.1} us/read)",
                s,
                s * 1e6 / reads.len() as f64
            );
        }
    }

    #[test]
    fn split_mass_stops_at_the_genus() {
        let cands = vec![c(3, 0.5, 1e-2, (1, 1)), c(8, 0.5, 1e-2, (2, 2))];
        // Both under genus A (1.0), neither species above 0.6 → A.
        assert!(row(final_tree(), 0.6, cands).starts_with("r\tA\t1.00000e0"));
    }
}
