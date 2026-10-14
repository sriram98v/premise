//! The legacy per-read tables, written from `(EM leaf, posterior, likelihood, position)`
//! candidates: `.matches` (the vertex each read is assigned to at `gamma`) and
//! `.posteriors` (the read's posterior on every EM leaf). `query` feeds them straight
//! from the fit and `analyze` from a report, so both produce the same bytes for the
//! same `gamma`.

use crate::tree::{assign_read, is_synthetic_root, label, SubtreeIndex, TaxTree};
use crate::utils::EMProb;
use anyhow::Result;
use phylo::node::NodeID;
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

/// Assigns reads on the final EM tree at a `gamma`; one classification feeds every table.
pub struct Classifier {
    tree: TaxTree,
    subtree: SubtreeIndex,
    gamma: EMProb,
    /// Posterior per EM leaf of the last read, in first-appearance order.
    hits: Vec<(NodeID, EMProb)>,
}

impl Classifier {
    /// `tree` is the final EM tree (collapsed vertices are leaves).
    pub fn new(tree: TaxTree, gamma: EMProb) -> Self {
        let subtree = SubtreeIndex::new(&tree);
        Self {
            tree,
            subtree,
            gamma,
            hits: Vec::new(),
        }
    }

    /// Classify one read from its candidates (in [`crate::SeqId`] order): `None` when
    /// the read is unclassified (no mass, or stopped at the synthetic root). Afterwards
    /// [`Classifier::hits`] holds its posterior per EM leaf.
    pub fn classify(&mut self, candidates: &[Candidate]) -> Option<Assignment> {
        self.hits.clear();
        for c in candidates {
            if c.posterior == 0.0 {
                continue;
            }
            match self.hits.iter_mut().find(|(l, _)| *l == c.leaf) {
                Some((_, p)) => *p += c.posterior,
                None => self.hits.push((c.leaf, c.posterior)),
            }
        }
        if self.hits.is_empty() {
            return None;
        }
        let (vertex, h) = assign_read(&self.tree, &self.hits, self.gamma);
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

    /// Posterior per EM leaf of the last classified read (pooled sequences summed).
    pub fn hits(&self) -> &[(NodeID, EMProb)] {
        &self.hits
    }

    pub fn label(&self, v: NodeID) -> &str {
        label(&self.tree, v)
    }
}

/// `.matches` and / or `.posteriors`, written read by read from one classification.
pub struct TableWriters<W: Write> {
    classifier: Classifier,
    matches: Option<W>,
    posteriors: Option<W>,
}

impl<W: Write> TableWriters<W> {
    /// Writes the headers. `label_header` is `RefID` or `TaxID`.
    pub fn new(
        label_header: &str,
        tree: TaxTree,
        gamma: EMProb,
        mut matches: Option<W>,
        mut posteriors: Option<W>,
    ) -> Result<Self> {
        if let Some(out) = matches.as_mut() {
            writeln!(
                out,
                "ReadID\t{}\tPosterior\tForward Position\tReverse Position",
                label_header
            )?;
        }
        if let Some(out) = posteriors.as_mut() {
            writeln!(out, "ReadID\t{}\tPosterior", label_header)?;
        }
        Ok(Self {
            classifier: Classifier::new(tree, gamma),
            matches,
            posteriors,
        })
    }

    /// One read's rows in every table.
    pub fn push(&mut self, read_id: &str, candidates: &[Candidate]) -> Result<()> {
        let assigned = self.classifier.classify(candidates);
        if let Some(out) = self.matches.as_mut() {
            match assigned {
                Some(a) => writeln!(
                    out,
                    "{}\t{}\t{:.5e}\t{}\t{}",
                    read_id,
                    self.classifier.label(a.vertex),
                    a.h,
                    a.pos.0,
                    a.pos.1
                )?,
                None => writeln!(out, "{}\tunclassified\t-\t-\t-", read_id)?,
            }
        }
        if let Some(out) = self.posteriors.as_mut() {
            if assigned.is_none() {
                writeln!(out, "{}\tunclassified\t-", read_id)?;
            } else {
                let mut rows: Vec<(&str, EMProb)> = self
                    .classifier
                    .hits()
                    .iter()
                    .map(|&(leaf, p)| (self.classifier.label(leaf), p))
                    .collect();
                rows.sort_by(|a, b| a.0.cmp(b.0));
                for (name, p) in rows {
                    writeln!(out, "{}\t{}\t{:.5e}", read_id, name, p)?;
                }
            }
        }
        Ok(())
    }

    /// Flush both tables and hand back their sinks `(matches, posteriors)`.
    pub fn finish(mut self) -> Result<(Option<W>, Option<W>)> {
        if let Some(out) = self.matches.as_mut() {
            out.flush()?;
        }
        if let Some(out) = self.posteriors.as_mut() {
            out.flush()?;
        }
        Ok((self.matches, self.posteriors))
    }
}

/// The `.matches` table alone.
pub struct MatchesWriter<W: Write>(TableWriters<W>);

impl<W: Write> MatchesWriter<W> {
    /// Writes the header. `tree` is the final EM tree (collapsed vertices are leaves).
    pub fn new(out: W, label_header: &str, tree: TaxTree, gamma: EMProb) -> Result<Self> {
        Ok(Self(TableWriters::new(
            label_header,
            tree,
            gamma,
            Some(out),
            None,
        )?))
    }

    /// Classify one read from its candidates (in [`crate::SeqId`] order).
    pub fn push(&mut self, read_id: &str, candidates: &[Candidate]) -> Result<()> {
        self.0.push(read_id, candidates)
    }

    pub fn finish(self) -> Result<W> {
        let (matches, _) = self.0.finish()?;
        Ok(matches.expect("the matches sink was given"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::report::fixtures::fitted;
    use crate::report::ReportTree;

    fn c(leaf: NodeID, posterior: EMProb, loglik: EMProb, pos: (usize, usize)) -> Candidate {
        Candidate {
            leaf,
            posterior,
            loglik,
            pos,
        }
    }

    fn final_tree() -> TaxTree {
        let (stored, tips, fit_tree, part, fit) = fitted();
        ReportTree::from_fit(&stored, &tips, &fit_tree, &part, &fit).to_tax_tree()
    }

    /// Both tables of a few reads: `(matches, posteriors)`.
    fn write_both(gamma: EMProb, reads: &[(&str, Vec<Candidate>)]) -> (String, String) {
        let mut w = TableWriters::new(
            "TaxID",
            final_tree(),
            gamma,
            Some(Vec::new()),
            Some(Vec::new()),
        )
        .unwrap();
        for (id, cands) in reads {
            w.push(id, cands).unwrap();
        }
        let (m, p) = w.finish().unwrap();
        (
            String::from_utf8(m.unwrap()).unwrap(),
            String::from_utf8(p.unwrap()).unwrap(),
        )
    }

    fn write(gamma: EMProb, reads: &[(&str, Vec<Candidate>)]) -> String {
        let mut w = MatchesWriter::new(Vec::new(), "TaxID", final_tree(), gamma).unwrap();
        for (id, cands) in reads {
            w.push(id, cands).unwrap();
        }
        String::from_utf8(w.finish().unwrap()).unwrap()
    }

    #[test]
    fn header_and_unclassified_rows() {
        let out = write(
            0.6,
            &[("r0", vec![]), ("r1", vec![c(3, 0.0, 1e-3, (1, 2))])],
        );
        assert_eq!(
            out,
            "ReadID\tTaxID\tPosterior\tForward Position\tReverse Position\nr0\tunclassified\t-\t-\t-\nr1\tunclassified\t-\t-\t-\n"
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
        let out = write(0.6, &[("r", cands.clone())]);
        assert_eq!(out.lines().nth(1).unwrap(), "r\tsp1\t7.00000e-1\t7\t8");
        // gamma 0.8: 3 has 0.7 < 0.8 → stop at genus A with 1.0; best position over all.
        let out = write(0.8, &[("r", cands.clone())]);
        assert_eq!(out.lines().nth(1).unwrap(), "r\tA\t1.00000e0\t9\t10");
        // gamma 0.7: 3 has exactly 0.7, which does not exceed it → the genus.
        let out = write(0.7, &[("r", cands)]);
        assert_eq!(out.lines().nth(1).unwrap(), "r\tA\t1.00000e0\t9\t10");
    }

    #[test]
    fn zero_posterior_rows_are_position_candidates_only() {
        // Only sequence 8 carries mass, but a zero-weight row under 8's subtree with a
        // higher likelihood still supplies the position (as `write_matches` did).
        let cands = vec![c(8, 0.0, 1e-2, (1, 1)), c(8, 1.0, 1e-5, (2, 2))];
        let out = write(0.6, &[("r", cands)]);
        assert_eq!(out.lines().nth(1).unwrap(), "r\ts2|x\t1.00000e0\t1\t1");
    }

    #[test]
    fn loglik_ties_go_to_the_later_row() {
        let cands = vec![c(8, 0.5, 1e-2, (1, 1)), c(8, 0.5, 1e-2, (2, 2))];
        let out = write(0.6, &[("r", cands)]);
        assert_eq!(out.lines().nth(1).unwrap(), "r\ts2|x\t1.00000e0\t2\t2");
    }

    #[test]
    fn split_mass_stops_at_the_root() {
        let cands = vec![c(3, 0.5, 1e-2, (1, 1)), c(8, 0.5, 1e-2, (2, 2))];
        // Both under genus A (1.0) → A; make the root split instead by giving mass to a
        // vertex outside A: there is none in the fitted tree, so A always wins here.
        let out = write(0.6, &[("r", cands)]);
        assert!(out.lines().nth(1).unwrap().starts_with("r\tA\t1.00000e0"));
    }

    #[test]
    fn posteriors_list_every_leaf_by_label_and_follow_the_matches_verdict() {
        let classified = vec![
            c(8, 0.3, 1e-7, (9, 10)),
            c(3, 0.35, 1e-9, (5, 6)),
            c(3, 0.35, 1e-8, (7, 8)),
        ];
        let (m, p) = write_both(
            0.6,
            &[
                ("r0", classified),
                ("r1", vec![]),
                ("r2", vec![c(8, 0.0, 1e-2, (1, 1))]),
            ],
        );
        assert_eq!(
            m,
            "ReadID\tTaxID\tPosterior\tForward Position\tReverse Position\n\
             r0\tsp1\t7.00000e-1\t7\t8\n\
             r1\tunclassified\t-\t-\t-\n\
             r2\tunclassified\t-\t-\t-\n"
        );
        // Pooled sequences 6, 7 are summed into species 3; rows sorted by label.
        assert_eq!(
            p,
            "ReadID\tTaxID\tPosterior\n\
             r0\ts2|x\t3.00000e-1\n\
             r0\tsp1\t7.00000e-1\n\
             r1\tunclassified\t-\n\
             r2\tunclassified\t-\n"
        );
    }

    #[test]
    fn posteriors_alone_need_no_matches_sink() {
        let mut w = TableWriters::new("RefID", final_tree(), 0.9, None, Some(Vec::new())).unwrap();
        w.push("r", &[c(8, 1.0, 1e-2, (1, 1))]).unwrap();
        let (m, p) = w.finish().unwrap();
        assert!(m.is_none());
        assert_eq!(
            String::from_utf8(p.unwrap()).unwrap(),
            "ReadID\tRefID\tPosterior\nr\ts2|x\t1.00000e0\n"
        );
    }
}
