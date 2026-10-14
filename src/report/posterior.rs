//! Posteriors from what the report stores: the alignment log-likelihoods of a read and
//! the fitted proportions. `query` and `analyze` both go through here, so a `.matches`
//! table is the same whichever wrote it.

use super::{AlignRow, ReportHead};
use crate::tree::NO_NODE;
use crate::utils::EMProb;
use crate::SeqId;
use phylo::node::NodeID;
use std::collections::HashMap;

/// Length of a table indexed by [`SeqId`] over REFMAP.
fn table_len(head: &ReportHead) -> usize {
    head.refmap
        .entries
        .last()
        .map_or(0, |(seq, _)| seq.index() + 1)
}

/// `π̃_s = π_ℓ / |P_ℓ|` indexed by [`SeqId`]: the proportion of the EM leaf owning the
/// sequence's tip, shared among the sequences the leaf pools; 0 once pruned or absent.
/// The same division the EM used (`contig_pi`), so the values are bitwise its own.
pub fn pi_tilde_by_seq(head: &ReportHead) -> Vec<EMProb> {
    let owner = head.tree.em_leaf_of();
    let contigs: HashMap<NodeID, usize> = head
        .tree
        .vertices
        .iter()
        .map(|v| (v.id, v.contigs))
        .collect();
    let props: HashMap<NodeID, EMProb> = head.props.iter().copied().collect();
    let mut table = vec![0.0; table_len(head)];
    for &(seq, tip) in &head.refmap.entries {
        let leaf = owner.get(&tip).copied().unwrap_or(NO_NODE);
        table[seq.index()] = match (leaf, contigs.get(&leaf)) {
            (NO_NODE, _) | (_, None) | (_, Some(0)) => 0.0,
            (_, Some(&n)) => props.get(&leaf).copied().unwrap_or(0.0) / n as EMProb,
        };
    }
    table
}

/// The EM leaf owning every sequence, indexed by [`SeqId`] ([`NO_NODE`] once pruned or
/// absent from REFMAP).
pub fn leaf_by_seq(head: &ReportHead) -> Vec<NodeID> {
    let owner = head.tree.em_leaf_of();
    let mut table = vec![NO_NODE; table_len(head)];
    for &(seq, tip) in &head.refmap.entries {
        table[seq.index()] = owner.get(&tip).copied().unwrap_or(NO_NODE);
    }
    table
}

/// A read's posterior on each of its rows: `π̃_s · exp(ll_s − m) / Σ`, with `m` the
/// row maximum (the EM's own scaling, `scaled_row`), summed in row order. Every
/// posterior is 0 when the read has no mass (all its sequences pruned).
pub fn row_posteriors(
    pi_tilde: impl Fn(SeqId) -> EMProb,
    rows: &[AlignRow],
    out: &mut Vec<EMProb>,
) {
    out.clear();
    if rows.is_empty() {
        return;
    }
    let m = rows
        .iter()
        .map(|r| r.loglik)
        .fold(EMProb::NEG_INFINITY, EMProb::max);
    let mut d: EMProb = 0.0;
    for r in rows {
        let e = (r.loglik - m).exp();
        let w = pi_tilde(r.seq) * e;
        out.push(w);
        d += w;
    }
    if d > 0.0 {
        for w in out.iter_mut() {
            *w /= d;
        }
    } else {
        out.iter_mut().for_each(|w| *w = 0.0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row(seq: u32, loglik: EMProb) -> AlignRow {
        AlignRow {
            seq: SeqId::new(seq),
            loglik,
            pos: (0, 0),
        }
    }

    #[test]
    fn single_row_is_one_and_no_mass_is_zero() {
        let mut out = Vec::new();
        row_posteriors(|_| 0.3, &[row(5, -900.0)], &mut out);
        assert_eq!(out, vec![1.0]);
        row_posteriors(|_| 0.0, &[row(1, -1.0), row(2, -2.0)], &mut out);
        assert_eq!(out, vec![0.0, 0.0]);
        row_posteriors(|_| 0.5, &[], &mut out);
        assert!(out.is_empty());
    }

    #[test]
    fn matches_the_hand_formula_and_ignores_scale() {
        let pi = |s: SeqId| [0.5, 0.25, 0.0][s.index()];
        let mut a = Vec::new();
        row_posteriors(pi, &[row(0, -1.0), row(1, 0.0), row(2, 3.0)], &mut a);
        // pruned row 2 contributes nothing but sets the scale
        let w0 = 0.5 * (-4.0f64).exp();
        let w1 = 0.25 * (-3.0f64).exp();
        assert_eq!(a, vec![w0 / (w0 + w1), w1 / (w0 + w1), 0.0]);
        // shifting every log-likelihood by a constant leaves the posteriors unchanged
        let mut b = Vec::new();
        row_posteriors(
            pi,
            &[row(0, -1001.0), row(1, -1000.0), row(2, -997.0)],
            &mut b,
        );
        assert_eq!(a, b);
    }

    #[test]
    fn pi_tilde_pools_collapsed_and_zeroes_pruned() {
        use crate::report::read::tests::sample_head;
        let head = sample_head();
        let pi = pi_tilde_by_seq(&head);
        assert_eq!(pi.len(), 10);
        // seqs 4, 9 sit at tips 6, 7 under collapsed species 3 (π 0.4 over 2 sequences)
        assert_eq!(pi[4], 0.4 / 2.0);
        assert_eq!(pi[9], 0.4 / 2.0);
        // seq 2 at active tip 8 (π 0.6, 1 sequence); seq 7 at pruned tip 9; 0 is unknown
        assert_eq!(pi[2], 0.6);
        assert_eq!(pi[7], 0.0);
        assert_eq!(pi[0], 0.0);
        let leaf = leaf_by_seq(&head);
        assert_eq!(leaf[4], 3);
        assert_eq!(leaf[2], 8);
        assert_eq!(leaf[7], NO_NODE);
        assert_eq!(leaf[0], NO_NODE);
    }
}
