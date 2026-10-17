//! Abundance estimation: the penalized EM over the taxonomy
//!
//! [`CsrLikelihood`] is the sparse read × sequence matrix of alignment likelihoods
//! `A = (a_is)`. The EM parameters live on a [`TaxTree`] whose leaves own the
//! sequences ([`Partition`]).

#[cfg(test)]
use crate::tree::assign_read;
#[cfg(test)]
use crate::tree::NO_NODE;
use crate::tree::{
    contig_pi, init_conditionals, init_uniform, leaf_marginals, m_step, penalty,
    prune_and_collapse, subtree_counts, Partition, TaxTree,
};
use crate::utils::{EMProb, QueryProgress, ReadAlignments, ReadIdx};
use crate::SeqId;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use phylo::node::NodeID;
use phylo::prelude::*;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::sync::atomic::Ordering;
use std::sync::Arc;

/// The read × sequence likelihood matrix, stored by row (CSR) and by column (CSC).
pub struct CsrLikelihood {
    /// Compact column → original [`SeqId`] of the sequence.
    refs: Vec<SeqId>,
    /// Original [`SeqId`] → compact column.
    ref_compact: HashMap<SeqId, u32>,
    /// CSR row → original [`ReadIdx`].
    reads: Vec<ReadIdx>,
    /// `row_of[read_idx]` = the row holding that read, or [`NO_ROW`] (read by the
    /// test-only [`Posteriors`]).
    #[cfg_attr(not(test), allow(dead_code))]
    row_of: Vec<u32>,
    /// `row_ptr[r]..row_ptr[r + 1]` bounds row `r`'s entries in `col`/`lik`.
    row_ptr: Vec<usize>,
    /// Compact column for each stored entry; sorted within a row.
    col: Vec<u32>,
    /// `a_is`, scaled so that the largest entry of every row is 1.
    lik: Vec<EMProb>,
    /// Log of the scale divided out of each row (its largest log-likelihood).
    log_scale: Vec<EMProb>,
    /// `col_ptr[s]..col_ptr[s + 1]` bounds column `s`'s entries in `row_idx`/`lik_t`.
    col_ptr: Vec<usize>,
    /// Row of each column-ordered entry; ascending within a column.
    row_idx: Vec<u32>,
    /// `a_is` in column order.
    lik_t: Vec<EMProb>,
}

/// Marker in [`CsrLikelihood::row_of`] for reads without a row.
const NO_ROW: u32 = u32::MAX;

/// Largest read → row table [`CsrBuilder::new`] allocates up front; a file-announced
/// read count beyond it grows the table as rows arrive.
const MAX_ROW_PREALLOC: usize = 1 << 16;

/// Scale one read's log-likelihoods by their maximum: `(exp(ll − max), max)`. Every
/// responsibility is unchanged and nothing underflows to 0.
fn scaled_row(lls: &[EMProb]) -> (Vec<EMProb>, EMProb) {
    let m = lls.iter().copied().fold(EMProb::NEG_INFINITY, EMProb::max);
    (lls.iter().map(|&l| (l - m).exp()).collect(), m)
}

/// Assembles a [`CsrLikelihood`] one read at a time from log-likelihoods, whether they
/// come from the aligner or from a report's ALIGNMENT_LIKELIHOODS.
pub struct CsrBuilder {
    reads: Vec<ReadIdx>,
    row_of: Vec<u32>,
    row_ptr: Vec<usize>,
    cols_seq: Vec<SeqId>,
    seen: HashSet<SeqId>,
    lik: Vec<EMProb>,
    log_scale: Vec<EMProb>,
    lls: Vec<EMProb>,
}

impl CsrBuilder {
    /// `n_reads` is a hint for the read → row table: an untrusted count (a report's
    /// `NREADS`) is only preallocated up to a cap, and a larger [`ReadIdx`] grows it.
    pub fn new(n_reads: usize) -> Self {
        Self {
            reads: Vec::new(),
            row_of: vec![NO_ROW; n_reads.min(MAX_ROW_PREALLOC)],
            row_ptr: vec![0],
            cols_seq: Vec::new(),
            seen: HashSet::new(),
            lik: Vec::new(),
            log_scale: Vec::new(),
            lls: Vec::new(),
        }
    }

    /// Reserve for the aligner's known sizes.
    fn with_capacity(n_reads: usize, nnz: usize) -> Self {
        let mut b = Self::new(0);
        b.row_of = vec![NO_ROW; n_reads];
        b.reads.reserve(n_reads);
        b.row_ptr.reserve(n_reads);
        b.cols_seq.reserve(nnz);
        b.lik.reserve(nnz);
        b.log_scale.reserve(n_reads);
        b
    }

    /// One read's `(sequence, log-likelihood)` entries in ascending sequence order;
    /// reads arrive in ascending [`ReadIdx`] order. Non-finite entries are dropped; a
    /// read left without any gets no row.
    pub fn push(&mut self, read_idx: ReadIdx, entries: impl IntoIterator<Item = (SeqId, EMProb)>) {
        debug_assert!(
            self.reads.last().is_none_or(|last| last.0 < read_idx.0),
            "reads must be pushed in ascending order"
        );
        self.lls.clear();
        for (seq, ll) in entries {
            if ll.is_finite() {
                self.cols_seq.push(seq);
                self.seen.insert(seq);
                self.lls.push(ll);
            }
        }
        if self.lls.is_empty() {
            return;
        }
        let (scaled, m) = scaled_row(&self.lls);
        self.lik.extend(scaled);
        self.log_scale.push(m);
        if read_idx.0 >= self.row_of.len() {
            self.row_of.resize(read_idx.0 + 1, NO_ROW);
        }
        self.row_of[read_idx.0] =
            u32::try_from(self.reads.len()).expect("more than u32::MAX reads with alignments");
        self.reads.push(read_idx);
        self.row_ptr.push(self.cols_seq.len());
    }

    /// Columns are the sorted, deduplicated sequences seen.
    pub fn finish(self) -> CsrLikelihood {
        let mut refs: Vec<SeqId> = self.seen.into_iter().collect();
        refs.sort_unstable_by_key(|r| r.0);
        let ref_compact: HashMap<SeqId, u32> = refs
            .iter()
            .enumerate()
            .map(|(i, r)| (*r, i as u32))
            .collect();
        let col = self.cols_seq.iter().map(|r| ref_compact[r]).collect();
        CsrLikelihood::from_rows(
            refs,
            ref_compact,
            self.reads,
            self.row_of,
            self.row_ptr,
            col,
            self.lik,
            self.log_scale,
        )
    }
}

impl CsrLikelihood {
    /// Build the matrix from the surviving alignments: one column per sequence with at
    /// least one finite log-likelihood, one row per read with such an alignment.
    pub fn build(aligns: &ReadAlignments<'_>) -> Self {
        let mut b = CsrBuilder::with_capacity(aligns.len(), aligns.nnz());
        for (read_idx, _, refs, als) in aligns.iter() {
            b.push(
                read_idx,
                refs.iter()
                    .zip(als)
                    .map(|(r, a)| (*r, *a.get_full_match_ll() as EMProb)),
            );
        }
        b.finish()
    }

    /// Assemble the CSC copy from the CSR parts.
    #[allow(clippy::too_many_arguments)]
    fn from_rows(
        refs: Vec<SeqId>,
        ref_compact: HashMap<SeqId, u32>,
        reads: Vec<ReadIdx>,
        row_of: Vec<u32>,
        row_ptr: Vec<usize>,
        col: Vec<u32>,
        lik: Vec<EMProb>,
        log_scale: Vec<EMProb>,
    ) -> Self {
        let n_refs = refs.len();
        let mut col_ptr = vec![0usize; n_refs + 1];
        for &c in &col {
            col_ptr[c as usize + 1] += 1;
        }
        for s in 0..n_refs {
            col_ptr[s + 1] += col_ptr[s];
        }
        let mut fill = col_ptr[..n_refs].to_vec();
        let mut row_idx = vec![0u32; col.len()];
        let mut lik_t = vec![0.0; col.len()];
        for r in 0..reads.len() {
            for k in row_ptr[r]..row_ptr[r + 1] {
                let s = col[k] as usize;
                row_idx[fill[s]] = r as u32;
                lik_t[fill[s]] = lik[k];
                fill[s] += 1;
            }
        }
        Self {
            refs,
            ref_compact,
            reads,
            row_of,
            row_ptr,
            col,
            lik,
            log_scale,
            col_ptr,
            row_idx,
            lik_t,
        }
    }

    pub fn n_reads(&self) -> usize {
        self.reads.len()
    }

    pub fn n_refs(&self) -> usize {
        self.refs.len()
    }

    /// The sequences with a column, in column order.
    pub fn refs(&self) -> &[SeqId] {
        &self.refs
    }

    /// The column of a sequence, if it has one.
    pub fn col_of(&self, seq: SeqId) -> Option<usize> {
        self.ref_compact.get(&seq).map(|&c| c as usize)
    }

    /// Position of the `(read, sequence)` entry
    #[cfg(test)]
    fn entry(&self, read_idx: ReadIdx, ref_idx: SeqId) -> Option<usize> {
        let range = self.row(read_idx)?;
        let c = *self.ref_compact.get(&ref_idx)?;
        let row = &self.col[range.clone()];
        row.binary_search(&c).ok().map(|i| range.start + i)
    }

    /// Per-column starting proportions
    pub fn initial_pi(&self) -> Vec<EMProb> {
        let best: Vec<(EMProb, EMProb)> = (0..self.n_reads())
            .into_par_iter()
            .map(|r| {
                let row = &self.lik[self.row_ptr[r]..self.row_ptr[r + 1]];
                let b = row.iter().copied().fold(EMProb::NEG_INFINITY, EMProb::max);
                let tied = row.iter().filter(|&&l| l == b).count();
                (b, 1.0 / tied as EMProb)
            })
            .collect();
        let counts: Vec<EMProb> = (0..self.n_refs())
            .into_par_iter()
            .map(|s| {
                (self.col_ptr[s]..self.col_ptr[s + 1])
                    .filter(|&k| self.lik_t[k] == best[self.row_idx[k] as usize].0)
                    .map(|k| best[self.row_idx[k] as usize].1)
                    .sum::<EMProb>()
            })
            .collect();
        let total: EMProb = counts.iter().sum();
        counts.iter().map(|&c| c / total).collect()
    }

    /// Normalizer of every read under per-column proportions `pi`.
    fn row_normalizers(&self, pi: &[EMProb]) -> Vec<EMProb> {
        (0..self.n_reads())
            .into_par_iter()
            .map(|r| {
                (self.row_ptr[r]..self.row_ptr[r + 1])
                    .map(|k| self.lik[k] * pi[self.col[k] as usize])
                    .sum::<EMProb>()
            })
            .collect()
    }

    /// E-step
    fn e_step(&self, pi: &[EMProb]) -> (EMProb, Vec<EMProb>, Vec<EMProb>) {
        let d = self.row_normalizers(pi);
        let ll: EMProb = d
            .iter()
            .zip(&self.log_scale)
            .filter(|(&di, _)| di > 0.0)
            .map(|(&di, &m)| di.ln() + m)
            .sum();
        let inv: Vec<EMProb> = d
            .iter()
            .map(|&di| if di > 0.0 { 1.0 / di } else { 0.0 })
            .collect();
        let m: Vec<EMProb> = (0..self.n_refs())
            .into_par_iter()
            .map(|s| {
                if pi[s] == 0.0 {
                    return 0.0;
                }
                pi[s]
                    * (self.col_ptr[s]..self.col_ptr[s + 1])
                        .map(|k| self.lik_t[k] * inv[self.row_idx[k] as usize])
                        .sum::<EMProb>()
            })
            .collect();
        (ll, m, d)
    }

    /// The row of a read, if it has one.
    #[cfg(test)]
    fn row(&self, read_idx: ReadIdx) -> Option<std::ops::Range<usize>> {
        let row = *self.row_of.get(read_idx.0)?;
        if row == NO_ROW {
            return None;
        }
        Some(self.row_ptr[row as usize]..self.row_ptr[row as usize + 1])
    }
}

/// Posteriors under the final proportions
#[cfg(test)]
pub struct Posteriors<'m> {
    pub(crate) matrix: &'m CsrLikelihood,
    /// `π̃_s` per column.
    pub(crate) pi: Vec<EMProb>,
    /// `D_i` per row.
    pub(crate) d: Vec<EMProb>,
}

#[cfg(test)]
impl Posteriors<'_> {
    fn weight(&self, k: usize, row: u32) -> EMProb {
        let d = self.d[row as usize];
        if d > 0.0 {
            self.pi[self.matrix.col[k] as usize] * self.matrix.lik[k] / d
        } else {
            0.0
        }
    }

    /// Posterior of one sequence for a read.
    pub fn get(&self, read_idx: ReadIdx, ref_idx: SeqId) -> EMProb {
        match (
            self.matrix.entry(read_idx, ref_idx),
            self.matrix.row_of.get(read_idx.0),
        ) {
            (Some(k), Some(&row)) => self.weight(k, row),
            _ => 0.0,
        }
    }

    /// The leaf owning a sequence, or [`NO_NODE`] when it has no column or was pruned.
    pub fn leaf_of(&self, ref_idx: SeqId, part: &Partition) -> NodeID {
        self.matrix
            .ref_compact
            .get(&ref_idx)
            .map_or(NO_NODE, |&c| part.leaf_of(c as usize))
    }

    /// A read's posterior per leaf of the tree, summed over the leaf's sequences
    pub fn per_leaf(&self, read_idx: ReadIdx, part: &Partition) -> Vec<(NodeID, EMProb)> {
        let Some(range) = self.matrix.row(read_idx) else {
            return Vec::new();
        };
        let row = self.matrix.row_of[read_idx.0];
        let mut hits: Vec<(NodeID, EMProb)> = Vec::new();
        for k in range {
            let leaf = part.leaf_of(self.matrix.col[k] as usize);
            if leaf == NO_NODE {
                continue;
            }
            let w = self.weight(k, row);
            if w == 0.0 {
                continue;
            }
            match hits.iter_mut().find(|(l, _)| *l == leaf) {
                Some((_, p)) => *p += w,
                None => hits.push((leaf, w)),
            }
        }
        hits
    }
}

/// How the conditionals start.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Init {
    /// Each read's best-scoring sequences share one count (ties split evenly); the counts
    /// summed up the tree give `q_{c|v} = N_c / N_v`.
    #[default]
    Split,
    /// `q_{c|v} = 1 / |ch(v)|` at every internal vertex.
    Uniform,
}

impl Init {
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "split" => Some(Self::Split),
            "uniform" => Some(Self::Uniform),
            _ => None,
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::Split => "split",
            Self::Uniform => "uniform",
        }
    }
}

/// Iterations the EM always runs (up to `num_iter`) before the convergence threshold may
/// stop it.
pub const MIN_EM_ITER: usize = 50;

/// Settings of [`estimate_proportions`].
pub struct EmParams {
    pub init: Init,
    pub num_iter: usize,
    /// Penalty weight; also the minimum expected reads a vertex needs to survive. 0 = plain EM.
    pub rho: EMProb,
    pub omega: EMProb,
    /// Relative change in the penalized log-likelihood that stops the iteration, once
    /// [`MIN_EM_ITER`] iterations have run.
    pub em_threshold: EMProb,
}

/// Result of [`estimate_proportions`].
pub struct EmFit {
    /// `π_v` per node id of the (pruned) tree.
    pub node_pi: Vec<EMProb>,
    /// `π̃_s` per column under the final tree.
    pub pi_tilde: Vec<EMProb>,
    /// Penalized log-likelihood after each iteration (first entry: the initial state).
    pub ll_trace: Vec<EMProb>,
    /// `(leaf, π_ℓ)` over the leaves of the tree at each iteration, aligned with
    /// `ll_trace`: entry 0 is the initial state, the last the final proportions.
    pub pi_trace: Vec<Vec<(NodeID, EMProb)>>,
    pub pruned: usize,
    pub collapsed: usize,
    /// Every child of the root fell below `ρ` at some iteration: the sample is one tip.
    pub root_collapsed: bool,
}

/// Run the penalized EM over `tree`, pruning it in place.
pub fn estimate_proportions(
    csr: &CsrLikelihood,
    tree: &mut TaxTree,
    part: &mut Partition,
    p: &EmParams,
    progress: Option<Arc<QueryProgress>>,
) -> EmFit {
    if csr.n_refs() == 0 {
        return EmFit {
            node_pi: Vec::new(),
            pi_tilde: Vec::new(),
            ll_trace: Vec::new(),
            pi_trace: Vec::new(),
            pruned: 0,
            collapsed: 0,
            root_collapsed: false,
        };
    }
    match p.init {
        Init::Split => init_conditionals(tree, part, &csr.initial_pi()),
        Init::Uniform => init_uniform(tree),
    }

    if let Some(ref pr) = progress {
        pr.phase.store(2, Ordering::Relaxed);
        pr.em_iter_total.store(p.num_iter as u64, Ordering::Relaxed);
        pr.em_iter_done.store(0, Ordering::Relaxed);
    }
    let pb = ProgressBar::with_draw_target(Some(p.num_iter as u64), ProgressDrawTarget::stderr());
    pb.set_style(ProgressStyle::with_template("Running EM: {spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] iteration {pos}/{len} ({eta}) ({msg})").unwrap());

    let penalized = |csr: &CsrLikelihood, tree: &TaxTree, part: &Partition| -> EMProb {
        let (ll, _, _) = csr.e_step(&contig_pi(part, &leaf_marginals(tree)));
        ll - penalty(tree, p.rho, p.omega)
    };
    let mut ll_trace = vec![penalized(csr, tree, part)];
    let mut pi_trace: Vec<Vec<(NodeID, EMProb)>> = Vec::new();
    let (mut n_pruned, mut n_collapsed, mut root_collapsed) = (0, 0, false);

    for i in 0..p.num_iter {
        let node_pi = leaf_marginals(tree);
        pi_trace.push(leaf_proportions(tree, &node_pi));
        let pi = contig_pi(part, &node_pi);
        let (ll, m, _) = csr.e_step(&pi);
        subtree_counts(tree, part, &m);
        let report = if p.rho > 0.0 {
            prune_and_collapse(tree, part, p.rho)
        } else {
            Default::default()
        };
        n_pruned += report.pruned;
        n_collapsed += report.collapsed;
        root_collapsed |= report.root_collapsed;
        m_step(tree, p.rho, p.omega);

        let lp = ll - penalty(tree, p.rho, p.omega);
        let prev = *ll_trace.last().expect("seeded");
        ll_trace.push(lp);
        let diff = lp - prev;
        pb.set_message(format!("{diff:.3e}"));
        pb.inc(1);
        if let Some(ref pr) = progress {
            pr.em_iter_done.store((i + 1) as u64, Ordering::Relaxed);
        }
        let scale = lp.abs().max(1.0);
        let converged = !report.changed() && diff.abs() <= p.em_threshold * scale;
        if i + 1 >= MIN_EM_ITER && converged {
            break;
        }
    }
    pb.finish_with_message(format!(
        "{} iterations, final penalized LL: {}",
        ll_trace.len() - 1,
        ll_trace.last().unwrap()
    ));

    let node_pi = leaf_marginals(tree);
    pi_trace.push(leaf_proportions(tree, &node_pi));
    let pi_tilde = contig_pi(part, &node_pi);
    EmFit {
        node_pi,
        pi_tilde,
        ll_trace,
        pi_trace,
        pruned: n_pruned,
        collapsed: n_collapsed,
        root_collapsed,
    }
}

/// `(leaf, π_ℓ)` for every leaf with mass, by node id.
fn leaf_proportions(tree: &TaxTree, node_pi: &[EMProb]) -> Vec<(NodeID, EMProb)> {
    let mut rows: Vec<(NodeID, EMProb)> = tree
        .get_leaf_ids()
        .filter_map(|l| {
            let p = node_pi.get(l).copied().unwrap_or(0.0);
            (p > 0.0).then_some((l, p))
        })
        .collect();
    rows.sort_unstable_by_key(|r| r.0);
    rows
}

/// Assign every read with a row to the deepest vertex whose subtree posterior is at least
/// `gamma`
#[cfg(test)]
pub fn assign_reads(
    csr: &CsrLikelihood,
    tree: &TaxTree,
    part: &Partition,
    posteriors: &Posteriors<'_>,
    gamma: EMProb,
) -> HashMap<ReadIdx, (NodeID, EMProb)> {
    csr.reads
        .par_iter()
        .filter_map(|&read_idx| {
            let hits = posteriors.per_leaf(read_idx, part);
            if hits.is_empty() {
                return None;
            }
            Some((read_idx, assign_read(tree, &hits, gamma)))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tree::{label, TaxNode};

    /// A matrix from dense rows of linear likelihoods (0.0 = no entry); sequence j is
    /// `SeqId(j)`. Rows are stored unscaled (`log_scale` 0).
    fn csr(rows: &[&[EMProb]]) -> CsrLikelihood {
        let logs: Vec<Vec<Option<EMProb>>> = rows
            .iter()
            .map(|r| r.iter().map(|&l| (l != 0.0).then(|| l.ln())).collect())
            .collect();
        let m = csr_log(&logs);
        // undo the scaling so the tests can reason in linear values
        let mut lik = Vec::with_capacity(m.lik.len());
        for r in 0..m.n_reads() {
            for k in m.row_ptr[r]..m.row_ptr[r + 1] {
                lik.push(m.lik[k] * m.log_scale[r].exp());
            }
        }
        let n = m.n_reads();
        CsrLikelihood::from_rows(
            m.refs,
            m.ref_compact,
            m.reads,
            m.row_of,
            m.row_ptr,
            m.col,
            lik,
            vec![0.0; n],
        )
    }

    /// A matrix from dense rows of log-likelihoods (`None` = no entry), scaled like `build`.
    fn csr_log(rows: &[Vec<Option<EMProb>>]) -> CsrLikelihood {
        let n_refs = rows[0].len();
        let refs: Vec<SeqId> = (0..n_refs as u32).map(SeqId).collect();
        let ref_compact = refs
            .iter()
            .enumerate()
            .map(|(i, r)| (*r, i as u32))
            .collect();
        let (mut row_ptr, mut col, mut lik, mut log_scale) =
            (vec![0], Vec::new(), Vec::new(), Vec::new());
        for row in rows {
            let mut lls = Vec::new();
            for (j, l) in row.iter().enumerate() {
                if let Some(l) = l {
                    col.push(j as u32);
                    lls.push(*l);
                }
            }
            let (scaled, m) = scaled_row(&lls);
            lik.extend(scaled);
            log_scale.push(m);
            row_ptr.push(col.len());
        }
        CsrLikelihood::from_rows(
            refs,
            ref_compact,
            (0..rows.len()).map(ReadIdx).collect(),
            (0..rows.len() as u32).collect(),
            row_ptr,
            col,
            lik,
            log_scale,
        )
    }

    fn node(id: usize, label: &str) -> TaxNode {
        let mut n = TaxNode::new(id);
        n.set_taxa(Some(label.to_string()));
        n.set_annotation(Some("-".into()));
        n.set_zeta(Some(0.0));
        n
    }

    /// Root 0 with one tip per sequence (tip id = column + 1).
    fn star(n_cols: usize) -> (TaxTree, Partition) {
        let mut t = TaxTree::new(0);
        t.set_node(node(0, "root"));
        for j in 0..n_cols {
            t.add_child(0, node(j + 1, &format!("c{j}")));
        }
        (t, Partition::new((1..=n_cols).collect()))
    }

    /// Root 0 → A(1){c0 (3), c1 (4)}, B(2){c2 (5)}.
    fn two_genera() -> (TaxTree, Partition) {
        let mut t = TaxTree::new(0);
        t.set_node(node(0, "root"));
        t.add_child(0, node(1, "A"));
        t.add_child(0, node(2, "B"));
        t.add_child(1, node(3, "sp1"));
        t.add_child(1, node(4, "sp2"));
        t.add_child(2, node(5, "sp3"));
        (t, Partition::new(vec![3, 4, 5]))
    }

    fn params(rho: EMProb) -> EmParams {
        EmParams {
            init: Init::Split,
            num_iter: 50,
            rho,
            omega: 1e-10,
            em_threshold: 1e-9,
        }
    }

    fn leaf_props(part: &Partition, pi: &[EMProb]) -> Vec<EMProb> {
        (0..part.n_cols())
            .map(|c| {
                let l = part.leaf_of(c);
                if l == NO_NODE {
                    0.0
                } else {
                    pi[l]
                }
            })
            .collect()
    }

    fn posteriors<'m>(m: &'m CsrLikelihood, fit: &EmFit) -> Posteriors<'m> {
        Posteriors {
            matrix: m,
            pi: fit.pi_tilde.clone(),
            d: m.row_normalizers(&fit.pi_tilde),
        }
    }

    #[test]
    fn csc_matches_the_dense_e_step() {
        let rows: Vec<&[EMProb]> = vec![
            &[0.9, 0.1, 0.0, 0.0],
            &[0.0, 0.4, 0.6, 0.0],
            &[0.2, 0.0, 0.3, 0.5],
        ];
        let m = csr(&rows);
        assert_eq!(m.col_ptr, vec![0, 2, 4, 6, 7]);
        assert_eq!(m.row_idx, vec![0, 2, 0, 1, 1, 2, 2]);
        let pi = [0.1, 0.2, 0.3, 0.4];
        let (ll, counts, d) = m.e_step(&pi);
        let mut want_ll = 0.0;
        let mut want = vec![0.0; 4];
        for (i, row) in rows.iter().enumerate() {
            let di: EMProb = row.iter().zip(pi).map(|(a, p)| a * p).sum();
            assert!((d[i] - di).abs() < 1e-15);
            want_ll += di.ln();
            for (s, a) in row.iter().enumerate() {
                want[s] += a * pi[s] / di;
            }
        }
        assert!((ll - want_ll).abs() < 1e-12);
        for (c, w) in counts.iter().zip(want) {
            assert!((c - w).abs() < 1e-12, "{counts:?}");
        }
        assert_eq!(m.entry(ReadIdx(2), SeqId(2)), Some(5));
        assert_eq!(m.entry(ReadIdx(0), SeqId(3)), None);
    }

    #[test]
    fn e_step_is_bitwise_identical_across_thread_counts() {
        let rows: Vec<Vec<EMProb>> = (0..500)
            .map(|i| {
                (0..7)
                    .map(|j| {
                        if (i + j) % 3 == 0 {
                            0.0
                        } else {
                            1.0 / (1.0 + ((i * 7 + j) % 11) as EMProb)
                        }
                    })
                    .collect()
            })
            .collect();
        let refs: Vec<&[EMProb]> = rows.iter().map(Vec::as_slice).collect();
        let m = csr(&refs);
        let pi: Vec<EMProb> = (1..=7).map(|j| j as EMProb / 28.0).collect();
        let run = |threads: usize| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| m.e_step(&pi))
        };
        let (ll1, m1, d1) = run(1);
        let (ll4, m4, d4) = run(4);
        assert_eq!(ll1.to_bits(), ll4.to_bits());
        assert_eq!(m1, m4);
        assert_eq!(d1, d4);
        assert_eq!(run(1).1, run(3).1);
    }

    #[test]
    fn row_scaling_survives_likelihoods_that_underflow_linear_space() {
        // both alignments at log-likelihood -800: exp underflows to 0 in f64
        let m = csr_log(&[
            vec![Some(-800.0), Some(-800.0)],
            vec![Some(-900.0), Some(-901.0)],
        ]);
        assert_eq!(m.lik, vec![1.0, 1.0, 1.0, (-1.0_f64).exp()]);
        let (ll, counts, d) = m.e_step(&[0.5, 0.5]);
        assert!(d.iter().all(|&x| x > 0.0));
        let want = -800.0 + (0.5 + 0.5_f64).ln() + (-900.0) + (0.5 + 0.5 * (-1.0_f64).exp()).ln();
        assert!((ll - want).abs() < 1e-9, "{ll} vs {want}");
        let r1 = 1.0 / (1.0 + (-1.0_f64).exp());
        assert!((counts[0] - (0.5 + r1)).abs() < 1e-12, "{counts:?}");
        assert!((counts[1] - (0.5 + 1.0 - r1)).abs() < 1e-12);
        assert_eq!(m.initial_pi(), vec![0.75, 0.25]);
    }

    #[test]
    fn star_em_drops_ref_with_fewer_expected_reads_than_rho() {
        let mut rows: Vec<&[EMProb]> = Vec::new();
        rows.extend(std::iter::repeat_n(&[1.0, 0.0, 0.0][..], 300));
        rows.extend(std::iter::repeat_n(&[0.0, 1.0, 0.0][..], 300));
        // Initially won by sequence 2, but only 100 reads: below rho.
        rows.extend(std::iter::repeat_n(&[0.0, 0.4, 0.6][..], 100));
        let m = csr(&rows);
        let (mut t, mut part) = star(3);

        let fit = estimate_proportions(&m, &mut t, &mut part, &params(150.0), None);
        assert_eq!(part.leaf_of(2), NO_NODE);
        assert!(fit.pruned >= 1);
        assert!(!fit.root_collapsed);
        let props = leaf_props(&part, &fit.node_pi);
        assert_eq!(props[2], 0.0);
        assert!(
            (props.iter().sum::<EMProb>() - 1.0).abs() < 1e-9,
            "{props:?}"
        );
        assert_eq!(fit.pi_tilde[2], 0.0);
        let post = posteriors(&m, &fit);
        assert_eq!(post.get(ReadIdx(650), SeqId(2)), 0.0);
        assert!((post.get(ReadIdx(650), SeqId(1)) - 1.0).abs() < 1e-12);
        let a = assign_reads(&m, &t, &part, &post, 0.6);
        assert!((600..700).all(|r| a[&ReadIdx(r)].0 == part.leaf_of(1)));
    }

    #[test]
    fn pi_trace_tracks_leaves_per_iteration() {
        let mut rows: Vec<&[EMProb]> = Vec::new();
        rows.extend(std::iter::repeat_n(&[1.0, 0.0, 0.0][..], 300));
        rows.extend(std::iter::repeat_n(&[0.0, 1.0, 0.0][..], 300));
        rows.extend(std::iter::repeat_n(&[0.0, 0.4, 0.6][..], 100));
        let m = csr(&rows);
        let (mut t, mut part) = star(3);
        let fit = estimate_proportions(&m, &mut t, &mut part, &params(150.0), None);

        assert_eq!(fit.pi_trace.len(), fit.ll_trace.len());
        assert!(fit.pi_trace.len() >= 2);
        // Entry 0: every initial leaf with mass; the last: the final leaves == node_pi.
        let first: Vec<NodeID> = fit.pi_trace[0].iter().map(|r| r.0).collect();
        assert_eq!(first, vec![1, 2, 3]);
        let last = fit.pi_trace.last().unwrap();
        let final_leaves: Vec<NodeID> = {
            let mut l: Vec<NodeID> = t.get_leaf_ids().collect();
            l.sort_unstable();
            l
        };
        assert_eq!(last.iter().map(|r| r.0).collect::<Vec<_>>(), final_leaves);
        for &(leaf, p) in last {
            assert_eq!(p, fit.node_pi[leaf]);
        }
        for entry in &fit.pi_trace {
            let total: EMProb = entry.iter().map(|r| r.1).sum();
            assert!((total - 1.0).abs() < 1e-9, "{entry:?}");
            assert!(entry.windows(2).all(|w| w[0].0 < w[1].0));
        }
        // Pruning only ever shrinks the leaf set.
        for w in fit.pi_trace.windows(2) {
            assert!(w[1].len() <= w[0].len());
        }
        assert!(last.len() < first.len());
    }

    /// The streaming `.matches` writer, fed posteriors recomputed from the stored
    /// log-likelihoods, must make exactly the decision `assign_reads` makes from the
    /// EM's own posteriors: same vertex, same `h` to the bit.
    #[test]
    fn matches_writer_agrees_with_assign_reads() {
        use crate::report::posterior::row_posteriors;
        use crate::report::tables::{Candidate, MatchesWriter};
        use crate::report::{AlignRow, ReportTree};
        use crate::tree::{is_synthetic_root, label};

        let mut rows: Vec<&[EMProb]> = Vec::new();
        rows.extend(std::iter::repeat_n(&[1.0, 0.0, 0.0][..], 40));
        rows.extend(std::iter::repeat_n(&[0.0, 1.0, 0.0][..], 30));
        rows.extend(std::iter::repeat_n(&[0.6, 0.4, 0.0][..], 20));
        rows.extend(std::iter::repeat_n(&[0.0, 0.5, 0.5][..], 20));
        rows.extend(std::iter::repeat_n(&[0.3, 0.3, 0.4][..], 10));
        rows.extend(std::iter::repeat_n(&[0.5, 0.5, 0.0][..], 10));
        // Scaled like `build`, so the EM's internal posteriors use the same arithmetic.
        let logs: Vec<Vec<Option<EMProb>>> = rows
            .iter()
            .map(|r| r.iter().map(|&l| (l != 0.0).then(|| l.ln())).collect())
            .collect();
        let m = csr_log(&logs);
        let (mut t, mut part) = two_genera();
        let stored = t.clone();
        let tips: Vec<NodeID> = (0..part.n_cols()).map(|c| part.leaf_of(c)).collect();
        let fit = estimate_proportions(&m, &mut t, &mut part, &params(12.0), None);
        let post = posteriors(&m, &fit);
        let report_tree = ReportTree::from_fit(&stored, &tips, &t, &part, &fit);
        let pi_tilde = |s: SeqId| fit.pi_tilde[s.index()];
        for gamma in [0.5, 0.6, 0.9, 0.99] {
            let oracle = assign_reads(&m, &t, &part, &post, gamma);
            let mut w = MatchesWriter::new(
                Vec::new(),
                "EM test",
                "RefID",
                report_tree.to_tax_tree(),
                gamma,
            )
            .unwrap();
            let mut posts = Vec::new();
            for (r, row) in logs.iter().enumerate() {
                let aligns: Vec<AlignRow> = row
                    .iter()
                    .enumerate()
                    .filter_map(|(s, l)| {
                        l.map(|loglik| AlignRow {
                            seq: SeqId(s as u32),
                            loglik,
                            pos: (r, s),
                        })
                    })
                    .collect();
                row_posteriors(pi_tilde, &aligns, &mut posts);
                let cands: Vec<Candidate> = aligns
                    .iter()
                    .zip(&posts)
                    .filter(|(a, _)| post.leaf_of(a.seq, &part) != NO_NODE)
                    .map(|(a, &posterior)| Candidate {
                        leaf: post.leaf_of(a.seq, &part),
                        posterior,
                        loglik: a.loglik,
                        pos: a.pos,
                    })
                    .collect();
                // The recomputed posterior is the EM's, bit for bit.
                for (a, &posterior) in aligns.iter().zip(&posts) {
                    assert_eq!(
                        posterior,
                        post.get(ReadIdx(r), a.seq),
                        "read {r} seq {:?}",
                        a.seq
                    );
                }
                w.push(&format!("r{r}"), &cands).unwrap();
            }
            let out = String::from_utf8(w.finish().unwrap()).unwrap();
            // The `#` comment line, then the column header.
            for (r, line) in out.lines().skip(2).enumerate() {
                let cols: Vec<&str> = line.split('\t').collect();
                match oracle.get(&ReadIdx(r)) {
                    Some(&(v, h)) if !is_synthetic_root(&t, v) => {
                        assert_eq!(cols[1], label(&t, v), "gamma {gamma} read {r}: {line}");
                        assert_eq!(cols[2], format!("{:.5e}", h), "gamma {gamma} read {r}");
                    }
                    _ => assert_eq!(cols[1], "unclassified", "gamma {gamma} read {r}: {line}"),
                }
            }
        }
    }

    #[test]
    fn csr_builder_skips_empty_reads_and_non_finite_entries() {
        // A huge announced count is not preallocated; the table grows with the rows.
        let b = CsrBuilder::new(usize::MAX / 2);
        assert_eq!(b.row_of.len(), MAX_ROW_PREALLOC);
        let mut b2 = CsrBuilder::new(2);
        b2.push(ReadIdx(5), [(SeqId(1), -1.0)]);
        assert_eq!(b2.finish().row_of.len(), 6);
        let mut b = CsrBuilder::new(2);
        b.push(ReadIdx(0), [(SeqId(7), -2.0), (SeqId(9), -1.0)]);
        b.push(ReadIdx(1), [(SeqId(9), EMProb::NEG_INFINITY)]);
        b.push(ReadIdx(3), [(SeqId(2), EMProb::NAN), (SeqId(7), -5.0)]);
        let m = b.finish();
        assert_eq!(m.refs(), &[SeqId(7), SeqId(9)]);
        assert_eq!(m.n_reads(), 2);
        assert_eq!(m.reads, vec![ReadIdx(0), ReadIdx(3)]);
        assert_eq!(m.row_of, vec![0, NO_ROW, NO_ROW, 1]);
        assert_eq!(m.row_ptr, vec![0, 2, 3]);
        assert_eq!(m.col, vec![0, 1, 0]);
        assert_eq!(m.log_scale, vec![-1.0, -5.0]);
        assert_eq!(m.lik, vec![(-1.0f64).exp(), 1.0, 1.0]);
        assert_eq!(m.col_ptr, vec![0, 2, 3]);
        assert_eq!(m.row_idx, vec![0, 1, 0]);
    }

    #[test]
    fn convergence_is_checked_only_after_the_minimum_iterations() {
        // Unambiguous reads: the split init is already the fixed point, so the
        // likelihood stops changing at once.
        let mut rows: Vec<&[EMProb]> = Vec::new();
        rows.extend(std::iter::repeat_n(&[1.0, 0.0][..], 30));
        rows.extend(std::iter::repeat_n(&[0.0, 1.0][..], 10));
        let m = csr(&rows);
        let run = |num_iter: usize| {
            let (mut t, mut part) = star(2);
            let p = EmParams {
                num_iter,
                em_threshold: 1e-3,
                ..params(0.0)
            };
            estimate_proportions(&m, &mut t, &mut part, &p, None)
                .ll_trace
                .len()
                - 1
        };
        assert_eq!(run(500), MIN_EM_ITER);
        // --iter stays the cap below the minimum.
        assert_eq!(run(20), 20);
    }

    #[test]
    fn rho_zero_star_em_matches_count_ratio() {
        let mut rows: Vec<&[EMProb]> = Vec::new();
        rows.extend(std::iter::repeat_n(&[1.0, 0.0][..], 30));
        rows.extend(std::iter::repeat_n(&[0.0, 1.0][..], 10));
        let m = csr(&rows);
        let (mut t, mut part) = star(2);
        let fit = estimate_proportions(&m, &mut t, &mut part, &params(0.0), None);
        let props = leaf_props(&part, &fit.node_pi);
        assert!(
            (props[0] - 0.75).abs() < 1e-9 && (props[1] - 0.25).abs() < 1e-9,
            "{props:?}"
        );
        assert_eq!((fit.pruned, fit.collapsed), (0, 0));
        // the trace is non-decreasing for the plain EM
        assert!(fit.ll_trace.windows(2).all(|w| w[1] >= w[0] - 1e-9));
    }

    #[test]
    fn hierarchical_em_prunes_weak_sibling() {
        let mut rows: Vec<&[EMProb]> = Vec::new();
        rows.extend(std::iter::repeat_n(&[1.0, 0.0, 0.0][..], 300));
        rows.extend(std::iter::repeat_n(&[0.0, 1.0, 0.0][..], 5));
        rows.extend(std::iter::repeat_n(&[0.0, 0.0, 1.0][..], 300));
        let m = csr(&rows);
        let (mut t, mut part) = two_genera();
        let fit = estimate_proportions(&m, &mut t, &mut part, &params(50.0), None);
        assert_eq!(part.leaf_of(1), NO_NODE, "sp2 pruned");
        assert!(t.get_node(4).is_none());
        assert!(t.get_node(1).is_some() && t.get_node(3).is_some());
        assert!((fit.node_pi[1] - 0.5).abs() < 1e-6, "A {}", fit.node_pi[1]);
        assert!((fit.node_pi[2] - 0.5).abs() < 1e-6, "B {}", fit.node_pi[2]);
        assert_eq!(label(&t, 1), "A");
        // the 5 reads of sp2 keep a row but no mass: unassigned
        let post = posteriors(&m, &fit);
        assert!(post.per_leaf(ReadIdx(302), &part).is_empty());
        let a = assign_reads(&m, &t, &part, &post, 0.6);
        assert!(!a.contains_key(&ReadIdx(302)));
        assert_eq!(a[&ReadIdx(0)].0, 3);
    }

    #[test]
    fn collapse_when_all_children_weak() {
        let mut rows: Vec<&[EMProb]> = Vec::new();
        rows.extend(std::iter::repeat_n(&[1.0, 0.0, 0.0, 0.0][..], 20));
        rows.extend(std::iter::repeat_n(&[0.0, 1.0, 0.0, 0.0][..], 20));
        rows.extend(std::iter::repeat_n(&[0.0, 0.0, 1.0, 0.0][..], 20));
        rows.extend(std::iter::repeat_n(&[0.0, 0.0, 0.0, 1.0][..], 300));
        let m = csr(&rows);
        // root 0 → A(1){3,4,5}, B(2){6}
        let mut t = TaxTree::new(0);
        t.set_node(node(0, "root"));
        t.add_child(0, node(1, "A"));
        t.add_child(0, node(2, "B"));
        for (id, l) in [(3, "a1"), (4, "a2"), (5, "a3")] {
            t.add_child(1, node(id, l));
        }
        t.add_child(2, node(6, "b1"));
        let mut part = Partition::new(vec![3, 4, 5, 6]);
        let fit = estimate_proportions(&m, &mut t, &mut part, &params(50.0), None);
        assert_eq!(fit.collapsed, 1);
        assert!(t.is_leaf(1), "A collapsed into a tip");
        assert_eq!(part.leaf_of(0), 1);
        assert_eq!(part.size_of(1), 3);
        // the penalty shrinks A below its count share, but it keeps mass
        let a = fit.node_pi[1];
        assert!(a > 0.0 && a < 60.0 / 360.0, "{a}");
        assert!((a + fit.node_pi[2] - 1.0).abs() < 1e-9);
        // sequences of A share its mass evenly in the mixture
        assert!((fit.pi_tilde[0] - fit.pi_tilde[1]).abs() < 1e-12);
        assert!((fit.pi_tilde[0] - fit.node_pi[1] / 3.0).abs() < 1e-12);
        // a read of a1 is assigned to the collapsed vertex with posterior 1
        let post = posteriors(&m, &fit);
        let hits = post.per_leaf(ReadIdx(0), &part);
        assert_eq!(hits, vec![(1, 1.0)]);
        let asg = assign_reads(&m, &t, &part, &post, 0.6);
        assert_eq!(asg[&ReadIdx(0)], (1, 1.0));
    }

    #[test]
    fn uniform_init_reaches_the_same_plain_fit() {
        let mut rows: Vec<&[EMProb]> = Vec::new();
        rows.extend(std::iter::repeat_n(&[1.0, 0.0][..], 30));
        rows.extend(std::iter::repeat_n(&[0.0, 1.0][..], 10));
        let m = csr(&rows);
        let (mut t, mut part) = star(2);
        let p = EmParams {
            init: Init::Uniform,
            ..params(0.0)
        };
        let fit = estimate_proportions(&m, &mut t, &mut part, &p, None);
        let props = leaf_props(&part, &fit.node_pi);
        assert!((props[0] - 0.75).abs() < 1e-9, "{props:?}");
        assert_eq!(Init::parse("uniform"), Some(Init::Uniform));
        assert_eq!(Init::parse("nope"), None);
    }

    #[test]
    fn root_collapse_is_reported() {
        let mut rows: Vec<&[EMProb]> = Vec::new();
        rows.extend(std::iter::repeat_n(&[1.0, 0.0][..], 10));
        rows.extend(std::iter::repeat_n(&[0.0, 1.0][..], 10));
        let m = csr(&rows);
        let (mut t, mut part) = star(2);
        let fit = estimate_proportions(&m, &mut t, &mut part, &params(50.0), None);
        assert!(fit.root_collapsed);
        assert!(t.is_leaf(0));
        assert_eq!(part.size_of(0), 2);
        assert!((fit.node_pi[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn initial_pi_splits_exact_ties_evenly() {
        // two identical sequences: every read ties, so neither may start at zero
        let m = csr(&[&[0.5, 0.5], &[0.2, 0.2], &[0.9, 0.9]]);
        assert_eq!(m.initial_pi(), vec![0.5, 0.5]);
    }

    #[test]
    fn initial_pi_counts_strict_winners_and_shares_ties() {
        // read 0 -> sequence 0; read 1 tied between 1 and 2; read 2 tied three ways
        let m = csr(&[&[0.9, 0.1, 0.0], &[0.1, 0.4, 0.4], &[0.3, 0.3, 0.3]]);
        let pi = m.initial_pi();
        let third = 1.0 / 3.0;
        let want = [
            (1.0 + third) / 3.0,
            (0.5 + third) / 3.0,
            (0.5 + third) / 3.0,
        ];
        for (p, w) in pi.iter().zip(want) {
            assert!((p - w).abs() < 1e-12, "{pi:?} vs {want:?}");
        }
        assert!((pi.iter().sum::<EMProb>() - 1.0).abs() < 1e-12);
    }
}
