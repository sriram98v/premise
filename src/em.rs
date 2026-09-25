//! Abundance estimation: the EM over the read × reference
//!
//! [`CsrLikelihood`] is the sparse matrix the EM iterates on;
//! [`get_proportions_par_sparse_l1_reg`] and [`get_proportions_par_sparse`] run the
//! penalized and unpenalized EM on it; [`refit_proportions_on_classified`] re-estimates
//! the proportions over the reads that end up classified.

use crate::utils::{EMProb, QueryProgress, ReadAlignments, ReadIdx};
use crate::SeqId;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::sync::atomic::Ordering;
use std::sync::Arc;

/// Compressed-sparse-row (CSR) view of the read × reference likelihood matrix.
pub struct CsrLikelihood {
    /// Compact reference id → original [`SeqId`].
    refs: Vec<SeqId>,
    /// Original [`SeqId`] → compact reference id.
    ref_compact: HashMap<SeqId, u32>,
    /// CSR row → original [`ReadIdx`].
    reads: Vec<ReadIdx>,
    /// `row_of[read_idx]` = the row holding that read, or [`NO_ROW`].
    row_of: Vec<u32>,
    /// `row_ptr[r]..row_ptr[r + 1]` bounds row `r`'s entries in `col`/`lik`.
    row_ptr: Vec<usize>,
    /// Compact reference id for each stored entry; sorted within a row.
    col: Vec<u32>,
    /// Emission likelihood P(read | ref) for each stored entry.
    lik: Vec<EMProb>,
}

/// Marker in [`CsrLikelihood::row_of`] for reads without a row.
const NO_ROW: u32 = u32::MAX;

impl CsrLikelihood {
    /// Build the matrix from the surviving alignments
    pub fn build(aligns: &ReadAlignments<'_>) -> Self {
        let mut reads = Vec::with_capacity(aligns.len());
        let mut row_of = vec![NO_ROW; aligns.len()];
        let mut row_ptr = Vec::with_capacity(aligns.len() + 1);
        let mut cols_seq: Vec<SeqId> = Vec::with_capacity(aligns.nnz());
        let mut lik = Vec::with_capacity(aligns.nnz());
        row_ptr.push(0);
        for (read_idx, _, refs, als) in aligns.iter() {
            let start = cols_seq.len();
            for (ref_idx, align) in refs.iter().zip(als) {
                let ll = align.get_full_match_ll().exp();
                if ll.is_finite() && ll != 0.0 {
                    cols_seq.push(*ref_idx);
                    lik.push(ll as EMProb);
                }
            }
            if cols_seq.len() > start {
                row_of[read_idx.0] = reads.len() as u32;
                reads.push(read_idx);
                row_ptr.push(cols_seq.len());
            }
        }
        let mut refs = cols_seq.clone();
        refs.sort_unstable_by_key(|r| r.0);
        refs.dedup();
        let ref_compact: HashMap<SeqId, u32> = refs
            .iter()
            .enumerate()
            .map(|(i, r)| (*r, i as u32))
            .collect();
        let col = cols_seq.iter().map(|r| ref_compact[r]).collect();
        Self {
            refs,
            ref_compact,
            reads,
            row_of,
            row_ptr,
            col,
            lik,
        }
    }

    /// Drop every row whose read is not in `keep`
    fn retain_reads(&self, keep: &HashSet<ReadIdx>) -> Self {
        let mut reads = Vec::with_capacity(keep.len());
        let mut row_of = vec![NO_ROW; self.row_of.len()];
        let mut row_ptr = Vec::with_capacity(keep.len() + 1);
        let mut col = Vec::new();
        let mut lik = Vec::new();
        row_ptr.push(0);
        for (r, read_idx) in self.reads.iter().enumerate() {
            if !keep.contains(read_idx) {
                continue;
            }
            let (s, e) = (self.row_ptr[r], self.row_ptr[r + 1]);
            col.extend_from_slice(&self.col[s..e]);
            lik.extend_from_slice(&self.lik[s..e]);
            row_of[read_idx.0] = reads.len() as u32;
            reads.push(*read_idx);
            row_ptr.push(col.len());
        }
        Self {
            refs: self.refs.clone(),
            ref_compact: self.ref_compact.clone(),
            reads,
            row_of,
            row_ptr,
            col,
            lik,
        }
    }

    pub fn n_reads(&self) -> usize {
        self.reads.len()
    }

    pub fn n_refs(&self) -> usize {
        self.refs.len()
    }

    /// Position of the `(read, reference)` entry
    fn entry(&self, read_idx: ReadIdx, ref_idx: SeqId) -> Option<usize> {
        let row = *self.row_of.get(read_idx.0)?;
        if row == NO_ROW {
            return None;
        }
        let c = *self.ref_compact.get(&ref_idx)?;
        let (s, e) = (self.row_ptr[row as usize], self.row_ptr[row as usize + 1]);
        self.col[s..e].binary_search(&c).ok().map(|k| s + k)
    }

    /// Initial proportions from each read's argmax emission likelihood
    fn initial_pi(&self) -> Vec<EMProb> {
        let counts: Vec<usize> = (0..self.n_reads())
            .into_par_iter()
            .fold(
                || vec![0usize; self.n_refs()],
                |mut acc, r| {
                    let (s, e) = (self.row_ptr[r], self.row_ptr[r + 1]);
                    if e > s {
                        let mut best = s;
                        for k in (s + 1)..e {
                            if EMProb::total_cmp(&self.lik[k], &self.lik[best]).is_gt() {
                                best = k;
                            }
                        }
                        acc[self.col[best] as usize] += 1;
                    }
                    acc
                },
            )
            .reduce(
                || vec![0usize; self.n_refs()],
                |mut a, b| {
                    for j in 0..a.len() {
                        a[j] += b[j];
                    }
                    a
                },
            );
        let total: usize = counts.iter().sum();
        counts
            .iter()
            .map(|&c| c as EMProb / total as EMProb)
            .collect()
    }

    /// E-step
    fn e_step(&self, pi: &[EMProb]) -> (EMProb, Vec<EMProb>) {
        (0..self.n_reads())
            .into_par_iter()
            .fold(
                || (0.0_f64, vec![0.0_f64; self.n_refs()]),
                |(mut ll, mut ej), r| {
                    let (s, e) = (self.row_ptr[r], self.row_ptr[r + 1]);
                    let mut denom = 0.0;
                    for k in s..e {
                        denom += self.lik[k] * pi[self.col[k] as usize];
                    }
                    if denom != 0.0 {
                        ll += denom.ln();
                        for k in s..e {
                            let w = self.lik[k] * pi[self.col[k] as usize] / denom;
                            if w.is_finite() {
                                ej[self.col[k] as usize] += w;
                            }
                        }
                    }
                    (ll, ej)
                },
            )
            .reduce(
                || (0.0_f64, vec![0.0_f64; self.n_refs()]),
                |(la, mut ea), (lb, eb)| {
                    for j in 0..ea.len() {
                        ea[j] += eb[j];
                    }
                    (la + lb, ea)
                },
            )
    }

    /// The posterior weight of every stored entry under proportions pi
    fn finalize(&self, pi: &[EMProb]) -> (HashMap<ReadIdx, SeqId>, Vec<EMProb>) {
        let mut weights = vec![0.0; self.lik.len()];
        let mut results: HashMap<ReadIdx, SeqId> = HashMap::with_capacity(self.n_reads());
        for r in 0..self.n_reads() {
            let (s, e) = (self.row_ptr[r], self.row_ptr[r + 1]);
            let mut denom = 0.0;
            for k in s..e {
                denom += self.lik[k] * pi[self.col[k] as usize];
            }
            if denom == 0.0 {
                continue;
            }
            let mut best_k = s;
            for k in s..e {
                let post = self.lik[k] * pi[self.col[k] as usize] / denom;
                weights[k] = if post.is_finite() { post } else { 0.0 };
                let cur = self.lik[k] * pi[self.col[k] as usize];
                let best = self.lik[best_k] * pi[self.col[best_k] as usize];
                if EMProb::total_cmp(&cur, &best).is_gt() {
                    best_k = k;
                }
            }
            results.insert(self.reads[r], self.refs[self.col[best_k] as usize]);
        }
        (results, weights)
    }
}

/// Posterior weights of the final EM proportions
pub struct Posteriors<'m> {
    pub(crate) matrix: &'m CsrLikelihood,
    pub(crate) weights: Vec<EMProb>,
}

impl Posteriors<'_> {
    pub fn get(&self, read_idx: ReadIdx, ref_idx: SeqId) -> EMProb {
        self.matrix
            .entry(read_idx, ref_idx)
            .map_or(0.0, |k| self.weights[k])
    }
}

/// Penalty weight used by the "unpenalized" EM path's closed-form M-step.
pub(crate) const UNPENALIZED_RHO: EMProb = 0.0;
/// Proportion floor used by the "unpenalized" EM path's closed-form M-step.
pub(crate) const UNPENALIZED_OMEGA: EMProb = 1e-20;

/// Run the unpenalized EM algorithm to estimate reference proportions.
pub fn get_proportions_par_sparse(
    csr: &CsrLikelihood,
    num_iter: usize,
    progress: Option<Arc<QueryProgress>>,
) -> (
    HashMap<ReadIdx, SeqId>,
    HashMap<SeqId, EMProb>,
    Vec<EMProb>,
    Vec<EMProb>,
) {
    if csr.n_refs() == 0 {
        return (HashMap::new(), HashMap::new(), Vec::new(), Vec::new());
    }
    let num_reads = csr.n_reads();
    let n_refs = csr.n_refs();

    let mut pi = csr.initial_pi();

    let pb = ProgressBar::with_draw_target(Some(num_iter as u64), ProgressDrawTarget::stderr());
    pb.set_style(ProgressStyle::with_template("Running EM: {spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent}% ({eta}) ({msg})").unwrap());

    let mut data_likelihoods: Vec<EMProb> = Vec::new();
    let (mut prev_data_loglikelihood, _) = csr.e_step(&pi);
    data_likelihoods.push(prev_data_loglikelihood);

    for i in 0..num_iter {
        let (data_loglikelihood, ej) = csr.e_step(&pi);

        // M-step: unpenalized closed-form update.
        let lambda_init = ej
            .iter()
            .map(|x| x - UNPENALIZED_RHO)
            .max_by(|f1, f2| EMProb::total_cmp(f1, f2))
            .unwrap();
        let lambda = _update_lambda(
            UNPENALIZED_RHO,
            UNPENALIZED_OMEGA,
            &ej,
            lambda_init,
            num_iter,
        );
        pi = (0..n_refs)
            .map(|j| {
                let tmp_pi = _update_pi(UNPENALIZED_RHO, UNPENALIZED_OMEGA, ej[j], lambda);
                if tmp_pi > 0.0 {
                    tmp_pi
                } else {
                    0.0
                }
            })
            .collect();

        let data_loglikelihood_diff = data_loglikelihood - prev_data_loglikelihood;
        prev_data_loglikelihood = data_loglikelihood;
        data_likelihoods.push(data_loglikelihood);

        pb.set_message(format!("{data_loglikelihood_diff:.3e}"));

        if data_loglikelihood_diff > 0.0 && data_loglikelihood_diff.abs() <= 1e-6 {
            break;
        }

        pb.inc(1);
        if let Some(ref p) = progress {
            p.em_iter_done.store((i + 1) as u64, Ordering::Relaxed);
        }
    }
    pb.finish_with_message(format!("Final data LL: {prev_data_loglikelihood}"));

    // Posteriors are built from the proportions produced by the final M-step
    let (results, w) = csr.finalize(&pi);

    let props: HashMap<SeqId, EMProb> = (0..n_refs)
        .filter(|&j| pi[j] * (num_reads as EMProb) > 1.0)
        .map(|j| (csr.refs[j], pi[j]))
        .collect();

    (results, props, w, data_likelihoods)
}

/// Re-estimate proportions over the classified reads alone, after reclassification.
pub fn refit_proportions_on_classified(
    csr: &CsrLikelihood,
    classified_reads: &HashSet<ReadIdx>,
    prev_props: &HashMap<SeqId, EMProb>,
    rho: EMProb,
    omega: EMProb,
    num_iter: usize,
) -> HashMap<SeqId, EMProb> {
    if classified_reads.is_empty() || csr.n_refs() == 0 {
        return HashMap::new();
    }
    let csr = csr.retain_reads(classified_reads);
    if csr.n_reads() == 0 {
        return HashMap::new();
    }

    let pi: Vec<EMProb> = csr
        .refs
        .iter()
        .map(|r| *prev_props.get(r).unwrap_or(&0.0))
        .collect();

    let (_, ej) = csr.e_step(&pi);

    let lambda_init = ej
        .iter()
        .map(|x| x - rho)
        .max_by(EMProb::total_cmp)
        .unwrap();
    let lambda = _update_lambda(rho, omega, &ej, lambda_init, num_iter);

    (0..csr.n_refs())
        .filter_map(|j| {
            let p = _update_pi(rho, omega, ej[j], lambda);
            if p > 0.0 {
                Some((csr.refs[j], p))
            } else {
                None
            }
        })
        .collect()
}

/// M-step update for a single reference proportion under the L1-regularized objective.
fn _update_pi(rho: EMProb, omega: EMProb, ej: EMProb, lambda: EMProb) -> EMProb {
    let phi = _compute_phi(lambda, omega, rho, ej);
    (-phi + (phi * phi + 4.0 * lambda * ej * omega).sqrt()) / (2.0 * lambda)
}

/// Find the Lagrange multiplier λ that enforces the simplex constraint $\sum\pi_j=1$.
fn _update_lambda(
    rho: EMProb,
    omega: EMProb,
    ejs: &[EMProb],
    lambda_init: EMProb,
    iterations: usize,
) -> EMProb {
    let mut lambda = lambda_init;
    for _ in 0..iterations {
        lambda -=
            (_compute_f(lambda, omega, rho, ejs)) / (_compute_deriv_f(lambda, omega, rho, ejs));
        if _compute_f(lambda, omega, rho, ejs) == 0.0 {
            break;
        }
    }
    lambda
}

/// Evaluate the constraint function f(λ) = Σ_j π_j(λ) − 1 used by Newton-Raphson.
fn _compute_f(lambda: EMProb, omega: EMProb, rho: EMProb, ejs: &[EMProb]) -> EMProb {
    ejs.iter()
        .map(|&ej| {
            let phi = _compute_phi(lambda, omega, rho, ej);
            -phi + (phi * phi + 4.0 * lambda * ej * omega).sqrt()
        })
        .sum::<EMProb>()
        - 2.0 * lambda
}

/// Evaluate the derivative f′(λ) = Σ_j ∂π_j/∂λ used by Newton-Raphson.
fn _compute_deriv_f(lambda: EMProb, omega: EMProb, rho: EMProb, ejs: &[EMProb]) -> EMProb {
    ejs.iter()
        .map(|&ej| {
            let phi = _compute_phi(lambda, omega, rho, ej);
            -omega
                + 0.5
                    * (1.0 / (phi * phi + 4.0 * lambda * ej * omega).sqrt())
                    * (2.0 * phi * omega + 4.0 * ej * omega)
        })
        .sum::<EMProb>()
        - 2.0
}

/// Compute the auxiliary scalar φ = λω + ρ − e_j used in the closed-form π update.
fn _compute_phi(lambda: EMProb, omega: EMProb, rho: EMProb, ej: EMProb) -> EMProb {
    return lambda * omega + rho - ej;
}

/// Run the L1-penalized EM algorithm to estimate reference proportions.
pub fn get_proportions_par_sparse_l1_reg(
    csr: &CsrLikelihood,
    num_iter: usize,
    rho: EMProb,
    omega: EMProb,
    em_threshold: EMProb,
    progress: Option<Arc<QueryProgress>>,
) -> (
    HashMap<ReadIdx, SeqId>,
    HashMap<SeqId, EMProb>,
    Vec<EMProb>,
    Vec<EMProb>,
) {
    if csr.n_refs() == 0 {
        return (HashMap::new(), HashMap::new(), Vec::new(), Vec::new());
    }
    let num_reads = csr.n_reads();
    let n_refs = csr.n_refs();

    let mut pi = csr.initial_pi();

    if let Some(ref p) = progress {
        p.phase.store(2, Ordering::Relaxed);
        p.em_iter_total.store(num_iter as u64, Ordering::Relaxed);
        p.em_iter_done.store(0, Ordering::Relaxed);
    }

    let pb = ProgressBar::with_draw_target(Some(num_iter as u64), ProgressDrawTarget::stderr());
    pb.set_style(ProgressStyle::with_template("Running EM: {spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent}% ({eta}) ({msg})").unwrap());

    let mut data_likelihoods: Vec<EMProb> = Vec::new();
    let (mut prev_data_loglikelihood, _) = csr.e_step(&pi);
    data_likelihoods.push(prev_data_loglikelihood);

    for i in 0..num_iter {
        let (data_loglikelihood, ej) = csr.e_step(&pi);

        let data_loglikelihood_diff = data_loglikelihood - prev_data_loglikelihood;
        prev_data_loglikelihood = data_loglikelihood;
        data_likelihoods.push(data_loglikelihood);

        // M-step: L1-penalized closed-form update solved via Newton-Raphson.
        let lambda_init = ej
            .iter()
            .map(|x| x - rho)
            .max_by(|f1, f2| EMProb::total_cmp(f1, f2))
            .unwrap();
        let lambda = _update_lambda(rho, omega, &ej, lambda_init, num_iter);
        pi = (0..n_refs)
            .map(|j| {
                let tmp_pi = _update_pi(rho, omega, ej[j], lambda);
                if tmp_pi > 0.0 {
                    tmp_pi
                } else {
                    0.0
                }
            })
            .collect();

        pb.set_message(format!("{data_loglikelihood_diff:.3e}"));

        if i > 20 && data_loglikelihood_diff > 0.0 && data_loglikelihood_diff <= em_threshold {
            break;
        }

        pb.inc(1);
        if let Some(ref p) = progress {
            p.em_iter_done.store((i + 1) as u64, Ordering::Relaxed);
        }
    }
    pb.finish_with_message(format!("Final data LL: {prev_data_loglikelihood}"));

    // Posteriors are built from the proportions produced by the final M-step
    let (results, w) = csr.finalize(&pi);

    let props: HashMap<SeqId, EMProb> = (0..n_refs)
        .filter(|&j| pi[j] * (num_reads as EMProb) > 1.0)
        .map(|j| (csr.refs[j], pi[j]))
        .collect();

    (results, props, w, data_likelihoods)
}
