use crate::utils::{
    error_prob, EMProb, MatchLikelihoods, QueryProgress, Read, ReadAlignment, ReadAlignments,
    ReadID, ReadPair,
};
use crate::{RefIndex, SeqId};
use anyhow::{Context, Result};
use bio::stats::{LogProb, Prob};
use haystackfm::alphabet::{self, encode_byte, iupac_bases, ALPHABET_SIZE, N, T};
use haystackfm::BidirInterval;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::sync::atomic::Ordering;
use std::sync::{Arc, LazyLock};

/// a read `N`, which stands for all four bases, scores this against every reference symbol
const LN_4: f64 = 2.0 * std::f64::consts::LN_2;

/// The pairing is an exact match if the read base is not `N` and the reference symbol is that base. Only these lie inside an SMEM.
#[inline]
fn exact(p: u8, b: u8) -> bool {
    p != N && p == b
}

/// The non-sentinel children of a cursor as `(text code, child interval)`. Slot 0
/// is the sentinel child, the occurrences that sit at a reference boundary; it is
/// never followed because a placement must lie entirely inside one reference.
#[inline]
fn kids(
    children: [Option<BidirInterval>; ALPHABET_SIZE],
) -> impl Iterator<Item = (u8, BidirInterval)> {
    children
        .into_iter()
        .enumerate()
        .skip(1)
        .filter_map(|(b, c)| c.map(|c| (b as u8, c)))
}

#[cfg(test)]
pub(crate) mod counters {
    use std::sync::atomic::AtomicUsize;
    pub static EXTENSIONS: AtomicUsize = AtomicUsize::new(0);
    pub static SMEMS: AtomicUsize = AtomicUsize::new(0);
    pub static WALK_NODES: AtomicUsize = AtomicUsize::new(0);
}
#[cfg(test)]
macro_rules! count {
    ($c:ident) => {
        counters::$c.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    };
}
#[cfg(not(test))]
macro_rules! count {
    ($c:ident) => {};
}

/// Log-likelihood terms for one Phred+33 quality byte.
#[derive(Clone, Copy)]
struct PhredRow {
    eps: f64,
    match_ll: f64,
    mismatch_ll: f64,
    /// `set_ll[j]` is the score of a read base that is one of the `j` bases of the
    /// reference symbol: `ln((1 - eps) + (j - 1) eps / 3) - ln j`.
    set_ll: [f64; 5],
}

static PHRED: LazyLock<[PhredRow; 256]> = LazyLock::new(|| {
    let mut table = [PhredRow {
        eps: 0.0,
        match_ll: 0.0,
        mismatch_ll: 0.0,
        set_ll: [0.0; 5],
    }; 256];
    for (q, row) in table.iter_mut().enumerate() {
        let e = *error_prob(q as u8);
        let mut set_ll = [0.0; 5];
        for (j, slot) in set_ll.iter_mut().enumerate().skip(1) {
            let j = j as f64;
            *slot = ((1.0 - e) + (j - 1.0) * e / 3.0).ln() - j.ln();
        }
        *row = PhredRow {
            eps: e,
            match_ll: (1.0 - e).ln(),
            mismatch_ll: (e / 3.0).ln(),
            set_ll,
        };
    }
    table
});

/// Per-position log-likelihood terms of one read, indexed by read position.
struct ScoreTables<'q> {
    quals: &'q [u8],
    /// `prefix_match[t]` = sum of the match scores of positions `..t`.
    prefix_match: Vec<f64>,
}

impl<'q> ScoreTables<'q> {
    fn new(quals: &'q [u8]) -> Self {
        let phred = &*PHRED;
        let mut prefix_match = Vec::with_capacity(quals.len() + 1);
        let mut acc = 0.0;
        prefix_match.push(acc);
        for &q in quals {
            acc += phred[q as usize].match_ll;
            prefix_match.push(acc);
        }
        Self {
            quals,
            prefix_match,
        }
    }

    #[inline]
    fn row(&self, t: usize) -> &'static PhredRow {
        &PHRED[self.quals[t] as usize]
    }

    /// Score of read symbol `p` (at position `t`) against reference symbol `b`.
    #[inline]
    fn score(&self, t: usize, p: u8, b: u8) -> f64 {
        if p == N {
            // All four bases against any reference base set
            -LN_4
        } else if p > T {
            // An ambiguity code in the read
            self.mixture(t, p, b)
        } else {
            let row = self.row(t);
            if b <= T {
                if p == b {
                    row.match_ll
                } else {
                    row.mismatch_ll
                }
            } else {
                let set = iupac_bases(b);
                if set.contains(&p) {
                    row.set_ll[set.len()]
                } else {
                    row.mismatch_ll
                }
            }
        }
    }

    /// The general symbol-pair mixture, for read symbols that are not a base or `N`.
    #[cold]
    fn mixture(&self, t: usize, p: u8, b: u8) -> f64 {
        let (read_set, ref_set) = (iupac_bases(p), iupac_bases(b));
        let n = (read_set.len() * ref_set.len()) as f64;
        if n == 0.0 {
            return 0.0;
        }
        let m = read_set.iter().filter(|x| ref_set.contains(x)).count() as f64;
        let e = self.row(t).eps;
        ((m / n) * (1.0 - e) + ((n - m) / n) * (e / 3.0)).ln()
    }

    /// Score of the exact run `start..end`.
    #[inline]
    fn run(&self, start: usize, end: usize) -> f64 {
        self.prefix_match[end] - self.prefix_match[start]
    }
}

/// Find the SMEM seeds of one read; returns the (SMEM start position, SMEM interval)
pub fn query_smems(
    fmidx: &RefIndex,
    seq: &[u8],
    mem_seed_length: usize,
    complement: bool,
) -> Result<Vec<(usize, BidirInterval)>> {
    let p = encode_read(seq, complement);
    Ok(smem_intervals(fmidx, &p, mem_seed_length))
}

/// Score every read-ref match for a single read as sum of all qualifying alignment
/// Returns log-likelihood of every reported match.
pub fn query_read(
    fmidx: &RefIndex,
    seq: &[u8],
    qual: &[u8],
    mem_seed_length: usize,
    complement: bool,
) -> Result<MatchLikelihoods> {
    if seq.len() != qual.len() {
        anyhow::bail!(
            "read has {} bases but {} quality scores",
            seq.len(),
            qual.len()
        );
    }
    let smems = query_smems(fmidx, seq, mem_seed_length, complement)?;
    let p = encode_read(seq, complement);
    let quals: Vec<u8> = match complement {
        true => qual.iter().rev().copied().collect(),
        false => qual.to_vec(),
    };

    let mut hits = extend_smems(fmidx, &p, &quals, &smems);
    hits.values_mut()
        .for_each(|by_pos| by_pos.retain(|_, ll| ll.exp() != 0.0));
    hits.retain(|_, by_pos| !by_pos.is_empty());

    Ok(MatchLikelihoods(hits))
}

/// Align every read pair against the FM-index in parallel and collect the per-read
/// likelihoods the EM runs on.
///
/// Returns the surviving alignments as one sparse read × reference table (the input of
/// [`crate::em::CsrLikelihood::build`]) and the set of every reference that received at
/// least one alignment.
pub fn query_fastq<'a>(
    fmidx: &RefIndex,
    read_pairs: &'a [ReadPair],
    mem_seed_length: usize,
    eps_1: LogProb,
    eps_2: LogProb,
    progress: Option<Arc<QueryProgress>>,
) -> Result<(ReadAlignments<'a>, HashSet<SeqId>)> {
    let pb =
        ProgressBar::with_draw_target(Some(read_pairs.len() as u64), ProgressDrawTarget::stderr());
    pb.set_style(ProgressStyle::with_template("Finding pairwise alignments: {spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent}% ({eta})").unwrap());

    if let Some(ref p) = progress {
        p.phase.store(1, Ordering::Relaxed);
        p.reads_total
            .store(read_pairs.len() as u64, Ordering::Relaxed);
    }

    let rows: Vec<(&ReadID, Vec<(SeqId, ReadAlignment)>)> = read_pairs
        .par_iter()
        .map(|pair| -> Result<_> {
            let row = align_pair(fmidx, pair, mem_seed_length, eps_1, eps_2)?;
            pb.inc(1);
            if let Some(ref p) = progress {
                p.reads_done.fetch_add(1, Ordering::Relaxed);
            }
            Ok((&pair.read_id, row))
        })
        .collect::<Result<Vec<_>>>()?;

    pb.finish_with_message("");

    let all_ref_ids: HashSet<SeqId> = rows
        .iter()
        .flat_map(|(_, row)| row.iter().map(|(ref_idx, _)| *ref_idx))
        .collect();
    Ok((ReadAlignments::from_rows(rows), all_ref_ids))
}

/// Score every qualifying match of a read using quality scores at seed length `k`
///
/// Returns log-likelihood of match
pub fn score_read(
    idx: &RefIndex,
    p: &[u8],
    quals: &[u8],
    k: usize,
) -> HashMap<SeqId, HashMap<usize, LogProb>> {
    let smems = smem_intervals(idx, p, k);
    extend_smems(idx, p, quals, &smems)
}

/// The read as alignment codes
fn encode_read(seq: &[u8], complement: bool) -> Vec<u8> {
    let encode = |b: u8| {
        encode_byte(b)
            .filter(|&c| c <= alphabet::N)
            .unwrap_or(alphabet::N)
    };
    match complement {
        true => bio::alphabets::dna::revcomp(seq)
            .into_iter()
            .map(encode)
            .collect(),
        false => seq.iter().map(|&b| encode(b)).collect(),
    }
}

/// Enumerates all SMEMs of length at least `k` with intervals
fn smem_intervals(idx: &RefIndex, p: &[u8], k: usize) -> Vec<(usize, BidirInterval)> {
    let n = p.len();
    if n == 0 {
        return Vec::new();
    }
    smems(idx, p, k.clamp(1, n))
}

/// Extends each SMEM over the whole read with score
fn extend_smems(
    idx: &RefIndex,
    p: &[u8],
    quals: &[u8],
    smems: &[(usize, BidirInterval)],
) -> HashMap<SeqId, HashMap<usize, LogProb>> {
    let n = p.len();
    if n == 0 {
        return HashMap::new();
    }
    let tables = ScoreTables::new(quals);
    let mut end_at = vec![None; n];
    for &(start, cursor) in smems {
        end_at[start] = Some(start + cursor.len as usize);
    }
    let mut walk = Walk {
        idx,
        p,
        tables: &tables,
        n,
        end_at,
        out: HashMap::new(),
    };
    for &(start, cursor) in smems {
        walk.seed(start, cursor);
    }
    walk.out
}

/// The cursor of the exact k-mer `p[at..at + k]` read straight from the index's lookup
/// tables (`k = idx.lookup_depth()`), standing in for `k` exact extensions from the full
/// interval. `None` when the index has no tables, the window runs past the read or holds
/// an `N`, or the k-mer does not occur; the caller then extends one symbol at a time.
#[inline]
fn lookup(idx: &RefIndex, p: &[u8], at: usize) -> Option<BidirInterval> {
    let k = idx.lookup_depth() as usize;
    if k == 0 {
        return None;
    }
    let kmer = p.get(at..at + k)?;
    if kmer.contains(&N) {
        return None;
    }
    idx.lookup_interval(kmer)
}

/// Finds all SMEM using pivot jumping of BWA-MEM
/// Runs in `O(n + Σ overlaps)`. A run that starts from the full interval first tries the
/// k-mer lookup tables, which give the same cursor as `k` successful exact extensions.
fn smems(idx: &RefIndex, p: &[u8], k: usize) -> Vec<(usize, BidirInterval)> {
    let n = p.len();
    let full = idx.full_interval();
    let depth = idx.lookup_depth() as usize;
    let mut out = Vec::new();
    let (mut i, mut e, mut x) = (0usize, 0usize, full);
    loop {
        if e == i {
            if let Some(y) = lookup(idx, p, e) {
                x = y;
                e += depth;
            }
        }
        while e < n && p[e] != N {
            count!(EXTENSIONS);
            match idx.extend_right(x, p[e]) {
                Some(y) => {
                    x = y;
                    e += 1;
                }
                None => break,
            }
        }
        if e - i >= k {
            count!(SMEMS);
            debug_assert_eq!(x.len as usize, e - i, "cursor length is the run length");
            out.push((i, x));
        }
        if e == n {
            break;
        }
        // Find the next SMEM start and its cursor.
        let restart = p[e] == N;
        let mut next = None;
        if !restart {
            // `y` spells `P[j..e+1)`; walk left while it keeps occurring. When the k-mer
            // ending at `e + 1` occurs, so does each of its suffixes, and the walk would
            // reach its start: take it from the tables and walk on from there.
            let jumped = match (e + 1).checked_sub(depth) {
                Some(j) if depth > 0 && j > i => lookup(idx, p, j).map(|y| (j, y)),
                _ => None,
            };
            let seed = jumped.or_else(|| {
                count!(EXTENSIONS);
                idx.extend_right(full, p[e]).map(|y| (e, y))
            });
            if let Some((mut j, mut y)) = seed {
                while j > i + 1 && p[j - 1] != N {
                    count!(EXTENSIONS);
                    match idx.extend_left(y, p[j - 1]) {
                        Some(z) => {
                            y = z;
                            j -= 1;
                        }
                        None => break,
                    }
                }
                next = Some((j, y));
            }
        }
        match next {
            Some((j, y)) => {
                i = j;
                x = y;
                e += 1;
            }
            None => {
                i = e + 1;
                e = i;
                x = full;
            }
        }
        if i + k > n {
            break;
        }
    }
    out
}

/// Depth-first cursor walk that extends each SMEM over the rest of the read.
struct Walk<'a> {
    idx: &'a RefIndex,
    p: &'a [u8],
    tables: &'a ScoreTables<'a>,
    n: usize,
    /// `end_at[s]` is the end of the reported SMEM starting at `s`, if any. A left
    /// branch that closes exactly such a run is owned by that SMEM.
    end_at: Vec<Option<usize>>,
    out: HashMap<SeqId, HashMap<usize, LogProb>>,
}

impl Walk<'_> {
    /// Split the occurrences of the SMEM starting at `s` (its cursor `cursor`) by
    /// their right flank and then by their left flank, scoring both flanking
    /// positions, and walk each branch.
    fn seed(&mut self, s: usize, cursor: BidirInterval) {
        let e = s + cursor.len as usize;
        let ll = self.tables.run(s, e);
        if e == self.n {
            self.after_right_flank(cursor, s, self.n, ll);
            return;
        }
        for (c, f) in kids(self.idx.children_right(&cursor)) {
            debug_assert!(!exact(self.p[e], c), "the SMEM's run would continue");
            let g = ll + self.tables.score(e, self.p[e], c);
            self.after_right_flank(f, s, e + 1, g);
        }
    }

    /// `f` spells `P[s..e)` with its right flank consumed; `right_start` is where the
    /// rightward walk resumes. Fork the left flank (a no-op filter for a true SMEM:
    /// no occurrence is preceded by `P[s - 1]`) and walk leftward.
    fn after_right_flank(&mut self, f: BidirInterval, s: usize, right_start: usize, g: f64) {
        if s == 0 {
            return self.walk_right(f, right_start, g);
        }
        for (a, c) in kids(self.idx.children_left(&f)) {
            debug_assert!(!exact(self.p[s - 1], a), "the SMEM was not left-maximal");
            let g = g + self.tables.score(s - 1, self.p[s - 1], a);
            self.walk_left(c, s as isize - 2, g, s - 1, right_start);
        }
    }

    /// Score positions `t, t-1, …, 0` leftward, then continue rightward from
    /// `right_start`. `run_end` is the exclusive end of the exact run currently being
    /// built (0 when there is none); a branch whose closed run is a reported SMEM is
    /// dropped, because that SMEM owns the diagonal.
    fn walk_left(
        &mut self,
        cur: BidirInterval,
        t: isize,
        g: f64,
        run_end: usize,
        right_start: usize,
    ) {
        if t < 0 {
            if run_end > 0 && self.end_at[0] == Some(run_end) {
                return;
            }
            return self.walk_right(cur, right_start, g);
        }
        let t = t as usize;
        let p = self.p[t];
        count!(WALK_NODES);
        for (b, c) in kids(self.idx.children_left(&cur)) {
            let g = g + self.tables.score(t, p, b);
            if exact(p, b) {
                // The exact run continues; its end is still `run_end`.
                self.walk_left(c, t as isize - 1, g, run_end, right_start);
            } else if self.end_at[t + 1] != Some(run_end) {
                // The run `[t + 1, run_end)` closes here and is not a reported SMEM.
                self.walk_left(c, t as isize - 1, g, t, right_start);
            }
        }
    }

    /// Score positions `t, t+1, …, n-1` rightward, then emit. Nothing is pruned: the
    /// leftmost SMEM run owns the diagonal.
    fn walk_right(&mut self, cur: BidirInterval, t: usize, g: f64) {
        if t == self.n {
            return self.emit(&cur, g);
        }
        let p = self.p[t];
        count!(WALK_NODES);
        for (b, c) in kids(self.idx.children_right(&cur)) {
            self.walk_right(c, t + 1, g + self.tables.score(t, p, b));
        }
    }

    /// `cur` spans the whole read, so every located offset is a diagonal and every
    /// window in the interval has score `g`.
    fn emit(&mut self, cur: &BidirInterval, g: f64) {
        for (id, d) in self.idx.locate_interval(cur) {
            let prev = self
                .out
                .entry(id)
                .or_default()
                .insert(d as usize, LogProb(g));
            debug_assert!(prev.is_none(), "placement ({id:?}, {d}) emitted twice");
        }
    }
}

/// Compute the combined log-probability for a read pair aligning to the same reference.
pub fn merge_read_pairs(
    forward: &HashMap<usize, LogProb>,
    reverse: &HashMap<usize, LogProb>,
    read_len: usize,
) -> (LogProb, LogProb, (usize, usize)) {
    let mut match_likelihood: EMProb = 0.0;
    let mut best_alignment_likelihood: LogProb = LogProb(f64::NEG_INFINITY);
    let mut best_alignment_positions: (usize, usize) = (0, 0);
    for (r1_start, r1_log_prob) in forward.iter() {
        let r1_end = r1_start + read_len;
        for (r2_rc_start, r2_log_prob) in reverse.iter() {
            let r2_rc_end = r2_rc_start + read_len;
            if r2_rc_end > *r1_start || r1_end > *r2_rc_start {
                let align_ll = r1_log_prob + r2_log_prob;
                if align_ll > best_alignment_likelihood {
                    best_alignment_likelihood = align_ll;
                    best_alignment_positions = (*r1_start, *r2_rc_start)
                }
                match_likelihood += (align_ll).exp();
            }
        }
    }
    (
        LogProb(match_likelihood.ln()),
        best_alignment_likelihood,
        best_alignment_positions,
    )
}

/// One row of [`query_fastq`]
fn align_pair(
    fmidx: &RefIndex,
    pair: &ReadPair,
    mem_seed_length: usize,
    eps_1: LogProb,
    eps_2: LogProb,
) -> Result<Vec<(SeqId, ReadAlignment)>> {
    let half_log_prob = LogProb::from(Prob(0.5));
    let id = pair.read_id.as_str();
    let query = |read: &Read, complement: bool, what: &str| {
        query_read(fmidx, read.seq(), read.qual(), mem_seed_length, complement)
            .with_context(|| format!("failed to align read '{id}' ({what})"))
    };
    let r1 = query(&pair.r1, false, "R1")?;
    let r1_rc = query(&pair.r1, true, "R1, reverse-complement")?;
    let r2 = query(&pair.r2, false, "R2")?;
    let r2_rc = query(&pair.r2, true, "R2, reverse-complement")?;

    let mut by_ref: HashMap<SeqId, ReadAlignment> = HashMap::new();
    // FR: R1 forward with R2 reverse-complemented.
    for (ref_idx, fwd) in r1.iter() {
        let Some(rev) = r2_rc.get(ref_idx) else {
            continue;
        };
        let (ll, best, pos) = merge_mates(fwd, rev);
        by_ref.insert(
            *ref_idx,
            ReadAlignment {
                full_match_ll: half_log_prob + ll,
                best_match_ll: best,
                pos,
                complement: false,
            },
        );
    }
    // RF: R2 forward with R1 reverse-complemented.
    for (ref_idx, fwd) in r2.iter() {
        let Some(rev) = r1_rc.get(ref_idx) else {
            continue;
        };
        let (ll, best, pos) = merge_mates(fwd, rev);
        let match_prob = half_log_prob + ll;
        match by_ref.get_mut(ref_idx) {
            Some(v) => {
                v.increment_full_match_ll(match_prob);
                if best > v.get_best_match_ll() {
                    v.update_best_match_ll(best);
                    v.update_pos(pos);
                    v.update_compl(true);
                }
            }
            None => {
                by_ref.insert(
                    *ref_idx,
                    ReadAlignment {
                        full_match_ll: match_prob,
                        best_match_ll: best,
                        pos,
                        complement: false,
                    },
                );
            }
        }
    }

    if let Some(max_likelihood) = by_ref
        .values()
        .max_by(|a, b| a.full_match_ll.total_cmp(&b.full_match_ll))
        .map(|v| v.best_match_ll)
    {
        by_ref
            .retain(|_, v| !(v.full_match_ll < eps_1 || v.full_match_ll - max_likelihood <= eps_2));
    }
    let mut row: Vec<(SeqId, ReadAlignment)> = by_ref.into_iter().collect();
    row.sort_unstable_by_key(|(ref_idx, _)| ref_idx.0);
    Ok(row)
}

/// Combine one mate's placements on a reference with the other mate's
fn merge_mates(
    forward: &HashMap<usize, LogProb>,
    reverse: &HashMap<usize, LogProb>,
) -> (LogProb, LogProb, (usize, usize)) {
    let (f_sum, f_best, f_pos) = summarize_placements(forward);
    let (r_sum, r_best, r_pos) = summarize_placements(reverse);
    (
        LogProb((f_sum * r_sum).ln()),
        f_best + r_best,
        (f_pos, r_pos),
    )
}

/// Linear-space sum of a mate's placement likelihoods, plus its best placement.
fn summarize_placements(placements: &HashMap<usize, LogProb>) -> (f64, LogProb, usize) {
    let mut sum = 0.0;
    let mut best = LogProb(f64::NEG_INFINITY);
    let mut pos = 0;
    for (&p, &ll) in placements {
        sum += ll.exp();
        if ll > best {
            best = ll;
            pos = p;
        }
    }
    (sum, best, pos)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::compute_match_log_prob;
    use crate::{build_index_from_bytes, load_index};
    use haystackfm::alphabet::{encode_char, A, C};

    const Q40: u8 = b'I';

    /// Shorthand for [`score_read`].
    fn sr(
        idx: &RefIndex,
        p: &[u8],
        quals: &[u8],
        k: usize,
    ) -> HashMap<SeqId, HashMap<usize, LogProb>> {
        score_read(idx, p, quals, k)
    }

    fn index_raw(refs: &[(&str, &str)]) -> RefIndex {
        index_raw_depth(refs, 0)
    }

    /// [`index_raw`] with k-mer lookup tables of depth `lookup_depth`.
    fn index_raw_depth(refs: &[(&str, &str)], lookup_depth: u32) -> RefIndex {
        use haystackfm::{DnaSequence, FmIndexConfig};
        let cfg = FmIndexConfig {
            sa_sample_rate: 1,
            use_gpu: false,
            lookup_depth,
            ..Default::default()
        };
        let dna: Vec<DnaSequence> = refs
            .iter()
            .map(|(h, s)| DnaSequence::from_str_with_header(s, h).unwrap())
            .collect();
        RefIndex::build_cpu(&dna, &cfg).unwrap()
    }

    fn index_from(refs: &[(&str, &str)]) -> RefIndex {
        let fasta: Vec<u8> = refs
            .iter()
            .flat_map(|(header, seq)| format!(">{header}\n{seq}\n").into_bytes())
            .collect();
        let (bytes, _) = build_index_from_bytes(&fasta).expect("index build");
        load_index(&bytes, "<test>").expect("index load")
    }

    fn enc(ascii: &str) -> Vec<u8> {
        ascii.chars().map(|c| encode_char(c).unwrap()).collect()
    }

    fn q40(n: usize) -> Vec<u8> {
        vec![Q40; n]
    }

    fn m() -> f64 {
        (1.0 - *error_prob(Q40)).ln()
    }

    fn mm() -> f64 {
        (*error_prob(Q40) / 3.0).ln()
    }

    fn set_m(j: u32) -> f64 {
        let e = *error_prob(Q40);
        let j = f64::from(j);
        ((1.0 - e) + (j - 1.0) * e / 3.0).ln() - j.ln()
    }

    fn nm() -> f64 {
        -LN_4
    }

    fn lcg(state: &mut u64) -> u64 {
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        *state >> 33
    }

    fn random_dna(len: usize, state: &mut u64) -> String {
        const B: [char; 4] = ['A', 'C', 'G', 'T'];
        (0..len).map(|_| B[lcg(state) as usize & 3]).collect()
    }

    fn brute_force_reach(idx: &RefIndex, read: &[u8]) -> Vec<usize> {
        let n = read.len();
        let texts: Vec<&[u8]> = (0..idx.num_sequences())
            .map(|j| idx.sequence(SeqId(j)).unwrap())
            .collect();
        (0..n)
            .map(|i| {
                let mut best = i;
                for text in &texts {
                    for start in 0..text.len() {
                        let mut e = i;
                        while e < n
                            && start + (e - i) < text.len()
                            && exact(read[e], text[start + (e - i)])
                        {
                            e += 1;
                        }
                        best = best.max(e);
                    }
                }
                best
            })
            .collect()
    }

    /// The SMEMs of length `>= k`: `[i, reach(i))` for `i == 0` and wherever `reach`
    /// increases.
    fn brute_force_smems(idx: &RefIndex, read: &[u8], k: usize) -> Vec<(usize, usize)> {
        let n = read.len();
        if n == 0 {
            return Vec::new();
        }
        let k = k.clamp(1, n);
        let reach = brute_force_reach(idx, read);
        (0..n)
            .filter(|&i| i == 0 || reach[i] > reach[i - 1])
            .map(|i| (i, reach[i]))
            .filter(|&(i, e)| e - i >= k)
            .collect()
    }

    /// Every occurrence of every SMEM as `(ref, run position, start, end)`, dropping
    /// those whose read placement would not fit inside the reference.
    fn brute_force_smem_occurrences(
        idx: &RefIndex,
        read: &[u8],
        k: usize,
    ) -> Vec<(u32, usize, usize, usize)> {
        let n = read.len();
        let mut out = Vec::new();
        for (s, e) in brute_force_smems(idx, read, k) {
            for j in 0..idx.num_sequences() {
                let text = idx.sequence(SeqId(j)).unwrap();
                for pos in 0..text.len().saturating_sub(e - s - 1) {
                    if (s..e).all(|t| exact(read[t], text[pos + t - s])) {
                        let Some(d) = pos.checked_sub(s) else {
                            continue;
                        };
                        if d + n <= text.len() {
                            out.push((j, pos, s, e));
                        }
                    }
                }
            }
        }
        out.sort_unstable();
        out
    }

    /// Every qualifying placement and its full-read log-likelihood: the diagonals
    /// carrying an occurrence of an SMEM of length `>= k`.
    fn brute_force(
        idx: &RefIndex,
        read: &[u8],
        quals: &[u8],
        k: usize,
    ) -> HashMap<SeqId, HashMap<usize, LogProb>> {
        let n = read.len();
        let mut out: HashMap<SeqId, HashMap<usize, LogProb>> = HashMap::new();
        for (j, pos, s, _) in brute_force_smem_occurrences(idx, read, k) {
            let id = SeqId(j);
            let text = idx.sequence(id).unwrap();
            let d = pos - s;
            let win = &text[d..d + n];
            out.entry(id)
                .or_default()
                .insert(d, compute_match_log_prob(read, quals, win));
        }
        out
    }

    /// Resolve [`smem_intervals`] to `(ref, run position, start, end)` per occurrence,
    /// checking each run's log-likelihood, and dropping the occurrences whose placement
    /// does not fit (which the sweep keeps but the walk cannot emit).
    fn smem_occurrences(
        idx: &RefIndex,
        read: &[u8],
        quals: &[u8],
        k: usize,
    ) -> Vec<(u32, usize, usize, usize)> {
        let n = read.len();
        let mut out = Vec::new();
        let tables = ScoreTables::new(quals);
        for (start, cursor) in smem_intervals(idx, read, k) {
            let end = start + cursor.len as usize;
            let run = &read[start..end];
            assert_close(
                tables.run(start, end),
                *compute_match_log_prob(run, &quals[start..end], run),
            );
            for (id, pos) in idx.locate_interval(&cursor) {
                let pos = pos as usize;
                let len = idx.sequence(id).unwrap().len();
                let Some(d) = pos.checked_sub(start) else {
                    continue;
                };
                if d + n <= len {
                    out.push((id.0, pos, start, end));
                }
            }
        }
        out.sort_unstable();
        out
    }

    fn assert_smems_match_oracle(idx: &RefIndex, read: &[u8], quals: &[u8], k: usize, ctx: &str) {
        assert_eq!(
            smem_occurrences(idx, read, quals, k),
            brute_force_smem_occurrences(idx, read, k),
            "SMEM occurrences differ: {ctx}"
        );
    }

    fn assert_placements_match_oracle(
        idx: &RefIndex,
        read: &[u8],
        quals: &[u8],
        k: usize,
        ctx: &str,
    ) {
        assert_same(
            &score_read(idx, read, quals, k),
            &brute_force(idx, read, quals, k),
            ctx,
        );
    }

    fn assert_same(
        got: &HashMap<SeqId, HashMap<usize, LogProb>>,
        want: &HashMap<SeqId, HashMap<usize, LogProb>>,
        ctx: &str,
    ) {
        let mut got_keys: Vec<_> = got
            .iter()
            .flat_map(|(id, m)| m.keys().map(move |d| (*id, *d)))
            .collect();
        let mut want_keys: Vec<_> = want
            .iter()
            .flat_map(|(id, m)| m.keys().map(move |d| (*id, *d)))
            .collect();
        got_keys.sort();
        want_keys.sort();
        assert_eq!(got_keys, want_keys, "placement sets differ: {ctx}");
        for (id, d) in got_keys {
            let a = *got[&id][&d];
            let b = *want[&id][&d];
            assert!(
                (a - b).abs() < 1e-9,
                "score differs at ({id:?}, {d}): {a} vs {b}: {ctx}"
            );
        }
    }

    fn single(got: &HashMap<SeqId, HashMap<usize, LogProb>>, id: u32) -> (usize, f64) {
        let m = &got[&SeqId(id)];
        assert_eq!(m.len(), 1, "expected one diagonal on ref {id}, got {m:?}");
        let (d, ll) = m.iter().next().unwrap();
        (*d, **ll)
    }

    fn assert_close(a: f64, b: f64) {
        assert!((a - b).abs() < 1e-9, "{a} vs {b}");
    }

    /// Reference `a` holds the read exactly, so `[0, n)` is the only SMEM and
    /// reference `b`'s six-base run is shadowed: it seeds nothing and `b` receives no
    /// placement. Shortening the read below `b`'s run makes both qualify again.
    /// The walk's own scorer is the same model as [`compute_match_log_prob`], which
    /// `utils.rs` shares with the rest of premise: every non-sentinel symbol pair, at
    /// several qualities, must agree exactly.
    #[test]
    fn score_tables_agree_with_compute_match_log_prob() {
        let quals: Vec<u8> = vec![b'!' + 2, b'5', Q40, b'~'];
        let tables = ScoreTables::new(&quals);
        for (t, &q) in quals.iter().enumerate() {
            for p in 1..ALPHABET_SIZE as u8 {
                for b in 1..ALPHABET_SIZE as u8 {
                    let want = *compute_match_log_prob(&[p], &[q], &[b]);
                    let got = tables.score(t, p, b);
                    assert!(
                        (got - want).abs() < 1e-12,
                        "read {p} vs ref {b} at Q{}: {got} vs {want}",
                        q - 33
                    );
                }
            }
            // The exact-run prefix sums are the plain-match case of the same model.
            assert_close(
                tables.run(t, t + 1),
                *compute_match_log_prob(&[A], &[q], &[A]),
            );
        }
    }

    #[test]
    fn sub_smem_run_is_shadowed() {
        let idx = index_from(&[("a", "ACGTACGTAC"), ("b", "ACGTACTTTT")]);
        let read = enc("ACGTACGTAC");
        let got = sr(&idx, &read, &q40(10), 6);

        assert_eq!(got.len(), 1, "{got:?}");
        let (d_a, ll_a) = single(&got, 0);
        assert_eq!(d_a, 0);
        assert_close(ll_a, 10.0 * m());
        assert_placements_match_oracle(&idx, &read, &q40(10), 6, "shadowed run");

        // `ACGTAC` is its own SMEM: both references qualify, `a` on two diagonals.
        let read = enc("ACGTAC");
        let got = sr(&idx, &read, &q40(6), 6);
        assert_eq!(got.len(), 2, "{got:?}");
        assert_close(*got[&SeqId(0)][&0], 6.0 * m());
        assert_close(*got[&SeqId(0)][&4], 6.0 * m());
        assert_close(*got[&SeqId(1)][&0], 6.0 * m());
        assert_placements_match_oracle(&idx, &read, &q40(6), 6, "shared SMEM");
    }

    #[test]
    fn alignment_with_two_qualifying_runs_scored_once() {
        let mut st = 7u64;
        let r = random_dna(60, &mut st);
        let idx = index_from(&[("r", &r)]);
        let mut read = enc(&r[..30]);
        read[15] = if read[15] == A { C } else { A };
        let got = sr(&idx, &read, &q40(30), 11);

        let (d, ll) = single(&got, 0);
        assert_eq!(d, 0);
        assert_close(ll, 29.0 * m() + mm());
        assert_same(&got, &brute_force(&idx, &read, &q40(30), 11), "two runs");
    }

    /// Every reference holding the window exactly receives its own placement. The
    /// reference that differs inside it does not: both of its runs are sub-matches of
    /// the full-length match in the others, so neither is an SMEM.
    #[test]
    fn shared_window_in_many_references_yields_entry_per_reference() {
        let mut st = 11u64;
        let block = random_dna(40, &mut st);
        let mut block_mut = block.clone();
        block_mut.replace_range(20..21, if &block[20..21] == "A" { "C" } else { "A" });
        let r0 = format!(
            "{}{}{}",
            random_dna(10, &mut st),
            block,
            random_dna(15, &mut st)
        );
        let r1 = format!(
            "{}{}{}",
            random_dna(20, &mut st),
            block,
            random_dna(15, &mut st)
        );
        let r2 = format!(
            "{}{}{}",
            random_dna(30, &mut st),
            block_mut,
            random_dna(15, &mut st)
        );
        let idx = index_from(&[("r0", &r0), ("r1", &r1), ("r2", &r2)]);
        let read = enc(&block);
        let got = sr(&idx, &read, &q40(40), 11);

        let (d0, ll0) = single(&got, 0);
        assert_eq!(d0, 10);
        assert_close(ll0, 40.0 * m());
        assert_eq!(single(&got, 1).0, 20);
        assert!(!got.contains_key(&SeqId(2)), "shadowed reference: {got:?}");
        assert_same(
            &got,
            &brute_force(&idx, &read, &q40(40), 11),
            "shared window",
        );
    }

    #[test]
    fn n_in_read_is_wildcard_scored_quarter() {
        let mut st = 3u64;
        let r = random_dna(60, &mut st);
        let idx = index_from(&[("r", &r)]);
        let mut read = enc(&r[..50]);
        read[15] = N;
        let got = sr(&idx, &read, &q40(50), 11);

        let (d, ll) = single(&got, 0);
        assert_eq!(d, 0);
        assert_close(ll, 49.0 * m() + nm());
    }

    #[test]
    fn n_in_reference_is_wildcard() {
        let mut st = 5u64;
        let mut r = random_dna(60, &mut st);
        r.replace_range(20..21, "N");
        let idx = index_from(&[("r", &r)]);
        let mut read_ascii = r[..50].to_string();
        read_ascii.replace_range(20..21, "A");
        let read = enc(&read_ascii);
        let got = sr(&idx, &read, &q40(50), 11);

        let (d, ll) = single(&got, 0);
        assert_eq!(d, 0);
        assert_close(ll, 49.0 * m() + nm());
        assert_same(&got, &brute_force(&idx, &read, &q40(50), 11), "ref N");
    }

    #[test]
    fn read_n_run_of_k_does_not_qualify_alone() {
        let mut st = 9u64;
        let r = random_dna(60, &mut st);
        let idx = index_from(&[("r", &r)]);

        let read = enc("NNNNNACGT");
        assert!(sr(&idx, &read, &q40(9), 5).is_empty());
        let read = enc("NNNNNNNN");
        assert!(sr(&idx, &read, &q40(8), 3).is_empty());
    }

    #[test]
    fn seed_start_after_read_n_forks_all_symbols() {
        let mut st = 13u64;
        let r = random_dna(60, &mut st);
        let idx = index_from(&[("r", &r)]);
        let mut read = enc(&r[..30]);
        read[5] = if read[5] == A { C } else { A };
        read[10] = N;
        let got = sr(&idx, &read, &q40(30), 10);

        let (d, ll) = single(&got, 0);
        assert_eq!(d, 0);
        assert_close(ll, 28.0 * m() + mm() + nm());
        assert_same(
            &got,
            &brute_force(&idx, &read, &q40(30), 10),
            "N before seed",
        );
    }

    #[test]
    fn empty_read_returns_nothing() {
        let idx = index_from(&[("a", "ACGTACGTAC")]);
        assert!(sr(&idx, &[], &[], 5).is_empty());
    }

    #[test]
    fn read_shorter_than_k_uses_read_length_as_k() {
        let mut st = 17u64;
        let r = random_dna(60, &mut st);
        let idx = index_from(&[("r", &r)]);
        let read = enc(&r[5..13]);
        let got = sr(&idx, &read, &q40(8), 11);
        let (d, ll) = single(&got, 0);
        assert_eq!(d, 5);
        assert_close(ll, 8.0 * m());

        let mut read = read;
        read[3] = if read[3] == A { C } else { A };
        assert!(sr(&idx, &read, &q40(8), 11).is_empty());
    }

    #[test]
    fn read_equal_to_k_requires_full_window() {
        let mut st = 19u64;
        let r = random_dna(60, &mut st);
        let idx = index_from(&[("r", &r)]);
        let read = enc(&r[5..16]);
        let got = sr(&idx, &read, &q40(11), 11);
        let (d, ll) = single(&got, 0);
        assert_eq!(d, 5);
        assert_close(ll, 11.0 * m());

        let mut read = read;
        read[0] = if read[0] == A { C } else { A };
        assert!(sr(&idx, &read, &q40(11), 11).is_empty());
    }

    /// A read `N` at position 0 means the first seed starts at 1 with `P[0] == N`,
    /// so `fork_left_maximal` must fork over all five symbols at `t == 0`.
    #[test]
    fn read_n_at_first_position() {
        let mut st = 23u64;
        let r = random_dna(60, &mut st);
        let idx = index_from(&[("r", &r)]);
        let mut read = enc(&r[..30]);
        read[0] = N;
        let got = sr(&idx, &read, &q40(30), 11);

        let (d, ll) = single(&got, 0);
        assert_eq!(d, 0);
        assert_close(ll, 29.0 * m() + nm());
        assert_same(&got, &brute_force(&idx, &read, &q40(30), 11), "N at 0");
    }

    /// A read `N` at the last position is scored by `walk_right` on its final
    /// step before `t == n` emits.
    #[test]
    fn read_n_at_last_position() {
        let mut st = 29u64;
        let r = random_dna(60, &mut st);
        let idx = index_from(&[("r", &r)]);
        let mut read = enc(&r[..30]);
        read[29] = N;
        let got = sr(&idx, &read, &q40(30), 11);

        let (d, ll) = single(&got, 0);
        assert_eq!(d, 0);
        assert_close(ll, 29.0 * m() + nm());
        assert_same(&got, &brute_force(&idx, &read, &q40(30), 11), "N at n-1");
    }

    /// References shorter than the read never yield a placement, even when a
    /// longer reference in the same index does.
    #[test]
    fn reference_shorter_than_read_is_skipped() {
        let mut st = 31u64;
        let r = random_dna(60, &mut st);
        let short = &r[5..25];
        let idx = index_from(&[("short", short), ("long", &r)]);
        let read = enc(&r[..30]);
        let got = sr(&idx, &read, &q40(30), 11);

        assert!(!got.contains_key(&SeqId(0)), "{got:?}");
        assert_eq!(single(&got, 1).0, 0);
        assert_same(&got, &brute_force(&idx, &read, &q40(30), 11), "short ref");
    }

    /// One random case: 2-3 references sharing a 30 bp block, ~2% reference `N`,
    /// a read sampled from one of them with substitutions and `N`s.
    type Case = (Vec<(String, String)>, Vec<u8>, Vec<u8>, usize);

    fn random_case(st: &mut u64) -> Case {
        random_case_with(st, &['N'])
    }

    /// Like [`random_case`], but every 50th reference symbol on average is drawn
    /// from `wild` instead of always being `N`.
    fn random_case_with(st: &mut u64, wild: &[char]) -> Case {
        let block = random_dna(30, st);
        let n_refs = 2 + (lcg(st) % 2) as usize;
        let refs: Vec<(String, String)> = (0..n_refs)
            .map(|j| {
                let total = 60 + (lcg(st) % 61) as usize;
                let at = (lcg(st) % (total - 30) as u64) as usize;
                let mut s: Vec<char> = format!(
                    "{}{}{}",
                    random_dna(at, st),
                    block,
                    random_dna(total - 30 - at, st)
                )
                .chars()
                .collect();
                for c in s.iter_mut() {
                    let w = lcg(st);
                    if w.is_multiple_of(50) {
                        *c = wild[(w / 50 % wild.len() as u64) as usize];
                    }
                }
                (format!("r{j}"), s.into_iter().collect())
            })
            .collect();
        let src = &refs[(lcg(st) % n_refs as u64) as usize].1;
        let len = 8 + (lcg(st) % 33) as usize;
        let start = (lcg(st) % (src.len() - len + 1) as u64) as usize;
        // Reads carry only ACGTN: `query_read` maps every other symbol to `N`.
        let mut read: Vec<char> = src[start..start + len]
            .chars()
            .map(|c| if "ACGT".contains(c) { c } else { 'N' })
            .collect();
        for _ in 0..(lcg(st) % 4) {
            let i = (lcg(st) % len as u64) as usize;
            read[i] = match read[i] {
                'A' => 'C',
                'C' => 'G',
                'G' => 'T',
                _ => 'A',
            };
        }
        for _ in 0..(lcg(st) % 3) {
            let i = (lcg(st) % len as u64) as usize;
            read[i] = 'N';
        }
        let quals: Vec<u8> = (0..len).map(|_| 33 + 2 + (lcg(st) % 39) as u8).collect();
        let k = [3usize, 5, 8, 11][(lcg(st) % 4) as usize];
        let read: String = read.into_iter().collect();
        (refs, enc(&read), quals, k)
    }

    /// Histogram of wildcard run lengths in a real index (`PREMISE_DEBUG_INDEX`).
    #[test]
    #[ignore]
    fn debug_wildcard_runs() {
        let path = std::env::var("PREMISE_DEBUG_INDEX").unwrap();
        let idx = load_index(&std::fs::read(&path).unwrap(), &path).unwrap();
        let mut hist: std::collections::BTreeMap<usize, usize> = Default::default();
        let mut codes = [0usize; ALPHABET_SIZE];
        for j in 0..idx.num_sequences() {
            let text = idx.sequence(SeqId(j)).unwrap();
            let mut run = 0;
            for &b in text {
                codes[b as usize] += 1;
                if b >= N {
                    run += 1;
                } else if run > 0 {
                    *hist.entry(run).or_default() += 1;
                    run = 0;
                }
            }
            if run > 0 {
                *hist.entry(run).or_default() += 1;
            }
        }
        eprintln!("runs (len: count): {hist:?}");
        eprintln!("code counts: {codes:?}");
    }

    /// Micro-benchmark against a real index: `PREMISE_DEBUG_INDEX` (.fmidx),
    /// `PREMISE_DEBUG_FASTQ` (reads), `PREMISE_DEBUG_K`, optional `PREMISE_DEBUG_N`
    /// (reads to take, default 200). Prints ms per oriented read.
    #[test]
    #[ignore]
    fn bench_full_index_reads() {
        use haystackfm::alphabet::encode_byte;
        use std::io::BufRead;
        let path = std::env::var("PREMISE_DEBUG_INDEX").unwrap();
        let idx = load_index(&std::fs::read(&path).unwrap(), &path).unwrap();
        let k: usize = std::env::var("PREMISE_DEBUG_K").unwrap().parse().unwrap();
        let take: usize = std::env::var("PREMISE_DEBUG_N")
            .map(|v| v.parse().unwrap())
            .unwrap_or(200);
        let fq = std::fs::File::open(std::env::var("PREMISE_DEBUG_FASTQ").unwrap()).unwrap();
        let reads: Vec<Vec<u8>> = std::io::BufReader::new(fq)
            .lines()
            .skip(1)
            .step_by(4)
            .take(take)
            .map(|l| {
                l.unwrap()
                    .bytes()
                    .map(|b| encode_byte(b).filter(|&c| c <= N).unwrap_or(N))
                    .collect()
            })
            .collect();
        let t0 = std::time::Instant::now();
        let mut placements = 0usize;
        for p in &reads {
            let quals = vec![Q40; p.len()];
            let rc: Vec<u8> = p
                .iter()
                .rev()
                .map(|&c| if c <= 4 { 5 - c } else { c })
                .collect();
            for r in [p, &rc] {
                placements += score_read(&idx, r, &quals, k)
                    .values()
                    .map(|m| m.len())
                    .sum::<usize>();
            }
        }
        let per = t0.elapsed().as_secs_f64() * 1e3 / (2 * reads.len()) as f64;
        eprintln!(
            "k={k}: {:.3} ms per oriented read, {placements} placements over {} reads",
            per,
            reads.len()
        );
        use std::sync::atomic::Ordering::Relaxed;
        eprintln!(
            "extensions {} smems {} | walk nodes {}",
            counters::EXTENSIONS.load(Relaxed),
            counters::SMEMS.load(Relaxed),
            counters::WALK_NODES.load(Relaxed)
        );
    }

    /// Half a match: the read base is one of the two bases of a reference code.
    fn hm() -> f64 {
        set_m(2)
    }

    /// An index that was not sanitized holds IUPAC codes above `N`. A code never
    /// lies inside a qualifying run, but the walk scores it by its base set: a read
    /// base in the set is a fraction of a match, one outside it a mismatch.
    #[test]
    fn iupac_code_in_reference_scored_by_base_set() {
        let mut st = 37u64;
        let mut r = random_dna(60, &mut st);
        r.replace_range(20..21, "Y");
        let idx = index_raw(&[("r", &r)]);
        assert!(
            idx.sequence(SeqId(0)).unwrap().contains(&7),
            "Y must be present as code 7"
        );
        for (base, flank) in [("C", hm()), ("A", mm())] {
            let mut read_ascii = r[..50].to_string();
            read_ascii.replace_range(20..21, base);
            let read = enc(&read_ascii);
            let got = sr(&idx, &read, &q40(50), 11);
            let (d, ll) = single(&got, 0);
            assert_eq!(d, 0);
            assert_close(ll, 49.0 * m() + flank);
            assert_placements_match_oracle(&idx, &read, &q40(50), 11, "IUPAC ref");
        }
    }

    /// A reference code splits a run: with `k` longer than either exact side the
    /// placement does not qualify; with `k` at the left side's length it does, and
    /// the code is scored as a flank.
    #[test]
    fn run_split_by_reference_code_needs_an_exact_side() {
        let mut st = 41u64;
        let mut r = random_dna(60, &mut st);
        r.replace_range(20..21, "Y");
        let idx = index_raw(&[("r", &r)]);
        let mut read_ascii = r[..50].to_string();
        read_ascii.replace_range(20..21, "C");
        let read = enc(&read_ascii);

        assert!(sr(&idx, &read, &q40(50), 30).is_empty());
        assert_placements_match_oracle(&idx, &read, &q40(50), 30, "split, k 30");

        let got = sr(&idx, &read, &q40(50), 20);
        let (d, ll) = single(&got, 0);
        assert_eq!(d, 0);
        assert_close(ll, 49.0 * m() + hm());
        assert_placements_match_oracle(&idx, &read, &q40(50), 20, "split, k 20");
    }

    /// The reference `N` sits just before the only qualifying run: it is the run's
    /// left flank, scored as a quarter match.
    #[test]
    fn wildcard_before_run_is_its_flank() {
        let mut st = 43u64;
        let mut r = random_dna(60, &mut st);
        r.replace_range(20..21, "N");
        let idx = index_raw(&[("r", &r)]);
        let mut read_ascii = r[20..55].to_string();
        read_ascii.replace_range(0..1, "A");
        let read = enc(&read_ascii);
        let got = sr(&idx, &read, &q40(35), 30);

        let (d, ll) = single(&got, 0);
        assert_eq!(d, 20);
        assert_close(ll, 34.0 * m() + nm());
        assert_placements_match_oracle(&idx, &read, &q40(35), 30, "wildcard before run");
    }

    /// A cluster of reference codes is walked symbol by symbol, each scored by its
    /// own base set.
    #[test]
    fn wildcard_cluster_is_walked() {
        let mut st = 47u64;
        let mut r = random_dna(60, &mut st);
        r.replace_range(20..23, "NYN");
        let idx = index_raw(&[("r", &r)]);
        let mut read_ascii = r[..50].to_string();
        read_ascii.replace_range(20..23, "ACG");
        let read = enc(&read_ascii);
        let got = sr(&idx, &read, &q40(50), 20);

        let (d, ll) = single(&got, 0);
        assert_eq!(d, 0);
        assert_close(ll, 47.0 * m() + 2.0 * nm() + hm());
        assert_placements_match_oracle(&idx, &read, &q40(50), 20, "wildcard cluster");
    }

    /// A reference code at every offset of a `k`-window splits it: a read that is
    /// exactly the window never qualifies, a longer read qualifies on the exact
    /// side and scores the code as a half match.
    #[test]
    fn wildcard_at_every_offset_splits_the_window() {
        let k = 11;
        for off in 0..k {
            let mut st = 53 + off as u64;
            let mut r = random_dna(40, &mut st);
            r.replace_range(10 + off..11 + off, "R");
            let idx = index_raw(&[("r", &r)]);

            let mut window = r[10..10 + k].to_string();
            window.replace_range(off..off + 1, "A");
            let read = enc(&window);
            assert!(sr(&idx, &read, &q40(k), k).is_empty(), "offset {off}");

            let mut long = r[5..35].to_string();
            long.replace_range(5 + off..6 + off, "A");
            let read = enc(&long);
            let got = sr(&idx, &read, &q40(30), k);
            assert_close(*got[&SeqId(0)][&5], 29.0 * m() + hm());
            let ctx = format!("wildcard at offset {off}");
            assert_placements_match_oracle(&idx, &read, &q40(30), k, &ctx);
            assert_smems_match_oracle(&idx, &read, &q40(30), k, &ctx);
        }
    }

    /// Copies of a reference that differ by one substitution inside the read carry
    /// runs that are sub-matches of the base's full-length match, so under SMEM
    /// seeding only the base qualifies. Removing the base makes every copy's runs
    /// SMEMs again.
    #[test]
    fn shadowed_copies_get_no_placement() {
        let mut st = 61u64;
        let base = random_dna(80, &mut st);
        let mut copies: Vec<(String, String)> = Vec::new();
        for c in 0..3 {
            let mut v = base.clone();
            let at = 30 + c;
            let sub = match &v[at..at + 1] {
                "A" => "C",
                "C" => "G",
                "G" => "T",
                _ => "A",
            };
            v.replace_range(at..at + 1, sub);
            copies.push((format!("copy{c}"), v));
        }
        let read = enc(&base[10..60]);
        let quals = q40(50);

        let mut with_base: Vec<(String, String)> = vec![("base".into(), base.clone())];
        with_base.extend(copies.iter().cloned());
        for (label, refs) in [("with base", with_base), ("copies only", copies.clone())] {
            let borrowed: Vec<(&str, &str)> =
                refs.iter().map(|(h, s)| (h.as_str(), s.as_str())).collect();
            let idx = index_from(&borrowed);
            for k in [5, 11, 20] {
                let ctx = format!("{label}, k {k}");
                assert_smems_match_oracle(&idx, &read, &quals, k, &ctx);
                assert_placements_match_oracle(&idx, &read, &quals, k, &ctx);
                let got = sr(&idx, &read, &quals, k);
                match label {
                    // Only the exact match qualifies; every copy's runs are sub-matches.
                    "with base" => assert_eq!(got.len(), 1, "{ctx}: {got:?}"),
                    // Without it the SMEMs are `[0, 22)` (only copy2 matches that far)
                    // and `[21, 50)` (only copy0): copy1's runs are shadowed by both.
                    _ => {
                        assert_eq!(got.len(), 2, "{ctx}: {got:?}");
                        assert!(!got.contains_key(&SeqId(1)), "{ctx}: {got:?}");
                    }
                }
            }
        }
    }

    /// A diagonal whose read carries two SMEM runs is emitted once, by the leftmost;
    /// a reference too short for the read gets nothing even though it holds the run.
    #[test]
    fn ownership_and_bounds() {
        let mut st = 67u64;
        let base = random_dna(80, &mut st);
        // Two references, each matching one half of the read exactly, so the read has
        // two SMEMs and every diagonal carrying both is owned by the leftmost.
        let mut left = base.clone();
        let mut right = base.clone();
        left.replace_range(45..46, if &base[45..46] == "A" { "C" } else { "A" });
        right.replace_range(25..26, if &base[25..26] == "A" { "C" } else { "A" });
        // A short reference holding only the tail of the read: the run fits, the read does not.
        let short = base[40..80].to_string();
        let idx = index_from(&[("left", &left), ("right", &right), ("short", &short)]);
        let read = enc(&base[10..60]);
        let quals = q40(50);
        for k in [5, 11, 14] {
            let ctx = format!("k {k}");
            assert_smems_match_oracle(&idx, &read, &quals, k, &ctx);
            assert_placements_match_oracle(&idx, &read, &quals, k, &ctx);
        }
        let got = sr(&idx, &read, &quals, 11);
        assert!(!got.contains_key(&SeqId(2)), "{got:?}");
        assert_close(*got[&SeqId(0)][&10], 49.0 * m() + mm());
        assert_close(*got[&SeqId(1)][&10], 49.0 * m() + mm());
    }

    /// Sweep edge cases: a symbol absent from the text restarts the cursor, `k == 1`,
    /// `k == n`, `k > n`, identical references (LCP widening), and runs touching a
    /// reference boundary on either side.
    #[test]
    fn sweep_edge_cases() {
        let mut st = 71u64;
        let base = random_dna(60, &mut st);
        let ac: String = base
            .chars()
            .map(|c| if c == 'G' || c == 'T' { 'A' } else { c })
            .collect();
        let idx = index_from(&[
            ("ac", &ac),
            ("dup1", &base),
            ("dup2", &base),
            ("head", &base[..40]),
            ("tail", &base[20..]),
        ]);
        let mut read_absent = enc(&ac[5..30]);
        read_absent[12] = haystackfm::alphabet::G;
        let read_mid = enc(&base[10..50]);
        for (read, ks, ctx) in [
            (read_absent, vec![3usize, 6, 12], "symbol absent"),
            (read_mid, vec![1, 11, 40, 41], "boundaries"),
        ] {
            let quals = q40(read.len());
            for k in ks {
                let ctx = format!("{ctx}, k {k}");
                assert_smems_match_oracle(&idx, &read, &quals, k, &ctx);
                assert_placements_match_oracle(&idx, &read, &quals, k, &ctx);
            }
        }
    }

    /// The SMEM enumeration agrees with haystackfm's `find_smems` on ACGT-only input
    /// (where its compatibility rule is equality), both in intervals and occurrences.
    #[test]
    fn smems_agree_with_haystackfm_find_smems() {
        let mut st = 73u64;
        for case in 0..50 {
            let (refs, mut read, quals, k) = random_case_with(&mut st, &['A']);
            // Read `N` matches everything under haystackfm's rule and nothing here.
            read.iter_mut().for_each(|c| {
                if *c == N {
                    *c = A;
                }
            });
            let borrowed: Vec<(&str, &str)> =
                refs.iter().map(|(h, s)| (h.as_str(), s.as_str())).collect();
            let idx = index_raw(&borrowed);
            let (idx, read) = (&idx, &read);
            let n = read.len();
            let hits = idx.find_smems(read, k, true);
            let ctx = format!("find_smems case {case}: refs {refs:?} read {read:?} k {k}");

            let want: Vec<(usize, usize)> =
                hits.iter().map(|h| (h.query_start, h.query_end)).collect();
            let mut got = brute_force_smems(idx, read, k);
            got.sort_unstable();
            assert_eq!(got, want, "oracle intervals vs find_smems: {ctx}");

            let mut want: Vec<(u32, usize, usize, usize)> = hits
                .iter()
                .flat_map(|h| {
                    let (s, e) = (h.query_start, h.query_end);
                    h.positions.iter().filter_map(move |&(id, pos)| {
                        let len = idx.sequence(id).unwrap().len();
                        let pos = pos as usize;
                        (pos >= s && pos - s + n <= len).then_some((id.0, pos, s, e))
                    })
                })
                .collect();
            want.sort_unstable();
            want.dedup();
            assert_eq!(
                brute_force_smem_occurrences(idx, read, k),
                want,
                "oracle occurrences vs find_smems: {ctx}"
            );
            assert_smems_match_oracle(idx, read, &quals, k, &ctx);
        }
    }

    /// Seeding from the k-mer lookup tables yields exactly the SMEMs (starts, cursors and
    /// lengths) of the symbol-by-symbol search, on IUPAC references and reads with `N`.
    #[test]
    fn lookup_seeding_matches_stepwise_smems() {
        const WILD: [char; 11] = ['N', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V'];
        let mut st = 0x100Cu64;
        let mut jumps = 0usize;
        for case in 0..200 {
            let (refs, read, _, k) = random_case_with(&mut st, &WILD);
            let borrowed: Vec<(&str, &str)> =
                refs.iter().map(|(h, s)| (h.as_str(), s.as_str())).collect();
            let plain = index_raw(&borrowed);
            let want = smem_intervals(&plain, &read, k);
            for depth in [1, 3, 5, 8] {
                let tabled = index_raw_depth(&borrowed, depth);
                assert_eq!(tabled.lookup_depth(), depth);
                jumps += (0..read.len())
                    .filter(|&at| lookup(&tabled, &read, at).is_some())
                    .count();
                assert_eq!(
                    smem_intervals(&tabled, &read, k),
                    want,
                    "case {case}, depth {depth}: refs {refs:?} read {read:?} k {k}"
                );
            }
        }
        assert!(jumps > 0, "the lookup tables were never hit");
    }

    /// Random SMEM enumeration against the oracle, N-only and all-IUPAC references.
    #[test]
    fn smems_vs_brute_force_random() {
        const WILD: [char; 11] = ['N', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V'];
        let mut st = 0x5EEDu64;
        for (label, wild) in [("N", &WILD[..1]), ("iupac", &WILD[..])] {
            for case in 0..200 {
                let (refs, read, quals, k) = random_case_with(&mut st, wild);
                let borrowed: Vec<(&str, &str)> =
                    refs.iter().map(|(h, s)| (h.as_str(), s.as_str())).collect();
                let idx = index_raw(&borrowed);
                let ctx = format!("{label} case {case}: refs {refs:?} read {read:?} k {k}");
                assert_smems_match_oracle(&idx, &read, &quals, k, &ctx);
            }
        }
    }

    #[test]
    fn differential_vs_brute_force_random() {
        let mut st = 0xC0FFEEu64;
        for case in 0..200 {
            let (refs, read, quals, k) = random_case(&mut st);
            let borrowed: Vec<(&str, &str)> =
                refs.iter().map(|(h, s)| (h.as_str(), s.as_str())).collect();
            let idx = index_from(&borrowed);
            let got = sr(&idx, &read, &quals, k);
            let ctx = format!("case {case}: refs {refs:?} read {read:?} k {k}");
            assert_placements_match_oracle(&idx, &read, &quals, k, &ctx);
            assert_smems_match_oracle(&idx, &read, &quals, k, &ctx);
            for (id, m) in &got {
                let len = idx.sequence(*id).unwrap().len();
                for d in m.keys() {
                    assert!(d + read.len() <= len, "diagonal out of bounds: {ctx}");
                }
            }
        }
    }

    /// The same differential on unsanitized indexes holding every IUPAC code, so
    /// seeding and walking are exercised across ambiguity codes above `N`.
    #[test]
    fn differential_vs_brute_force_random_iupac() {
        const WILD: [char; 11] = ['N', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V'];
        let mut st = 0xBADC0DEu64;
        let mut wild_positions = 0usize;
        for case in 0..200 {
            let (refs, read, quals, k) = random_case_with(&mut st, &WILD);
            let borrowed: Vec<(&str, &str)> =
                refs.iter().map(|(h, s)| (h.as_str(), s.as_str())).collect();
            let idx = index_raw(&borrowed);
            wild_positions += (0..idx.num_sequences())
                .map(|j| {
                    idx.sequence(SeqId(j))
                        .unwrap()
                        .iter()
                        .filter(|&&b| b >= N)
                        .count()
                })
                .sum::<usize>();
            let ctx = format!("iupac case {case}: refs {refs:?} read {read:?} k {k}");
            assert_placements_match_oracle(&idx, &read, &quals, k, &ctx);
        }
        assert!(
            wild_positions > 100,
            "too few wildcards drawn: {wild_positions}"
        );
    }
}

#[cfg(test)]
mod query_tests {
    use super::*;
    use crate::io::{build_index_from_bytes, load_index};
    use bio::io::fastq;
    // Single-line version of the fixture ref_A sequence (156 bp).
    const TEST_FASTA: &[u8] =
        b">ref_A\nAGCTAGCTAGCTAGCTTACGATCGATCGAATCGAATCGATCGATCGATCGATCGATCGAATCGATCGATCGAATCGATCGATCGATCGAATCGATCGATCGAATCGATCGATCGAATCGATCGATCGAATCGATCGATCGAATCGATCGATCGAAT\n";

    fn build_test_index(fasta: &[u8]) -> RefIndex {
        let (bytes, _) =
            build_index_from_bytes(fasta).expect("build_index_from_bytes failed in test");
        load_index(&bytes, "<test>").expect("load_index failed in test")
    }

    /// [`query_read`] on a FASTQ record.
    fn qr(
        fmidx: &RefIndex,
        rec: &fastq::Record,
        k: usize,
        complement: bool,
    ) -> Result<MatchLikelihoods> {
        query_read(fmidx, rec.seq(), rec.qual(), k, complement)
    }

    fn make_record(id: &str, seq: &[u8]) -> fastq::Record {
        fastq::Record::with_attrs(id, None, seq, &vec![b'I'; seq.len()])
    }

    /// `merge_mates` (the pipeline's separable merge) must agree with
    /// `merge_read_pairs` (the cross product with its orientation test) on every
    /// input: the orientation test admits every placement pair for reads of positive
    /// length, so the sum factorises and the best pair is the best of each mate.
    #[test]
    fn separable_merge_matches_cross_product_merge() {
        let mut state = 0x1234_5678_9abc_def0u64;
        let mut next = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            state >> 33
        };
        for _ in 0..500 {
            let n_f = 1 + (next() % 5) as usize;
            let n_r = 1 + (next() % 5) as usize;
            let read_len = 1 + (next() % 300) as usize;
            let mut fwd: HashMap<usize, LogProb> = HashMap::new();
            let mut rev: HashMap<usize, LogProb> = HashMap::new();
            for _ in 0..n_f {
                fwd.insert(
                    (next() % 2000) as usize,
                    LogProb(-((next() % 400) as f64) / 7.0),
                );
            }
            for _ in 0..n_r {
                rev.insert(
                    (next() % 2000) as usize,
                    LogProb(-((next() % 400) as f64) / 7.0),
                );
            }
            let (want_ll, want_best, want_pos) = merge_read_pairs(&fwd, &rev, read_len);
            let (got_ll, got_best, got_pos) = merge_mates(&fwd, &rev);
            assert!(
                (*want_ll - *got_ll).abs() <= 1e-9 * want_ll.abs().max(1.0),
                "summed likelihood: {want_ll:?} vs {got_ll:?}"
            );
            assert_eq!(want_best, got_best, "best pair log-likelihood");
            // Positions agree unless the best score is tied between placements.
            let tied = fwd
                .values()
                .filter(|&&l| l == want_best - rev[&want_pos.1])
                .count()
                > 1
                || rev
                    .values()
                    .filter(|&&l| l == want_best - fwd[&want_pos.0])
                    .count()
                    > 1;
            if !tied {
                assert_eq!(want_pos, got_pos, "best pair positions");
            }
        }
    }

    /// A read matching the final bases of a reference aligns at
    /// `ref_pos + read_len == ref_len`, which the right-edge guard must accept.
    /// The bound was `<`, which silently discarded every alignment flush against
    /// the reference's 3' end.
    #[test]
    fn query_read_aligns_read_flush_against_reference_end() {
        let fmidx = build_test_index(TEST_FASTA);
        let ref_len = fmidx.sequence(SeqId::new(0)).unwrap().len();
        let read_seq = b"ATCGATCGAATCGATCGATCGAATCGATCGATCGAATCGATCGATCGAAT";
        let expected_pos = ref_len - read_seq.len();

        let hits = qr(&fmidx, &make_record("r", read_seq), 5, false).expect("query_read failed");

        let positions = hits
            .get(&SeqId::new(0))
            .expect("no alignments to ref_A for a read taken verbatim from its 3' end");
        assert!(
            positions.contains_key(&expected_pos),
            "expected an alignment at ref_pos {expected_pos} (flush with ref_len {ref_len}); \
             got positions: {:?}",
            positions.keys().collect::<Vec<_>>()
        );
    }

    /// BUG: N at position 15, kmer_size=16.
    ///
    /// All kmers at positions 0-15 contain N and are filtered.  The first valid
    /// kmer is at actual read position 16, but enumerate() (applied after filter)
    /// assigns it index 0 → MEMPos.read_start = 0 → ref_pos = 16 - 0 = 16 (wrong).
    ///
    /// Correct behaviour: alignment at ref position 0.
    /// Bug behaviour:     alignment at ref position 16, or read unaligned because
    ///                    matching probability at that wrong position ≈ 0.
    #[test]
    fn query_read_n_at_15_produces_alignment_at_position_0() {
        let fmidx = build_test_index(TEST_FASTA);
        let record = make_record("r", b"AGCTAGCTAGCTAGCNTACGATCGATCGAATCGAATCGATCGATCGATCG");
        let hits = qr(&fmidx, &record, 5, false).expect("query_read failed");

        let positions = hits.get(&SeqId::new(0)).unwrap_or_else(|| {
            panic!(
                "no alignments to ref_A — N-filtering bug likely caused alignment \
                 probability at the wrong position to underflow to 0"
            )
        });

        assert!(
            positions.contains_key(&0),
            "expected alignment at ref position 0 (read = ref_A[0..50] with N at pos 15);\n\
             got positions: {:?}\n\
             BUG: enumerate() after filter indexed first valid kmer (actual pos 16) as 0, \
             placing ref_pos at 16 instead of 0.",
            positions.keys().collect::<Vec<_>>()
        );
    }

    /// BUG: N at position 25, kmer_size=16.
    ///
    /// Kmers at positions 10-25 are filtered.  Valid kmers before N (0-9) get
    /// correct enumerate indices 0-9.  Valid kmers after N (26-34) get indices
    /// 10-18 instead of 26-34, so ref_pos = 26 - 10 = 16 (wrong; correct is 0).
    #[test]
    fn query_read_n_at_25_produces_alignment_at_position_0() {
        let fmidx = build_test_index(TEST_FASTA);
        let record = make_record("r", b"AGCTAGCTAGCTAGCTTACGATCGANCGAATCGAATCGATCGATCGATCG");
        let hits = qr(&fmidx, &record, 5, false).expect("query_read failed");

        let positions = hits
            .get(&SeqId::new(0))
            .expect("expected at least one alignment to ref_A");

        assert!(
            positions.contains_key(&0),
            "expected alignment at ref position 0; got positions: {:?}\n\
             BUG: post-N kmers (26-34) got enumerate indices 10-18, producing a \
             spurious hit at position 16 with no correct hit at position 0.",
            positions.keys().collect::<Vec<_>>()
        );
    }

    /// read_seq = revcomp(record.seq()).  N at record position 34 → N at position
    /// 50-1-34 = 15 of read_seq, triggering the same kmer-filtering shift.
    #[test]
    fn query_read_complement_n_at_record_pos_34_produces_alignment_at_zero() {
        let fmidx = build_test_index(TEST_FASTA);
        let ref_a_prefix: &[u8] = b"AGCTAGCTAGCTAGCTTACGATCGATCGAATCGAATCGATCGATCGATCG";
        let rc = bio::alphabets::dna::revcomp(ref_a_prefix);
        // N at record pos 34 → pos 15 of revcomp(record), triggering the bug.
        let mut seq_with_n = rc.clone();
        seq_with_n[34] = b'N';

        let record = make_record("r_rc", &seq_with_n);
        let hits = qr(&fmidx, &record, 5, true).expect("query_read (complement) failed");

        let positions = hits.get(&SeqId::new(0)).unwrap_or_else(|| {
            panic!(
                "no alignments to ref_A via complement path — N-filtering bug may have \
                 caused alignment probability at the wrong position to underflow"
            )
        });

        assert!(
            positions.contains_key(&0),
            "expected alignment at ref position 0 via complement path; \
             got positions: {:?}\n\
             BUG: same enumerate() shift applies in the complement path.",
            positions.keys().collect::<Vec<_>>()
        );
    }

    /// IUPAC letters in a read encode to codes above `N`, which no sanitized
    /// reference contains; the wrapper maps them to `N` so they score as a
    /// quarter-match instead of silently killing the placement.
    #[test]
    fn iupac_letters_in_read_treated_as_n() {
        let fmidx = build_test_index(TEST_FASTA);
        let with_r = make_record("r", b"AGCTAGCTAGCTAGCRTACGATCGATCGAATCGAATCGATCGATCGATCG");
        let with_n = make_record("n", b"AGCTAGCTAGCTAGCNTACGATCGATCGAATCGAATCGATCGATCGATCG");
        let hits_r = qr(&fmidx, &with_r, 11, false).unwrap();
        let hits_n = qr(&fmidx, &with_n, 11, false).unwrap();
        let ll_r = hits_r[&SeqId::new(0)][&0];
        let ll_n = hits_n[&SeqId::new(0)][&0];
        assert!(
            (*ll_r - *ll_n).abs() < 1e-12,
            "R should score like N: {ll_r:?} vs {ll_n:?}"
        );
    }
}
