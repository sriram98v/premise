//! Types and helpers shared by the alignment ([`crate::align`]), EM ([`crate::em`]) and
//! I/O ([`crate::io`]) stages

use crate::SeqId;
use bio::io::fastq;
use bio::stats::{LogProb, Prob};
use haystackfm::{alphabet, SymbolSet};
use itertools::izip;
use std::collections::HashMap;
use std::ffi::OsStr;
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::LazyLock;
use std::time::{SystemTime, UNIX_EPOCH};

/// Floating-point type used for all EM probabilities and log-probabilities in linear space.
pub type EMProb = f64;

/// Newtype wrapper around a raw read identifier string.
#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct ReadID(pub(crate) String);

impl std::ops::Deref for ReadID {
    type Target = String;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl std::fmt::Display for ReadID {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

/// Dense integer index assigned to each read for use as a [`ReadAlignments`] row.
#[derive(Debug, Clone, Hash, PartialEq, Eq, Copy)]
pub struct ReadIdx(pub(crate) usize);

impl ReadIdx {
    /// Construct a dense read index. Primarily for tests and benchmarks.
    pub fn new(n: usize) -> Self {
        Self(n)
    }
}

impl std::ops::Deref for ReadIdx {
    type Target = usize;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

pub struct MatchLikelihoods(pub(crate) HashMap<SeqId, HashMap<usize, LogProb>>);

impl MatchLikelihoods {
    /// Create an empty `MatchLikelihoods` map.
    pub fn new() -> Self {
        Self(HashMap::new())
    }
}

impl std::ops::Deref for MatchLikelihoods {
    type Target = HashMap<SeqId, HashMap<usize, LogProb>>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl std::ops::DerefMut for MatchLikelihoods {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

#[derive(Debug, Clone, Copy)]
pub struct ReadAlignment {
    pub(crate) full_match_ll: LogProb,
    pub(crate) best_match_ll: LogProb,
    pub(crate) pos: (usize, usize),
    pub(crate) complement: bool,
}

impl ReadAlignment {
    pub fn increment_full_match_ll(&mut self, match_prob: LogProb) {
        self.full_match_ll = LogProb::from(Prob(self.full_match_ll.exp() + match_prob.exp()));
    }

    pub fn update_best_match_ll(&mut self, match_ll: LogProb) {
        self.best_match_ll = match_ll;
    }

    pub fn update_pos(&mut self, pos: (usize, usize)) {
        self.pos = pos;
    }

    pub fn update_compl(&mut self, compl: bool) {
        self.complement = compl;
    }

    pub fn get_full_match_ll(&self) -> LogProb {
        self.full_match_ll
    }

    pub fn get_best_match_ll(&self) -> LogProb {
        self.best_match_ll
    }

    pub fn get_pos(&self) -> (usize, usize) {
        self.pos
    }

    pub fn get_compl(&self) -> bool {
        self.complement
    }
}

/// Alignment results of every read that kept at least one reference
#[derive(Debug)]
pub struct ReadAlignments<'a> {
    reads: Vec<&'a ReadID>,
    row_ptr: Vec<usize>,
    refs: Vec<SeqId>,
    aligns: Vec<ReadAlignment>,
    index: HashMap<&'a ReadID, ReadIdx>,
}

impl<'a> ReadAlignments<'a> {
    /// Assemble from one row per read pair, each sorted by `SeqId`. Reads whose row
    /// is empty (no reference survived the cutoffs) get no [`ReadIdx`].
    pub(crate) fn from_rows(rows: Vec<(&'a ReadID, Vec<(SeqId, ReadAlignment)>)>) -> Self {
        let n_rows = rows.iter().filter(|(_, row)| !row.is_empty()).count();
        let nnz: usize = rows.iter().map(|(_, row)| row.len()).sum();
        let mut reads = Vec::with_capacity(n_rows);
        let mut row_ptr = Vec::with_capacity(n_rows + 1);
        let mut refs = Vec::with_capacity(nnz);
        let mut aligns = Vec::with_capacity(nnz);
        let mut index = HashMap::with_capacity(n_rows);
        row_ptr.push(0);
        for (read_id, row) in rows {
            if row.is_empty() {
                continue;
            }
            debug_assert!(row.windows(2).all(|w| w[0].0 .0 < w[1].0 .0));
            index.insert(read_id, ReadIdx(reads.len()));
            reads.push(read_id);
            for (ref_idx, align) in row {
                refs.push(ref_idx);
                aligns.push(align);
            }
            row_ptr.push(refs.len());
        }
        Self {
            reads,
            row_ptr,
            refs,
            aligns,
            index,
        }
    }

    /// Number of reads with at least one alignment.
    pub fn len(&self) -> usize {
        self.reads.len()
    }

    pub fn is_empty(&self) -> bool {
        self.reads.is_empty()
    }

    /// Number of stored (read, reference) entries.
    pub fn nnz(&self) -> usize {
        self.refs.len()
    }

    /// The row index of a read, if it has any alignment.
    pub fn read_idx(&self, read_id: &ReadID) -> Option<ReadIdx> {
        self.index.get(read_id).copied()
    }

    pub fn contains(&self, read_id: &ReadID) -> bool {
        self.index.contains_key(read_id)
    }

    pub fn read_id(&self, read_idx: ReadIdx) -> &'a ReadID {
        self.reads[read_idx.0]
    }

    /// The references a read aligned to (sorted by id) and the matching alignments.
    pub fn row(&self, read_idx: ReadIdx) -> (&[SeqId], &[ReadAlignment]) {
        let (s, e) = (self.row_ptr[read_idx.0], self.row_ptr[read_idx.0 + 1]);
        (&self.refs[s..e], &self.aligns[s..e])
    }

    /// The alignment of a read against one reference, if any.
    pub fn get(&self, read_idx: ReadIdx, ref_idx: SeqId) -> Option<&ReadAlignment> {
        let (refs, aligns) = self.row(read_idx);
        refs.binary_search_by_key(&ref_idx.0, |r| r.0)
            .ok()
            .map(|k| &aligns[k])
    }

    /// Every row as `(read index, read id, references, alignments)`.
    pub fn iter(&self) -> impl Iterator<Item = (ReadIdx, &'a ReadID, &[SeqId], &[ReadAlignment])> {
        (0..self.reads.len()).map(move |r| {
            let (refs, aligns) = self.row(ReadIdx(r));
            (ReadIdx(r), self.reads[r], refs, aligns)
        })
    }
}

/// One read's bases and Phred+33 qualities, as two boxed byte slices
#[derive(Debug, Clone)]
pub struct Read {
    seq: Box<[u8]>,
    qual: Box<[u8]>,
}

impl Read {
    pub fn new(seq: &[u8], qual: &[u8]) -> Self {
        Self {
            seq: seq.into(),
            qual: qual.into(),
        }
    }

    pub fn seq(&self) -> &[u8] {
        &self.seq
    }

    pub fn qual(&self) -> &[u8] {
        &self.qual
    }

    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }
}

impl From<fastq::Record> for Read {
    fn from(record: fastq::Record) -> Self {
        Read::new(record.seq(), record.qual())
    }
}

/// A matched pair of forward (R1) and reverse (R2) reads sharing a common read ID.
pub struct ReadPair {
    pub(crate) read_id: ReadID,
    pub(crate) r1: Read,
    pub(crate) r2: Read,
}

impl ReadPair {
    /// Construct a read pair from an ID and its two mate FASTQ records.
    pub fn new(read_id: &str, r1: fastq::Record, r2: fastq::Record) -> Self {
        Self {
            read_id: ReadID(read_id.to_string()),
            r1: r1.into(),
            r2: r2.into(),
        }
    }
}

/// Real-time progress counters for a running web query, exposed via `/api/query/progress`.
pub struct QueryProgress {
    pub phase: AtomicU64,
    pub reads_done: AtomicU64,
    pub reads_total: AtomicU64,
    pub em_iter_done: AtomicU64,
    pub em_iter_total: AtomicU64,
    pub started_ms: u64,
}

impl QueryProgress {
    pub fn new() -> Self {
        let started_ms = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map(|d| d.as_millis() as u64)
            .unwrap_or(0);
        Self {
            phase: AtomicU64::new(0),
            reads_done: AtomicU64::new(0),
            reads_total: AtomicU64::new(0),
            em_iter_done: AtomicU64::new(0),
            em_iter_total: AtomicU64::new(0),
            started_ms,
        }
    }

    pub fn to_json(&self) -> String {
        let now_ms = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map(|d| d.as_millis() as u64)
            .unwrap_or(0);
        let elapsed_ms = now_ms.saturating_sub(self.started_ms);
        format!(
            r#"{{"phase":{},"reads_done":{},"reads_total":{},"em_iter_done":{},"em_iter_total":{},"elapsed_ms":{}}}"#,
            self.phase.load(Ordering::Relaxed),
            self.reads_done.load(Ordering::Relaxed),
            self.reads_total.load(Ordering::Relaxed),
            self.em_iter_done.load(Ordering::Relaxed),
            self.em_iter_total.load(Ordering::Relaxed),
            elapsed_ms,
        )
    }
}

/// The `{A, C, G, T}` bases an IUPAC symbol stands for; empty for the sentinel and for
/// codes outside the alphabet.
fn iupac_set(code: u8) -> SymbolSet {
    SymbolSet::from_codes(alphabet::iupac_bases(code))
}

/// Per-position mixture weights `(w_match, w_mismatch_scaled)` for every ordered
/// `(read symbol, reference symbol)` pair, flattened as `read_code * 16 + ref_code`.
static AMBIG_WEIGHTS: LazyLock<[(f64, f64); 256]> = LazyLock::new(|| {
    let mut table = [(0.0f64, 0.0f64); 256];
    for read_code in 0..alphabet::ALPHABET_SIZE {
        for ref_code in 0..alphabet::ALPHABET_SIZE {
            let read_set = iupac_set(read_code as u8);
            let ref_set = iupac_set(ref_code as u8);
            let n = read_set.len() * ref_set.len();
            if n == 0 {
                continue;
            }
            let m = read_set.intersection(ref_set).len();
            let total = f64::from(n);
            table[read_code * alphabet::ALPHABET_SIZE + ref_code] = (
                f64::from(m) / total,
                (f64::from(n - m) / total) * (1_f64 / 3_f64),
            );
        }
    }
    table
});

/// Compute probability of match given ref is true source
pub fn compute_match_log_prob(
    q_seq: &[u8],
    quality_score_vec: &[u8],
    aligned_ref_seq: &[u8],
) -> LogProb {
    let mut match_log_likelihood = 0_f64;
    for (read_char, reference_char, quality_score) in
        izip!(q_seq.iter(), aligned_ref_seq.iter(), quality_score_vec)
    {
        let idx = (*read_char as usize) * alphabet::ALPHABET_SIZE + (*reference_char as usize);
        let (w_match, w_mismatch) = match AMBIG_WEIGHTS.get(idx) {
            Some(&weights) => weights,
            // Read or reference symbol is outside the alphabet entirely.
            None => continue,
        };
        // Sentinel or unknown code on either side: no base set, nothing to condition on.
        if w_match == 0.0 && w_mismatch == 0.0 {
            continue;
        }
        let error_prob = *error_prob(*quality_score);
        match_log_likelihood += (w_match * (1_f64 - error_prob) + w_mismatch * error_prob).ln();
    }
    LogProb(match_log_likelihood)
}

/// Precomputed lookup table
static ERROR_PROB_LUT: LazyLock<[f64; 256]> = LazyLock::new(|| {
    let mut table = [0.0f64; 256];
    for (i, slot) in table.iter_mut().enumerate() {
        *slot = 10_f64.powf(-((i as f64 - 33.0) / 10_f64));
    }
    table
});

/// Convert a Phred+33 quality byte to its linear-space base-call error probability.
pub fn error_prob(q: u8) -> Prob {
    Prob(ERROR_PROB_LUT[q as usize])
}

/// Return the DNA complement of a sequence (A↔T, C↔G).
pub fn complement(q_seq: Vec<char>) -> Vec<char> {
    q_seq
        .iter()
        .map(|x| match x {
            'T' => 'A',
            'C' => 'G',
            'A' => 'T',
            'G' => 'C',
            'N' => 'N',
            _ => 'E',
        })
        .collect()
}

/// Extract the file extension from a filename, returning `None` if there is none.
pub fn get_extension_from_filename(filename: &str) -> Option<&str> {
    Path::new(filename).extension().and_then(OsStr::to_str)
}
#[cfg(test)]
mod tests {
    use super::*;

    const Q40: u8 = b'I';

    /// Bases in the alphabet's code space, which is what `compute_match_log_prob` compares.
    fn bases(ascii: &str) -> Vec<u8> {
        ascii
            .chars()
            .map(|c| alphabet::encode_char(c).expect("test fixture uses IUPAC bases"))
            .collect()
    }

    /// Score a single position, so the mixture expectations below stay readable.
    fn score_one(read: &str, reference: &str) -> f64 {
        *compute_match_log_prob(&bases(read), &[Q40], &bases(reference))
    }

    /// The mixture a symbol pair should produce: `m` of `n` base pairs agree.
    fn expected_mixture(m: u32, n: u32) -> f64 {
        let e = *error_prob(Q40);
        let (m, n) = (f64::from(m), f64::from(n));
        ((m / n) * (1.0 - e) + ((n - m) / n) * (e * (1.0 / 3.0))).ln()
    }

    #[test]
    fn r_in_reference_is_half_match_half_mismatch() {
        // R = {A, G}: a read A is a match under half the prior, a mismatch under the other.
        let e = *error_prob(Q40);
        let expected = (0.5 * (1.0 - e) + 0.5 * (e * (1.0 / 3.0))).ln();

        assert_eq!(
            score_one("A", "R"),
            expected,
            "read A vs reference R should be 0.5*P(match) + 0.5*P(mismatch)"
        );
        assert_eq!(
            score_one("G", "R"),
            expected,
            "read G vs reference R should score the same as read A"
        );
    }

    #[test]
    fn two_base_codes_score_one_half() {
        // (code, a base inside its set, a base outside it)
        let cases = [
            ("R", "A", "C"),
            ("Y", "C", "A"),
            ("S", "G", "A"),
            ("W", "T", "C"),
            ("K", "G", "A"),
            ("M", "C", "G"),
        ];
        for (code, inside, outside) in cases {
            assert_eq!(
                score_one(inside, code),
                expected_mixture(1, 2),
                "read {inside} vs reference {code} should be the 1-of-2 mixture"
            );
            assert_eq!(
                score_one(outside, code),
                expected_mixture(0, 2),
                "read {outside} is in no branch of {code}: full mismatch"
            );
        }
    }

    #[test]
    fn three_base_codes_score_one_third() {
        let cases = [
            ("B", "C", "A"),
            ("D", "A", "C"),
            ("H", "A", "G"),
            ("V", "A", "T"),
        ];
        for (code, inside, outside) in cases {
            assert_eq!(
                score_one(inside, code),
                expected_mixture(1, 3),
                "read {inside} vs reference {code} should be the 1-of-3 mixture"
            );
            assert_eq!(
                score_one(outside, code),
                expected_mixture(0, 3),
                "read {outside} is in no branch of {code}: full mismatch"
            );
        }
    }

    #[test]
    fn n_against_any_base_scores_one_quarter() {
        // N spreads the prior over all four bases, so the error term cancels exactly:
        // (1*(1-e) + 3*(e/3)) / 4 == 1/4, whatever the quality.
        let quarter = 0.25_f64.ln();
        for base in ["A", "C", "G", "T"] {
            assert!(
                (score_one(base, "N") - quarter).abs() < 1e-12,
                "reference N vs read {base} should be ln(0.25)"
            );
            assert!(
                (score_one("N", base) - quarter).abs() < 1e-12,
                "read N vs reference {base} should be ln(0.25)"
            );
        }
        assert!(
            (score_one("N", "N") - quarter).abs() < 1e-12,
            "N on both sides should still be ln(0.25)"
        );
    }

    #[test]
    fn ambiguity_on_both_sides_uses_the_product_prior() {
        // R vs R: 4 ordered base pairs, 2 of which agree (A/A and G/G).
        assert_eq!(score_one("R", "R"), expected_mixture(2, 4));
        // R = {A,G} vs Y = {C,T}: disjoint, so every pair is a mismatch.
        assert_eq!(score_one("R", "Y"), expected_mixture(0, 4));
        // R = {A,G} vs V = {A,C,G}: 6 pairs, A/A and G/G agree.
        assert_eq!(score_one("R", "V"), expected_mixture(2, 6));
    }

    #[test]
    fn exact_acgt_scoring_is_bit_identical() {
        let read = bases("AACGTACGT");
        let reference = bases("ACCGTAGGT");
        let quals = [Q40, b'#', b'5', Q40, b'J', b'!', Q40, b'2', b'F'];

        // The pre-mixture implementation, verbatim.
        let mut expected = 0_f64;
        for (r, f, q) in izip!(read.iter(), reference.iter(), quals.iter()) {
            expected += match r == f {
                true => (1_f64 - *error_prob(*q)).ln(),
                false => (*error_prob(*q) * (1_f64 / 3_f64)).ln(),
            };
        }

        assert_eq!(
            *compute_match_log_prob(&read, &quals, &reference),
            expected,
            "pure-ACGT scoring must not move by even one ULP"
        );
    }

    #[test]
    fn weights_match_iupac_definition() {
        // Independent oracle: the IUPAC code table written out by hand, so a change in
        // haystackfm's `iupac_bases` cannot silently move the likelihood model.
        const A: u8 = 0b0001;
        const C: u8 = 0b0010;
        const G: u8 = 0b0100;
        const T: u8 = 0b1000;
        let mask = |code: u8| match code {
            alphabet::A => A,
            alphabet::C => C,
            alphabet::G => G,
            alphabet::T => T,
            alphabet::R => A | G,
            alphabet::Y => C | T,
            alphabet::S => G | C,
            alphabet::W => A | T,
            alphabet::K => G | T,
            alphabet::M => A | C,
            alphabet::B => C | G | T,
            alphabet::D => A | G | T,
            alphabet::H => A | C | T,
            alphabet::V => A | C | G,
            alphabet::N => A | C | G | T,
            _ => 0,
        };
        for read_code in 0..alphabet::ALPHABET_SIZE as u8 {
            for ref_code in 0..alphabet::ALPHABET_SIZE as u8 {
                let (read_mask, ref_mask) = (mask(read_code), mask(ref_code));
                let n = read_mask.count_ones() * ref_mask.count_ones();
                let m = (read_mask & ref_mask).count_ones();
                let idx = read_code as usize * alphabet::ALPHABET_SIZE + ref_code as usize;
                if n == 0 {
                    assert_eq!(
                        AMBIG_WEIGHTS[idx],
                        (0.0, 0.0),
                        "read {read_code} vs ref {ref_code}"
                    );
                    continue;
                }

                let idx = read_code as usize * alphabet::ALPHABET_SIZE + ref_code as usize;
                let (w_match, w_mismatch) = AMBIG_WEIGHTS[idx];

                assert_eq!(
                    w_match,
                    f64::from(m) / f64::from(n),
                    "match weight wrong for read {read_code} vs ref {ref_code}"
                );
                assert_eq!(
                    w_mismatch,
                    (f64::from(n - m) / f64::from(n)) * (1_f64 / 3_f64),
                    "mismatch weight wrong for read {read_code} vs ref {ref_code}"
                );
            }
        }
    }

    #[test]
    fn sentinel_code_is_uninformative() {
        let quals = [Q40];
        assert_eq!(
            *compute_match_log_prob(&[alphabet::SENTINEL], &quals, &bases("A")),
            0.0,
            "sentinel in the read has no base set to condition on"
        );
        assert_eq!(
            *compute_match_log_prob(&bases("A"), &quals, &[alphabet::SENTINEL]),
            0.0,
            "sentinel in the reference has no base set to condition on"
        );
        assert_eq!(
            *compute_match_log_prob(&[200], &quals, &bases("A")),
            0.0,
            "a code outside the alphabet is skipped rather than indexing out of bounds"
        );
    }

    #[test]
    fn complement_n_maps_to_n() {
        let seq = vec!['A', 'N', 'C', 'G', 'T', 'N'];
        let result = complement(seq);
        assert_eq!(result, vec!['T', 'N', 'G', 'C', 'A', 'N']);
    }

    #[test]
    fn complement_standard_bases() {
        let seq = vec!['A', 'T', 'C', 'G'];
        assert_eq!(complement(seq), vec!['T', 'A', 'G', 'C']);
    }
}
