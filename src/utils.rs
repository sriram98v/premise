use bio::stats::{LogProb, Prob};
use haystackfm::alphabet;
use itertools::izip;
use std::ffi::OsStr;
use std::path::Path;
use std::sync::LazyLock;

/// Floating-point type used for all EM probabilities and log-probabilities in linear space.
///
/// Values in alignment likelihood maps are stored as raw `f64`; the [`bio::stats::LogProb`]
/// newtype is used at boundaries where log-space arithmetic is required.
pub type EMProb = f64;

/// The `{A, C, G, T}` set an IUPAC symbol stands for, as a bitmask
/// (bit 0 = A, bit 1 = C, bit 2 = G, bit 3 = T).
///
/// Indexed by haystackfm alphabet code, not ASCII. The sentinel and any code
/// outside the 16-symbol alphabet map to the empty set.
const fn iupac_mask(code: u8) -> u8 {
    const A: u8 = 0b0001;
    const C: u8 = 0b0010;
    const G: u8 = 0b0100;
    const T: u8 = 0b1000;
    match code {
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
    }
}

/// Per-position mixture weights `(w_match, w_mismatch_scaled)` for every ordered
/// `(read symbol, reference symbol)` pair, flattened as `read_code * 16 + ref_code`.
///
/// An IUPAC symbol denotes a uniform prior over the bases it stands for, so a pair
/// of symbols with base sets `S_r` and `S_R` gives
///
/// ```text
/// P = [ m * (1 - e) + (n - m) * (e / 3) ] / n     n = |S_r| * |S_R|,  m = |S_r n S_R|
/// ```
///
/// Reference `R` = {A, G} against read `A` is the `n = 2, m = 1` case:
/// `0.5 * P(match) + 0.5 * P(mismatch)`. Exact ACGT symbols give `n = 1`, which
/// collapses to the plain match/mismatch pair.
///
/// The `1/3` mismatch factor is folded into `w_mismatch_scaled` here rather than
/// applied in the scoring loop. That keeps the exact-mismatch case a single
/// `(1/3) * e` product, bit-identical to the pre-mixture implementation — writing
/// `w * (e / 3.0)` in the loop would not be, since `1/3` has no exact `f64`.
///
/// Pairs involving the sentinel or an unknown code have `n = 0` and are stored as
/// `(0.0, 0.0)`, which the scoring loop reads as "uninformative, skip".
static AMBIG_WEIGHTS: LazyLock<[(f64, f64); 256]> = LazyLock::new(|| {
    let mut table = [(0.0f64, 0.0f64); 256];
    for read_code in 0..alphabet::ALPHABET_SIZE {
        for ref_code in 0..alphabet::ALPHABET_SIZE {
            let read_mask = iupac_mask(read_code as u8);
            let ref_mask = iupac_mask(ref_code as u8);
            let n = read_mask.count_ones() * ref_mask.count_ones();
            if n == 0 {
                continue;
            }
            let m = (read_mask & ref_mask).count_ones();
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
///
/// Both sequences are in haystackfm's alphabet code space (`A = 1`, `C = 2`, `G = 3`,
/// `T = 4`, `N = 5`, `R`..`V` = 6..15), not ASCII — the index stores and serves
/// reference bases that way, and reads are encoded once per orientation before seeding.
///
/// Ambiguity codes on either side are scored as a uniform mixture over the bases they
/// stand for; see [`AMBIG_WEIGHTS`]. `N` is simply the four-base case and contributes
/// `ln(0.25)` per position rather than being treated as a free perfect match.
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

/// Precomputed lookup table mapping every possible Phred+33 quality byte to its
/// linear-space base-call error probability `10^(-((q - 33) / 10))`.
///
/// Built once on first use. On the alignment hot path `error_prob` is called
/// once per base of every diagonal of every read (×4 orientations); memoizing
/// the `powf` here removes it entirely from that loop. For valid quality bytes
/// (`q >= 33`) the table values are bit-identical to the direct formula.
static ERROR_PROB_LUT: LazyLock<[f64; 256]> = LazyLock::new(|| {
    let mut table = [0.0f64; 256];
    for (i, slot) in table.iter_mut().enumerate() {
        *slot = 10_f64.powf(-((i as f64 - 33.0) / 10_f64));
    }
    table
});

/// Convert a Phred+33 quality byte to its linear-space base-call error probability.
///
/// Uses the standard formula: P(error) = 10^(-(Q - 33) / 10), served from a
/// precomputed 256-entry lookup table ([`ERROR_PROB_LUT`]).
pub fn error_prob(q: u8) -> Prob {
    Prob(ERROR_PROB_LUT[q as usize])
}

/// Return the DNA complement of a sequence (A↔T, C↔G).
///
/// N bases are mapped to `'N'` (ambiguous complement of ambiguous).
/// Other non-ACGTN characters are mapped to `'E'` to signal an unexpected base.
/// Note: this does not reverse the sequence; for reverse-complement use
/// [`bio::alphabets::dna::revcomp`] instead.
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
    /// Written as ASCII here for readability and encoded once, so the fixtures stay legible.
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
    fn weights_match_haystackfm_iupac_bases() {
        // Independent oracle: recompute the intersection from haystackfm's own table
        // rather than from our bitmask, so a typo in `iupac_mask` cannot hide.
        for read_code in 1..alphabet::ALPHABET_SIZE as u8 {
            for ref_code in 1..alphabet::ALPHABET_SIZE as u8 {
                let read_set = alphabet::iupac_bases(read_code);
                let ref_set = alphabet::iupac_bases(ref_code);
                let n = (read_set.len() * ref_set.len()) as u32;
                let m = read_set.iter().filter(|b| ref_set.contains(b)).count() as u32;

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
