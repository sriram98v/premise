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

/// Compute probability of match given ref is true source
///
/// Both sequences are in haystackfm's alphabet code space (`A = 1`, `C = 2`, `G = 3`,
/// `T = 4`, `N = 5`), not ASCII — the index stores and serves reference bases that way,
/// and reads are encoded once per orientation before seeding.
///
/// N bases in the read or reference are treated as uninformative and contribute
/// 0.0 to the log-likelihood (i.e., probability 1.0 for that position).
pub fn compute_match_log_prob(
    q_seq: &[u8],
    quality_score_vec: &[u8],
    aligned_ref_seq: &[u8],
) -> LogProb {
    let mut match_log_likelihood = 0_f64;
    for (read_char, reference_char, quality_score) in
        izip!(q_seq.iter(), aligned_ref_seq.iter(), quality_score_vec)
    {
        if *read_char == alphabet::N || *reference_char == alphabet::N {
            // N is ambiguous — skip this position (multiply likelihood by 1.0)
            continue;
        }
        match read_char == reference_char {
            true => match_log_likelihood += (1_f64 - *error_prob(*quality_score)).ln(),
            false => match_log_likelihood += (*error_prob(*quality_score) * (1_f64 / 3_f64)).ln(),
        }
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

    // Quality byte for Phred 40 ('I' = 73 in ASCII, Q = 73 - 33 = 40)
    const Q40: u8 = b'I';

    /// Bases in the alphabet's code space, which is what `compute_match_log_prob` compares.
    /// Written as ASCII here for readability and encoded once, so the fixtures stay legible.
    fn bases(ascii: &str) -> Vec<u8> {
        ascii
            .chars()
            .map(|c| alphabet::encode_char(c).expect("test fixture uses IUPAC bases"))
            .collect()
    }

    #[test]
    fn n_in_read_skipped_contributes_zero_to_log_prob() {
        // A read that is identical to the reference except position 0 is N
        // should score the same as if that position were absent.
        let read = bases("NACGT");
        let reference = bases("AACGT");
        let quals = [Q40; 5];

        let with_n = compute_match_log_prob(&read, &quals, &reference);

        // Expected: only the 4 non-N positions contribute, all matches
        let expected = compute_match_log_prob(&bases("ACGT"), &[Q40; 4], &bases("ACGT"));

        assert!(
            (*with_n - *expected).abs() < 1e-10,
            "N in read should be skipped: got {with_n:?}, expected {expected:?}"
        );
    }

    #[test]
    fn n_in_reference_skipped_contributes_zero_to_log_prob() {
        let read = bases("AACGT");
        let reference = bases("NACGT");
        let quals = [Q40; 5];

        let with_n = compute_match_log_prob(&read, &quals, &reference);
        let expected = compute_match_log_prob(&bases("ACGT"), &[Q40; 4], &bases("ACGT"));

        assert!(
            (*with_n - *expected).abs() < 1e-10,
            "N in reference should be skipped: got {with_n:?}, expected {expected:?}"
        );
    }

    #[test]
    fn all_n_read_returns_zero_log_prob() {
        let read = bases("NNNNN");
        let reference = bases("AACGT");
        let quals = [Q40; 5];

        let result = compute_match_log_prob(&read, &quals, &reference);
        assert_eq!(
            *result, 0.0,
            "All-N read should yield log-prob 0.0 (prob 1.0)"
        );
    }

    #[test]
    fn no_n_read_behavior_unchanged() {
        // Regression: reads without N should behave identically to before.
        let read = bases("AACGT");
        let reference = bases("AACGT");
        let quals = [Q40; 5];

        let result = compute_match_log_prob(&read, &quals, &reference);
        // All matches at Q40: each position contributes ln(1 - 10^-4) ≈ -1e-4
        assert!(
            *result < 0.0 && *result > -1.0,
            "Perfect match log-prob should be slightly negative: {result:?}"
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
