use bio::stats::{LogProb, Prob};
use itertools::izip;
use std::path::Path;
use std::{io::Write, fs};
use std::collections::HashMap;
use std::ffi::OsStr;

use crate::EMProb;

/// Compute probability of match given ref is true source
pub fn compute_match_log_prob(q_seq: &[u8], quality_score_vec: &[u8], aligned_ref_seq: &[u8]) -> LogProb{
    let mut match_log_likelihood = 0_f64;
    for (read_char,reference_char,quality_score) in izip!(q_seq.iter(), aligned_ref_seq.iter(), quality_score_vec){
        match read_char==reference_char{
            true => match_log_likelihood += ((1_f64-*error_prob(*quality_score))).ln(),
            false => match_log_likelihood += (*error_prob(*quality_score)*(1_f64/3_f64)).ln(),
        }
    }
    LogProb(match_log_likelihood)
}

/// Compute minimum length of match required for a partial match.
pub fn kmer_length(seq_len: usize, percent_mismatch: EMProb)->usize{
    let num_mismatches: usize = (seq_len as EMProb * (percent_mismatch/100 as EMProb)).floor() as usize;
    seq_len/(num_mismatches+1)
}

/// Count the number of positions at which `read_seq` and `ref_seq` differ.
pub fn num_mismatches(read_seq: &[u8], ref_seq: &[u8])->usize{
    read_seq.iter().zip(ref_seq.iter()).map(|(a,b)| (a!=b) as usize).sum()
}

/// Write alignment matches to a TSV file or stdout.
///
/// If `outpath` is `Some`, the results are written to that path (overwriting any existing file).
/// If `outpath` is `None`, results are printed to stdout.
/// Each row contains: Read_ID, Ref_ID, hit position, and log-likelihood.
pub fn write_matches(outpath: Option<&String>, matches: &HashMap<String, Vec<(String, usize, f32)>>)->std::io::Result<()>
{
    match outpath{
        Some(path) => {
            if Path::new(path).exists(){
                fs::remove_file(path)?;
            }
            let mut file_ref = fs::OpenOptions::new().create_new(true).append(true).open(path).expect("Invalid path for result file!");
            let result_header:String = "Read_ID\tRef_ID\thit position\tLog-Likelihood\n".to_string();
            file_ref.write_all(result_header.as_bytes()).expect("write failed");
            matches.iter()
                .filter(|(_seq_id, all_matches)| {
                    all_matches[0].0.is_empty()
                })
                .for_each(|(seq_id, all_matches)| {
                    all_matches.iter()
                        .for_each(|(read_id, hit_pos, match_score)| {
                            let out_string:String = format!("{}\t{}\t{}\t{}\n", seq_id, read_id, hit_pos, match_score);
                            file_ref.write_all(out_string.as_bytes()).expect("write failed");
                        });
                    
            });
        },
        None => {
            let result_header:String = "Read_ID\tRef_ID\tposition\tLog-Likelihood\n".to_string();
            print!("{}", result_header);
            matches.iter()
                .for_each(|(seq_id, all_matches)| {
                    all_matches.iter()
                        .for_each(|(read_id, hit_pos, match_score)| {
                            let out_string:String = format!("{}\t{}\t{}\t{}\n", seq_id, read_id, hit_pos, match_score);
                            println!("{}", out_string);
                        });
        });
        }
    }
    Ok(())
}

/// Convert a Phred+33 quality byte to its linear-space base-call error probability.
///
/// Uses the standard formula: P(error) = 10^(-(Q - 33) / 10).
pub fn error_prob(q: u8)->Prob{
    Prob(10_f64.powf(-((q-33) as f64/10_f64)))
}

/// Return the DNA complement of a sequence (A↔T, C↔G).
///
/// Non-ACGT characters are mapped to `'E'` to signal an unexpected base.
/// Note: this does not reverse the sequence; for reverse-complement use
/// [`bio::alphabets::dna::revcomp`] instead.
pub fn complement(q_seq: Vec<char>)->Vec<char>{
    q_seq.iter()
        .map(|x|{
            match x{
                'T' => 'A',
                'C' => 'G',
                'A' => 'T',
                'G' => 'C',
                _ => 'E',
            }
        })
        .collect()
}

/// Extract the file extension from a filename, returning `None` if there is none.
pub fn get_extension_from_filename(filename: &str) -> Option<&str> {
    Path::new(filename)
        .extension()
        .and_then(OsStr::to_str)
}