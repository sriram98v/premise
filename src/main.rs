extern crate clap;
pub mod utils;
use bio::alphabets;
use bio::alphabets::Alphabet;
use bio::bio_types::sequence::SequenceRead;
use bio::io::fasta;
use bio::stats::{LogProb, Prob};
use clap::{Arg, ArgAction, Command, arg};
use dashmap::{DashMap, DashSet};
use indicatif::{ProgressBar, ProgressDrawTarget};
use itertools::Itertools;
use num::{Float, Zero};
use utils::*;
use core::f64;
use std::collections::{BTreeSet, HashSet, VecDeque};
use std::fmt::Debug;
use std::thread;
use std::{collections::HashMap, fs::File, io::BufReader};
use bio::io::fastq;
use indicatif::ProgressStyle;
use rayon::prelude::*;
use anyhow::Result;
use genedex::{FmIndexConfig, alphabet, FmIndexFlat64};
use flate2::read::GzDecoder;
use chrono::Local;
use savefile::prelude::*;
use std::time::Instant;
use std::io::Cursor;
use tiny_http::{Server, Response, Header, StatusCode};
use genedex::text_with_rank_support::{FlatTextWithRankSupport, Block64};
use std::cmp;

/// Newtype wrapper around a raw read identifier string.
///
/// Used as a map key throughout the alignment and EM pipeline. `Deref`s to
/// `String` for convenient string operations.
#[derive(Debug, Clone, Hash, PartialEq, Eq, Savefile, PartialOrd, Ord)]
pub struct ReadID(String);

impl std::ops::Deref for ReadID {
    type Target = String;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}


/// Newtype wrapper around a raw reference sequence identifier string.
///
/// Stored in the FM-index alongside sequence data and used to label output rows.
/// `Deref`s to `String` for convenient string operations.
#[derive(Debug, Clone, Hash, PartialEq, Eq, Savefile)]
pub struct RefID(String);

impl std::ops::Deref for RefID {
    type Target = String;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

/// Dense integer index assigned to each read for use as a `SparseArray` row key.
///
/// Indices are assigned in iteration order when reads are loaded and are stable
/// for the lifetime of a single run. `Deref`s to `usize`.
#[derive(Debug, Clone, Hash, PartialEq, Eq, Copy, Savefile)]
pub struct ReadIdx(usize);

impl std::ops::Deref for ReadIdx {
    type Target = usize;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}


/// Dense integer index assigned to each reference sequence for use as a `SparseArray` column key.
///
/// Indices correspond to the order in which sequences appear in the FASTA file
/// and are stored in the serialized FM-index. `Deref`s to `usize`.
#[derive(Debug, Clone, Hash, PartialEq, Eq, Copy, Savefile)]
pub struct RefIdx(usize);

impl std::ops::Deref for RefIdx {
    type Target = usize;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

/// Serializable container for a built FM-index and its associated metadata.
///
/// Bundles the raw [`FmIndexFlat64`] with bidirectional index↔ID mappings and
/// the original reference sequences, so that a single `.fmidx` file contains
/// everything needed for alignment without re-reading the FASTA.
#[derive(Savefile)]
pub struct IOFMIndex{
    fmidx: FmIndexFlat64<i64>,
    idx_to_id: HashMap<RefIdx, RefID>,
    id_to_idx: HashMap<RefID, RefIdx>,
    idx_to_seq: HashMap<RefIdx, Vec<u8>>
}

/// Floating-point type used for all EM probabilities and log-probabilities in linear space.
///
/// Values in alignment likelihood maps are stored as raw `f64`; the [`bio::stats::LogProb`]
/// newtype is used at boundaries where log-space arithmetic is required.
type EMProb = f64;

/// Per-read alignment likelihoods: maps each reference index to a map of
/// alignment start positions → log-probability of the read originating from
/// that position on that reference.
pub struct MatchLikelihoods(HashMap<RefIdx, HashMap<usize, LogProb>>);

impl MatchLikelihoods {
    /// Create an empty `MatchLikelihoods` map.
    pub fn new() -> Self {
        Self(HashMap::new())
    }
}

impl std::ops::Deref for MatchLikelihoods {
    type Target = HashMap<RefIdx, HashMap<usize, LogProb>>;
    fn deref(&self) -> &Self::Target { &self.0 }
}

impl std::ops::DerefMut for MatchLikelihoods {
    fn deref_mut(&mut self) -> &mut Self::Target { &mut self.0 }
}

/// Collection of alignment likelihoods for all reads in a dataset.
///
/// Maps each read ID (borrowed from the input slice) to its per-reference
/// likelihood deques. The inner [`DashMap`] allows concurrent writes during
/// parallel alignment.
pub struct ReadAlignments<'a>(DashMap<&'a ReadID, HashMap<RefIdx, VecDeque<LogProb>>>);

impl<'a> ReadAlignments<'a> {
    /// Create an empty `ReadAlignments` map.
    pub fn new() -> Self {
        Self(DashMap::new())
    }
}

impl<'a> std::ops::Deref for ReadAlignments<'a> {
    type Target = DashMap<&'a ReadID, HashMap<RefIdx, VecDeque<LogProb>>>;
    fn deref(&self) -> &Self::Target { &self.0 }
}

impl<'a> std::ops::DerefMut for ReadAlignments<'a> {
    fn deref_mut(&mut self) -> &mut Self::Target { &mut self.0 }
}

/// Half-open interval describing where a single MEM seed aligns on both the read and a reference.
///
/// All positions are 0-based and half-open (`start` inclusive, `end` exclusive).
/// Multiple `MEMPos` values for the same read–reference pair are merged into
/// longer alignments by [`merge_kmer_matches`].
#[derive(Debug, Clone, Copy, Hash, PartialEq, Eq, Savefile)]
pub struct MEMPos{
    pub ref_start: usize,
    pub ref_end: usize,
    pub read_start: usize,
    pub read_end: usize,
}

impl MEMPos{
    /// Returns `true` if `kmer` starts at or after this MEM on both the read and reference,
    /// indicating the two seeds could belong to the same collinear chain.
    fn is_consecutive(&self, kmer: &MEMPos)->bool{
        if kmer.read_start>=self.read_start && kmer.ref_start>=self.ref_start{
            return true;
        }
        false
    }
    /// Returns `true` if `kmer` overlaps this MEM with the same offset on both axes,
    /// meaning the two seeds represent the same underlying alignment and can be merged.
    fn is_overlapping(&self, kmer: &MEMPos)->bool{
        if kmer.read_start-self.read_start == kmer.ref_start-self.ref_start && kmer.read_end-self.read_end == kmer.ref_end-self.ref_end{
            return true;
        }
        false
    }
    /// Extend this MEM to span the union of its coordinates and those of `kmer`.
    fn merge(&mut self, kmer: &MEMPos){
        self.read_start = cmp::min(kmer.read_start, self.read_start);
        self.read_end = cmp::max(kmer.read_end, self.read_end);
        self.ref_start = cmp::min(kmer.ref_start, self.ref_start);
        self.ref_end = cmp::max(kmer.ref_end, self.ref_end);
    }
}

/// A matched pair of forward (R1) and reverse (R2) FASTQ records sharing a common read ID.
///
/// Created by joining the R1 and R2 FASTQ files on read ID before alignment,
/// and passed as a slice to [`process_read_pairs`].
pub struct ReadPair{
    read_id: ReadID,
    r1: fastq::Record,
    r2: fastq::Record,
}

/// Dictionary of Keys (DoK) format for sparse array to store likelihoods
#[derive(Debug, Clone, Default)]
pub struct SparseArray<T: Float + Zero + Copy + Send + Sync + Debug>{
    pub values: DashMap<(ReadIdx,RefIdx), T>,
    pub read_idxs_ref_map: DashMap<ReadIdx, HashSet<RefIdx>>,
    pub ref_idxs_read_map: DashMap<RefIdx, HashSet<ReadIdx>>,
}

impl<T: Float + Zero + Copy + Send + Sync + Debug> SparseArray<T>{
    /// Insert a likelihood value for the (read, reference) pair, updating both index maps.
    fn insert(&mut self, read_idx: ReadIdx, ref_idx: RefIdx, val: T){
        self.values.insert((read_idx, ref_idx), val);
        self.read_idxs_ref_map.entry(read_idx).or_default().insert(ref_idx);
        self.ref_idxs_read_map.entry(ref_idx).or_default().insert(read_idx);
    }

    /// Return the stored value for a (read, reference) pair, or zero if absent.
    fn get(&self, index: &(ReadIdx, RefIdx)) -> T {
        if self.values.contains_key(index){
            *self.values.get(index).unwrap()
        }
        else {
            T::zero()
        }
    }

    /// Return the set of all read indices that have at least one stored entry.
    fn get_read_idxs(&self)->HashSet<ReadIdx>{
        self.read_idxs_ref_map.iter().map(|x| *x.key()).collect()
    }

    /// Return the set of all reference indices that have at least one stored entry.
    fn get_ref_idxs(&self)->HashSet<RefIdx>{
        self.ref_idxs_read_map.iter().map(|x| *x.key()).collect()
    }

    /// Return all (reference index, value) pairs stored for a given read.
    fn get_all_read_hits_idx(&self, read_idx: &ReadIdx)->Vec<(RefIdx, T)>{
        let ref_idxs: &HashSet<RefIdx> = &self.read_idxs_ref_map.get(read_idx).unwrap();
        let mut out: Vec<(RefIdx, T)> = Vec::with_capacity(ref_idxs.len());

        for ref_idx in ref_idxs{
            if let Some(x) = self.values.get(&(*read_idx, *ref_idx)) {
                out.push((*ref_idx, *x));
            }
        }
        out
    }

    /// Return all stored values for a given reference (one per read that aligns to it).
    fn get_all_ref_hits(&self, ref_idx: &RefIdx)->Vec<T>{
        let read_idxs: &HashSet<ReadIdx> = &self.ref_idxs_read_map.get(ref_idx).unwrap();
        let mut out: Vec<T> = Vec::with_capacity(read_idxs.len());


        for read_idx in read_idxs{
            if let Some(x) = self.values.get(&(*read_idx, *ref_idx)) {
                out.push(*x);
            }
        }
        out
    }

    /// Return the total number of distinct reads with stored entries.
    fn num_reads(&self)->usize{
        self.get_read_idxs().len()
    }
}

/// Merge overlapping or consecutive MEM seeds into a single, extended MEM.
///
/// Iterates the sorted deque and attempts to absorb each seed into a running
/// MEM if it is both consecutive and overlapping with it; otherwise starts a
/// new running MEM.
fn merge_kmer_matches(kmer_matches: &VecDeque<MEMPos>)->Vec<MEMPos>{
    let mut out_vec = Vec::with_capacity(kmer_matches.len());
    out_vec.push(kmer_matches[0]);
    for kmer_match in kmer_matches{
        for running_mem in out_vec.iter_mut(){
            if running_mem.is_consecutive(kmer_match) && running_mem.is_overlapping(kmer_match){
                    running_mem.merge(kmer_match);
            }
        }
    }

    out_vec
}

/// Filters the matches found for different kmers and removes repeated alignments.
fn clean_kmer_matches(fmidx: &FmIndexFlat64<i64>, refs: &HashMap<RefIdx, Vec<u8>>, record: &fastq::Record, percent_mismatch: &EMProb, complement: bool)->(HashMap<RefIdx, HashSet<usize>>, HashMap<RefIdx, Vec<MEMPos>>){
    let mut match_positions: HashMap<RefIdx, HashSet<usize>> = HashMap::new();
    let mut mems: HashMap<RefIdx, VecDeque<MEMPos>> = HashMap::new(); // Contains all the MEMs between a read and a reference
    let read_seq = match complement {
        true => bio::alphabets::dna::revcomp(record.seq()),
        false => record.seq().to_vec(),
    };
    let kmer_size = kmer_length(read_seq.len(), *percent_mismatch);

    fmidx.locate_many(
        read_seq.windows(kmer_size).filter(|kmer| {
                let read_alphabet = Alphabet::new(*kmer);
                let dna_alphabet = alphabets::dna::alphabet();

                read_alphabet.intersection(&dna_alphabet)==read_alphabet
            })
        )
        .zip(record.seq().windows(kmer_size))
        .enumerate()
        .for_each(|(kmer_start_read, (hits, _read_kmer))| {

            hits.into_iter()
                .for_each(|ref_entry| {
                    let seq_idx: RefIdx = RefIdx(ref_entry.text_id);
                    let _seq = refs.get(&seq_idx).unwrap();
                    let ref_start_pos_kmer = ref_entry.position;

                    mems.entry(seq_idx)
                        .or_default()
                        .push_back(
                                MEMPos { ref_start: ref_start_pos_kmer, 
                                    ref_end: ref_start_pos_kmer+kmer_size,
                                    read_start: kmer_start_read, 
                                    read_end: kmer_start_read+kmer_size }
                        );

                    if ref_start_pos_kmer>=kmer_start_read{
                        // start position of alignment in reference for read
                        let align_start = ref_start_pos_kmer-kmer_start_read;
                        match_positions.entry(seq_idx).and_modify(|align_pos| {align_pos.insert(align_start);}).or_default().insert(align_start);
                    }


                })
        });
    
    let merged_mems: HashMap<RefIdx, Vec<MEMPos>> = mems.into_iter().map(|(k, v)| (k, merge_kmer_matches(&v))).collect();

    (match_positions, merged_mems)
}

/// Aligns a single read to each of the references
/// Returns a pair of Hashmaps. The first maps the read to its best alignment to each reference (reference_name, (alignment_start_pos, likelihood of alignment)).
/// The second returns the sum of likelihoods of all alignments to each reference.(reference_name, (sum of likelihood of all alignments)). 
fn query_read(fmidx: &FmIndexFlat64<i64>, refs: &HashMap<RefIdx, Vec<u8>>, record: &fastq::Record, percent_mismatch: &EMProb, complement: bool)->Result<MatchLikelihoods>{

    let read_len = record.seq().len();
    let read_seq = match complement {
        true => bio::alphabets::dna::revcomp(record.seq()),
        false => record.seq().to_vec(),
    };
    let read_qual = match complement {
        true => record.qual().iter().rev().cloned().collect_vec(),
        false => record.qual().to_vec(),
    };

    let mut match_likelihood = MatchLikelihoods::new();

    let (_other_matches, mems) = clean_kmer_matches(fmidx, refs, record, percent_mismatch, complement);
    
    mems.into_iter().for_each(|hit| {
        let ref_id = hit.0;
        let ref_seq = refs.get(&ref_id).unwrap();
        let ref_len = ref_seq.len();

        
        for mem in hit.1.iter(){
            let read_pos = mem.read_start;

            if mem.ref_start<read_pos{
                continue;
            }
            let ref_pos = mem.ref_start-mem.read_start;
            if ref_pos+read_len<ref_len{
                let ref_match_seg = &ref_seq[ref_pos..ref_pos+read_len];
                    let match_log_prob = compute_match_log_prob(&read_seq, &read_qual, ref_match_seg);

                    if match_log_prob.exp()!=0.0{
                        // update match score
                        match_likelihood.entry(ref_id)
                            .and_modify(|e| {
                                e.insert(ref_pos, match_log_prob);
                            })
                            .or_default()
                                .insert(ref_pos, match_log_prob);
                    }
            }
        }
    });

    Ok(match_likelihood)
}

/// Compute the combined log-probability for a read pair aligning to the same reference.
///
/// Sums the linear-space product of each (R1, R2_rc) alignment position pair where
/// the R2 reverse-complement ends after the R1 start position (i.e. the pair is
/// in a valid FR orientation), then returns the result in log space.
fn merge_read_pairs(forward: &HashMap<usize, LogProb>, reverse: &HashMap<usize, LogProb>, read_len: usize) -> LogProb{
    let mut match_likelihood: EMProb = 0.0;
    for (r1_start, r1_log_prob) in forward.iter(){
        for (r2_rc_start, r2_log_prob) in reverse.iter(){
            let r2_rc_end = r2_rc_start + read_len;
            if r2_rc_end>*r1_start{
                match_likelihood += (r1_log_prob+r2_log_prob).exp();
            }
        }
    }
    LogProb(match_likelihood.ln())
}

/// Align all read pairs against the FM-index in parallel and collect per-read likelihoods.
///
/// For each `ReadPair`, both orientations (FR and RF) are queried via [`query_read`].
/// Alignment likelihoods that fall below the absolute threshold `eps_1` or are more
/// than `eps_2` log-units below the best alignment for that read are discarded.
/// Returns a map of read ID → per-reference likelihood deques, and the set of all
/// references that received at least one alignment.
fn process_read_pairs<'a>(fmidx: &FmIndexFlat64<i64>,
                    refs: &HashMap<RefIdx, Vec<u8>>,
                    read_pairs: &'a[ReadPair],
                    percent_mismatch: &EMProb,
                    eps_1: LogProb,
                    eps_2: LogProb)-> Result<(ReadAlignments<'a>, DashSet<RefIdx>)>
{
    let out_aligns = ReadAlignments::new();

    let pb = ProgressBar::with_draw_target(Some(read_pairs.len() as u64), ProgressDrawTarget::stderr());
    pb.set_style(ProgressStyle::with_template("Finding pairwise alignments: {spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent}% ({eta})").unwrap());

    let all_ref_ids: DashSet<RefIdx> = DashSet::new();

    let half_log_prob = LogProb::from(Prob(0.5));

    read_pairs
        .par_iter()
        .map(|read_pair| {
            let r1_rec = &read_pair.r1;
            let r2_rec = &read_pair.r2;

            let r1_len = r1_rec.len();
            let r2_len = r2_rec.len();

            let mut match_likelihoods: HashMap<RefIdx, LogProb> = HashMap::new();

            let r1_match_likelihoods = query_read(fmidx, refs, r1_rec, percent_mismatch, false).unwrap();
            let r1_rc_match_likelihoods = query_read(fmidx, refs, r1_rec, percent_mismatch, true).unwrap();

            let r2_match_likelihoods = query_read(fmidx, refs, r2_rec, percent_mismatch, false).unwrap();
            let r2_rc_match_likelihoods = query_read(fmidx, refs, r2_rec, percent_mismatch, true).unwrap();


            r1_match_likelihoods.keys().filter(|k| r2_rc_match_likelihoods.contains_key(k)).for_each(|x| {
                let match_prob = half_log_prob + merge_read_pairs(r1_match_likelihoods.get(x).unwrap(), r2_rc_match_likelihoods.get(x).unwrap(), r2_len);
                match_likelihoods.insert(*x, match_prob);
            });
            
            r1_rc_match_likelihoods.keys().filter(|k| r2_match_likelihoods.contains_key(k)).for_each(|x| {
                let match_prob = half_log_prob + merge_read_pairs(r2_match_likelihoods.get(x).unwrap(), r1_rc_match_likelihoods.get(x).unwrap(), r1_len);
                match_likelihoods.entry(*x).and_modify(|v| *v = LogProb::from(Prob(v.exp() + match_prob.exp()))).or_insert(match_prob);
            });

            let all_keys = match_likelihoods.keys().cloned().collect_vec();
            
            if match_likelihoods.len()>=1{
                let max_likelihood = match_likelihoods.values().max_by(|a, b| a.total_cmp(&b)).unwrap().clone();
                for k in all_keys.iter(){
                    if *match_likelihoods.get(k).unwrap()<eps_1 || match_likelihoods.get(k).unwrap()-max_likelihood<=eps_2{
                        match_likelihoods.remove(k);
                    }
                }
            }

            (&read_pair.read_id, match_likelihoods)
        })
        .for_each(|(read_id, match_likelihoods)| {
            match_likelihoods.iter()
                .for_each(|x| {
                    all_ref_ids.insert(*x.0);

                    out_aligns.entry(read_id)
                        .or_default()
                        .entry(*x.0).or_default().push_back(*x.1);

                });
            pb.inc(1);
        });

    pb.finish_with_message("");
    
    Ok((out_aligns, all_ref_ids))
}

/// Run the unpenalized EM algorithm to estimate reference proportions.
///
/// Initializes proportions by plurality vote, then iterates the E- and M-steps
/// for up to `num_iter` rounds, stopping early when the change in data
/// log-likelihood drops below 1e-6. Returns the per-read MAP assignments,
/// the filtered proportion map, the final posterior weight matrix, and the
/// data log-likelihood at each iteration.
fn get_proportions_par_sparse(ll_array: &SparseArray<EMProb>, num_iter: usize)->(HashMap<ReadIdx, RefIdx>, HashMap<RefIdx, EMProb>, SparseArray<EMProb>, Vec<EMProb>){
    let num_reads = ll_array.num_reads();

    let read_idxs = ll_array.get_read_idxs();
    let _ref_idxs = ll_array.get_ref_idxs();

    let mut props: Vec<HashMap<RefIdx, EMProb>> = vec![HashMap::new();num_iter+1];
    let w: SparseArray<EMProb> = SparseArray::default();

    let initial_props: DashMap<RefIdx, usize> = DashMap::new();

    read_idxs.par_iter().for_each(|read_idx| {
        let row_argmax = ll_array.get_all_read_hits_idx(read_idx).iter().max_by(|&(_, f1), &(_, f2)| {
            EMProb::total_cmp(f1, f2)
        }).unwrap().0;
        initial_props.entry(row_argmax).and_modify(|e| { *e += 1 }).or_insert(1);
    });

    let total = initial_props.iter().map(|x| *x.value()).sum::<usize>();

    for (ref_idx, count) in initial_props.into_iter(){
        props[0].insert(ref_idx, (count as EMProb)/(total as EMProb));
    }

    let pb = ProgressBar::with_draw_target(Some(num_iter as u64), ProgressDrawTarget::stderr());
    pb.set_style(ProgressStyle::with_template("Running EM: {spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent}% ({eta}) ({msg})").unwrap());

    let clik: DashMap<ReadIdx, EMProb> = DashMap::new();

    read_idxs
        .iter()
        .par_bridge()
        .for_each(|read_idx| {
            let pr_rspr_s = ll_array.get_all_read_hits_idx(read_idx)
                .into_iter()
                .map(|(ref_idx, val)| val*props[0].get(&ref_idx).unwrap_or(&(0 as EMProb))).sum::<EMProb>();

            clik.entry(*read_idx).and_modify(|e| *e += pr_rspr_s).or_insert(pr_rspr_s);
        });

    
    let mut prev_data_loglikelihood = clik.iter().filter(|x| *x.value()!=0.0).map(|x| x.value().ln()).sum::<EMProb>();

    let mut data_likelihoods: Vec<EMProb> = Vec::new();

    data_likelihoods.push(prev_data_loglikelihood);


    for i in 0..num_iter {

        read_idxs.iter()
            .par_bridge()
            .for_each(|x| {
                for (ref_idx, val) in ll_array.get_all_read_hits_idx(x){
                    w.values.insert((*x, ref_idx), val*props[i].get(&ref_idx).unwrap_or(&(0 as EMProb)));
                    w.read_idxs_ref_map.entry(*x).or_default().insert(ref_idx);
                    w.ref_idxs_read_map.entry(ref_idx).or_default().insert(*x);
                }
            });

        let clik: DashMap<ReadIdx, EMProb> = DashMap::new();
        
        read_idxs
            .iter()
            .par_bridge()
            .for_each(|read_idx| {
                let pr_rspr_s = ll_array.get_all_read_hits_idx(read_idx)
                    .into_iter()
                    .map(|(ref_idx, val)| val*props[i].get(&ref_idx).unwrap_or(&(0 as EMProb))).sum::<EMProb>();

                clik.entry(*read_idx).and_modify(|e| *e += pr_rspr_s).or_insert(pr_rspr_s);
            });

        read_idxs.iter()
            .par_bridge()
            .for_each(|x| {
                for (ref_idx, val) in ll_array.get_all_read_hits_idx(x){
                    let new_val = (val*props[i].get(&ref_idx).unwrap_or(&(0 as EMProb)))/ *clik.get(x).unwrap();
                    match new_val.is_finite(){
                        true => w.values.insert((*x, ref_idx),new_val),
                        _ => w.values.insert((*x, ref_idx),0.0),
                    };
                }
            });

        let ejs: HashMap<RefIdx, EMProb> = ll_array.get_ref_idxs()
            .par_iter()
            .map(|ref_idx| {
                let vals = w.get_all_ref_hits(ref_idx);
                (*ref_idx, vals.iter().sum::<EMProb>()) 
            }).collect();

        let lambda_init = ejs.values()
            .map(|x| x-20.0).max_by(|f1, f2| {
            EMProb::total_cmp(f1, f2)
        }).unwrap();

        let lambda = _update_lambda(20.0, 1e-20, &ejs, lambda_init, num_iter);

        props[i+1] = ll_array.get_ref_idxs()
            .par_iter()
            .map(|ref_idx| {
                let tmp_pi = _update_pi(20.0, 1e-20, *ejs.get(ref_idx).unwrap(), lambda);
                if tmp_pi>0.0{
                    (*ref_idx, tmp_pi)
                }
                else{
                    (*ref_idx, 0.0)
                }
            })
            .collect();

        let data_loglikelihood = clik.iter().par_bridge().filter(|x| *x.value()!=0.0).map(|x| x.value().ln()).sum::<EMProb>();

        let data_loglikelihood_diff= data_loglikelihood-prev_data_loglikelihood;

        prev_data_loglikelihood = data_loglikelihood;

        data_likelihoods.push(prev_data_loglikelihood);

        pb.set_message(format!("{data_loglikelihood_diff:.3e}"));

        if data_loglikelihood_diff>0.0 && data_loglikelihood_diff.abs()<=1e-6{
            props[num_iter] = props[i+1].clone();
            break;
        }


        pb.inc(1);

    }
    pb.finish_with_message(format!("Final data LL: {prev_data_loglikelihood}"));

    let results: HashMap<ReadIdx, RefIdx> = read_idxs.iter().par_bridge().map(|read_idx| {
        let row_argmax = w.get_all_read_hits_idx(read_idx).iter().max_by(|&(_, f1), &(_, f2)| {
            EMProb::total_cmp(f1, f2)
        }).unwrap().0;
        (*read_idx, row_argmax)
    }).collect();

    (results, props[num_iter].iter().filter(|(_,v)| **v*(num_reads as EMProb)>1.0).map(|(k,v)|(*k,*v)).collect(), w, data_likelihoods)
}

/// M-step update for a single reference proportion under the L1-regularized objective.
///
/// Computes the closed-form solution for π_j given the current Lagrange multiplier
/// `lambda`, the expected count `ej`, penalty weight `rho`, and floor `omega`.
fn _update_pi(rho: EMProb, omega: EMProb, ej: EMProb, lambda: EMProb)-> EMProb{
    let phi = _compute_phi(lambda, omega, rho, ej);
    (-phi + (phi*phi - 4.0*lambda*ej*omega).sqrt())/(2.0*lambda)
}

/// Find the Lagrange multiplier λ that enforces the simplex constraint Σπ_j = 1.
///
/// Uses Newton-Raphson starting from `lambda_init`, iterating for up to
/// `iterations` steps or until the constraint function equals zero.
fn _update_lambda(rho: EMProb, omega: EMProb, ejs: &HashMap<RefIdx, EMProb>, lambda_init: EMProb, iterations: usize)-> EMProb{
    let mut lambda = lambda_init;
    for _ in 0..iterations{
        // dbg!(lambda, _compute_f(lambda, omega, rho, ejs), _compute_deriv_f(lambda, omega, rho, ejs));
        lambda -= (_compute_f(lambda, omega, rho, ejs))/(_compute_deriv_f(lambda, omega, rho, ejs)); 
        if _compute_f(lambda, omega, rho, ejs)==0.0{
            break;
        }
    };
    // dbg!(lambda, _compute_f(lambda, omega, rho, ejs), _compute_deriv_f(lambda, omega, rho, ejs));
    lambda
}

/// Evaluate the constraint function f(λ) = Σ_j π_j(λ) − 1 used by Newton-Raphson.
///
/// A root of this function gives the λ for which the estimated proportions sum to one.
fn _compute_f(lambda: EMProb, omega: EMProb, rho: EMProb, ejs: &HashMap<RefIdx, EMProb>)-> EMProb{
    ejs.keys()
        .map(|ref_idx| {
            let ej = *ejs.get(ref_idx).unwrap();
            let phi = _compute_phi(lambda, omega, rho, ej);
            -phi + (phi*phi - 4.0*lambda*ej*omega).sqrt()
        }).sum::<EMProb>() - 2.0*lambda
}

/// Evaluate the derivative f′(λ) = Σ_j ∂π_j/∂λ used by Newton-Raphson.
fn _compute_deriv_f(lambda: EMProb, omega: EMProb, rho: EMProb, ejs: &HashMap<RefIdx, EMProb>)-> EMProb{
    ejs.keys()
        .map(|ref_idx| {
            let ej = *ejs.get(ref_idx).unwrap();
            let phi = _compute_phi(lambda, omega, rho, ej);
            omega + 0.5*(1.0/(phi*phi - 4.0*lambda*ej*omega).sqrt())*(2.0*phi*omega-4.0*ej*omega)
        }).sum::<EMProb>() - 2.0
}

/// Compute the auxiliary scalar φ = λω + ρ − e_j used in the closed-form π update.
fn _compute_phi(lambda: EMProb, omega: EMProb, rho: EMProb, ej: EMProb)-> EMProb{
    return lambda*omega+rho-ej
}

/// Return the L1 norm (sum of absolute values) of a proportion map.
fn _l1_norm(props: &HashMap<RefIdx, EMProb>)->EMProb{
    props.values().sum()
}

/// Return the squared L2 norm (sum of squares) of a proportion map.
fn _l2_norm_sq(props: &HashMap<RefIdx, EMProb>)->EMProb{
    props.values().map(|x| x*x).sum()
}

/// Compute the log-barrier penalty Σ_j ln(1 + π_j / ω) for a proportion map.
fn _penalty(props: &HashMap<RefIdx, EMProb>, omega: EMProb)->EMProb{
    props.values().map(|x| (1.0+(x/omega)).ln()).sum()
}

/// Run the L1-penalized EM algorithm to estimate reference proportions.
///
/// Identical in structure to [`get_proportions_par_sparse`] but applies a
/// sparsity-inducing penalty with weight `rho` and floor `omega`. The Lagrange
/// multiplier enforcing the simplex constraint is solved via Newton-Raphson at
/// each M-step. Convergence is declared when the improvement in data
/// log-likelihood is positive but ≤ `em_threshold` (after a 20-iteration burn-in).
fn get_proportions_par_sparse_l1_reg(ll_array: &SparseArray<EMProb>, num_iter: usize, rho: EMProb, omega: EMProb, em_threshold: EMProb)->(HashMap<ReadIdx, RefIdx>, HashMap<RefIdx, EMProb>, SparseArray<EMProb>, Vec<EMProb>){
    let num_reads = ll_array.num_reads();

    let read_idxs = ll_array.get_read_idxs();
    let _ref_idxs = ll_array.get_ref_idxs();

    let mut props: Vec<HashMap<RefIdx, EMProb>> = vec![HashMap::new();num_iter+1];
    let w: SparseArray<EMProb> = SparseArray::default();

    let initial_props: DashMap<RefIdx, usize> = DashMap::new();

    read_idxs.par_iter().for_each(|read_idx| {
        let row_argmax = ll_array.get_all_read_hits_idx(read_idx).iter().max_by(|&(_, f1), &(_, f2)| {
            EMProb::total_cmp(f1, f2)
        }).unwrap().0;
        initial_props.entry(row_argmax).and_modify(|e| { *e += 1 }).or_insert(1);
    });

    let total = initial_props.iter().map(|x| *x.value()).sum::<usize>();

    for (ref_idx, count) in initial_props.into_iter(){
        props[0].insert(ref_idx, (count as EMProb)/(total as EMProb));
    }

    let pb = ProgressBar::with_draw_target(Some(num_iter as u64), ProgressDrawTarget::stderr());
    pb.set_style(ProgressStyle::with_template("Running EM: {spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent}% ({eta}) ({msg})").unwrap());

    let clik: DashMap<ReadIdx, EMProb> = DashMap::new();

    read_idxs
        .iter()
        .par_bridge()
        .for_each(|read_idx| {
            let pr_rspr_s = ll_array.get_all_read_hits_idx(read_idx)
                .into_iter()
                .map(|(ref_idx, val)| val*props[0].get(&ref_idx).unwrap_or(&(0 as EMProb))).sum::<EMProb>();

            clik.entry(*read_idx).and_modify(|e| *e += pr_rspr_s).or_insert(pr_rspr_s);
        });

    
    let mut prev_data_loglikelihood = clik.iter().filter(|x| *x.value()!=0.0).map(|x| x.value().ln()).sum::<EMProb>();

    let mut data_likelihoods: Vec<EMProb> = Vec::new();

    data_likelihoods.push(prev_data_loglikelihood);

    for i in 0..num_iter {

        read_idxs.iter()
            .par_bridge()
            .for_each(|x| {
                for (ref_idx, val) in ll_array.get_all_read_hits_idx(x){
                    w.values.insert((*x, ref_idx), val*props[i].get(&ref_idx).unwrap_or(&(0 as EMProb)));
                    w.read_idxs_ref_map.entry(*x).or_default().insert(ref_idx);
                    w.ref_idxs_read_map.entry(ref_idx).or_default().insert(*x);
                }
            });

        let clik: DashMap<ReadIdx, EMProb> = DashMap::new();
        
        read_idxs
            .iter()
            .par_bridge()
            .for_each(|read_idx| {
                let pr_rspr_s = ll_array.get_all_read_hits_idx(read_idx)
                    .into_iter()
                    .map(|(ref_idx, val)| val*props[i].get(&ref_idx).unwrap_or(&(0 as EMProb))).sum::<EMProb>();

                clik.entry(*read_idx).and_modify(|e| *e += pr_rspr_s).or_insert(pr_rspr_s);
            });

        read_idxs.iter()
            .par_bridge()
            .for_each(|x| {
                for (ref_idx, val) in ll_array.get_all_read_hits_idx(x){
                    let new_val = (val*props[i].get(&ref_idx).unwrap_or(&(0 as EMProb)))/ *clik.get(x).unwrap();
                    match new_val.is_finite(){
                        true => w.values.insert((*x, ref_idx),new_val),
                        _ => w.values.insert((*x, ref_idx),0.0),
                    };
                }
            });

        let data_loglikelihood = clik.iter().par_bridge().filter(|x| *x.value()!=0.0).map(|x| x.value().ln()).sum::<EMProb>();

        let data_loglikelihood_diff= data_loglikelihood-prev_data_loglikelihood;

        prev_data_loglikelihood = data_loglikelihood;

        data_likelihoods.push(prev_data_loglikelihood);

        let ejs: HashMap<RefIdx, EMProb> = ll_array.get_ref_idxs()
            .par_iter()
            .map(|ref_idx| {
                let vals = w.get_all_ref_hits(ref_idx);
                (*ref_idx, vals.iter().sum::<EMProb>()) 
            }).collect();

        let lambda_init = ejs.values()
            .map(|x| x-rho).max_by(|f1, f2| {
            EMProb::total_cmp(f1, f2)
        }).unwrap();

        let lambda = _update_lambda(rho, omega, &ejs, lambda_init, num_iter);

        props[i+1] = ll_array.get_ref_idxs()
            .par_iter()
            .map(|ref_idx| {
                let tmp_pi = _update_pi(rho, 1e-20, *ejs.get(ref_idx).unwrap(), lambda);
                if tmp_pi>0.0{
                    (*ref_idx, tmp_pi)
                }
                else{
                    (*ref_idx, 0.0)
                }
            })
            .collect();

        pb.set_message(format!("{data_loglikelihood_diff:.3e}"));

        if i>20 && data_loglikelihood_diff>0.0 && data_loglikelihood_diff<=em_threshold{
            props[num_iter] = props[i+1].clone();
            break;
        }


        pb.inc(1);

    }
    pb.finish_with_message(format!("Final data LL: {prev_data_loglikelihood}"));

    let results: HashMap<ReadIdx, RefIdx> = read_idxs.iter().par_bridge().map(|read_idx| {
        let row_argmax = w.get_all_read_hits_idx(read_idx).iter().max_by(|&(_, f1), &(_, f2)| {
            EMProb::total_cmp(f1, f2)
        }).unwrap().0;
        (*read_idx, row_argmax)
    }).collect();

    (results, props[num_iter].iter().filter(|(_,v)| **v*(num_reads as EMProb)>1.0).map(|(k,v)|(*k,*v)).collect(), w, data_likelihoods)
}



/// Build a serialized FM-index from raw FASTA bytes.
///
/// Parses all records from `fasta_data`, constructs a [`FmIndexFlat64`] over the
/// concatenated sequences using an IUPAC-tolerant DNA alphabet, serializes the
/// index together with its ID↔index mappings into a byte vector, and returns
/// that vector alongside a human-readable log string.
fn build_index_from_bytes(fasta_data: &[u8]) -> Result<(Vec<u8>, String)> {
    let cursor = Cursor::new(fasta_data);
    let reader = BufReader::new(cursor);
    let records = fasta::Reader::new(reader).records();

    let mut ref_ids: HashMap<RefID, RefIdx> = HashMap::new();
    let mut ref_ids_rev: HashMap<RefIdx, RefID> = HashMap::new();
    let mut refs: HashMap<RefIdx, Vec<u8>> = HashMap::new();
    let mut refs_texts: Vec<Vec<u8>> = vec![];

    for (idx, result) in records.enumerate() {
        let record = result.map_err(|e| anyhow::anyhow!("FASTA parse error: {}", e))?;
        ref_ids.insert(RefID(record.id().to_string()), RefIdx(idx));
        ref_ids_rev.insert(RefIdx(idx), RefID(record.id().to_string()));
        refs.insert(RefIdx(idx), record.seq().to_vec());
        refs_texts.push(record.seq().to_vec());
    }

    if refs_texts.is_empty() {
        return Err(anyhow::anyhow!("No sequences found in FASTA file"));
    }

    let log_str = format!("Timestamp: {}\nNum references: {}",
                Local::now().format("%Y-%m-%d %H:%M:%S"),
                ref_ids.len(),
            );
    println!("{}", log_str);

    let dna_alphabet = alphabet::ascii_dna_iupac_as_dna_with_n();
    let fmidx: FmIndexFlat64<_> = FmIndexConfig::<i64, FlatTextWithRankSupport<i64, Block64>>::new()
        .suffix_array_sampling_rate(1)
        .lookup_table_depth(13)
        .construct_index(
            refs_texts.iter().map(|x| str::from_utf8(x).unwrap()).collect_vec(),
            dna_alphabet,
        );

    let io_struct = IOFMIndex {
        fmidx,
        idx_to_id: ref_ids_rev,
        id_to_idx: ref_ids,
        idx_to_seq: refs,
    };

    let bytes = save_to_mem(0, &io_struct).map_err(|e| anyhow::anyhow!("Serialization error: {}", e))?;
    Ok((bytes, log_str))
}

/// HTTP handler for `POST /api/build`.
///
/// Reads the raw FASTA body from the request, delegates to [`build_index_from_bytes`],
/// and responds with the serialized `.fmidx` binary. The log string is included in
/// the `X-Premise-Log` response header (newlines replaced by ` | `).
fn handle_api_build(mut request: tiny_http::Request) {
    let result: Result<(Vec<u8>, String)> = (|| {
        let mut body = Vec::new();
        request.as_reader().read_to_end(&mut body)?;
        if body.is_empty() {
            return Err(anyhow::anyhow!("Empty request body"));
        }
        build_index_from_bytes(&body)
    })();

    match result {
        Ok((data, log_str)) => {
            println!("{}", log_str);
            let log_header_val = log_str.replace('\n', " | ");
            let ct = Header::from_bytes(b"Content-Type", b"application/octet-stream").unwrap();
            let cd = Header::from_bytes(b"Content-Disposition", b"attachment; filename=\"output.fmidx\"").unwrap();
            let lg = Header::from_bytes(b"X-Premise-Log", log_header_val.as_bytes())
                .unwrap_or_else(|_| Header::from_bytes(b"X-Premise-Log", b"").unwrap());
            let _ = request.respond(Response::from_data(data).with_header(ct).with_header(cd).with_header(lg));
        }
        Err(e) => {
            let _ = request.respond(
                Response::from_string(e.to_string()).with_status_code(StatusCode(500)),
            );
        }
    }
}

// ─── Align helpers ────────────────────────────────────────────────────────────

/// Server-side state for a single upload/run session.
///
/// When the GUI uploads files it receives a session ID; subsequent requests use
/// that ID to locate the temporary directory containing the index and reads.
/// `r1_ext` and `r2_ext` record the file extensions of the uploaded reads so
/// that the correct filenames can be reconstructed at run time.
struct AlignSession {
    dir: std::path::PathBuf,
    r1_ext: Option<String>,
    r2_ext: Option<String>,
}

/// Generate a unique session ID based on the current system time in nanoseconds.
fn new_session_id() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    format!("{:x}", SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos())
}

/// Parse the query-string portion of a URL into a key→value map.
///
/// Splits on `?`, then on `&`, then on the first `=` of each pair.
/// Missing values default to an empty string.
fn parse_qs(url: &str) -> HashMap<String, String> {
    url.split('?')
        .nth(1)
        .unwrap_or("")
        .split('&')
        .filter_map(|pair| {
            let mut it = pair.splitn(2, '=');
            let k = it.next()?.to_string();
            let v = it.next().unwrap_or("").to_string();
            Some((k, v))
        })
        .collect()
}

/// Load forward (R1) reads from a FASTQ or gzipped FASTQ file into a read-ID map.
///
/// Strips a trailing `/1` suffix from read IDs (common in paired-end naming
/// conventions) so that R1 and R2 IDs can be matched by bare read name.
fn load_fastq_forward(path: &str) -> Result<HashMap<ReadID, fastq::Record>> {
    match get_extension_from_filename(path) {
        Some("gz") => {
            let f = File::open(path)?;
            Ok(fastq::Reader::from_bufread(BufReader::new(GzDecoder::new(f)))
                .records()
                .filter_map(|x| x.ok())
                .map(|rec| (ReadID(rec.id().to_string()), rec))
                .collect())
        }
        Some("fastq") | Some("fq") => {
            let f = File::open(path)?;
            Ok(fastq::Reader::from_bufread(BufReader::new(f))
                .records()
                .filter_map(|x| x.ok())
                .map(|rec| (ReadID(rec.id().strip_suffix("/1").unwrap_or(rec.id()).to_string()), rec))
                .collect())
        }
        _ => Err(anyhow::anyhow!("Unsupported R1 file type: {}", path)),
    }
}

/// Load reverse (R2) reads from a FASTQ or gzipped FASTQ file into a read-ID map.
///
/// Strips a trailing `/2` suffix from read IDs so that IDs match those in the R1 map.
fn load_fastq_reverse(path: &str) -> Result<HashMap<ReadID, fastq::Record>> {
    match get_extension_from_filename(path) {
        Some("gz") => {
            let f = File::open(path)?;
            Ok(fastq::Reader::from_bufread(BufReader::new(GzDecoder::new(f)))
                .records()
                .filter_map(|x| x.ok())
                .map(|rec| (ReadID(rec.id().to_string()), rec))
                .collect())
        }
        Some("fastq") | Some("fq") => {
            let f = File::open(path)?;
            Ok(fastq::Reader::from_bufread(BufReader::new(f))
                .records()
                .filter_map(|x| x.ok())
                .map(|rec| (ReadID(rec.id().strip_suffix("/2").unwrap_or(rec.id()).to_string()), rec))
                .collect())
        }
        _ => Err(anyhow::anyhow!("Unsupported R2 file type: {}", path)),
    }
}

/// Run pairwise alignment of reads against the reference index and return a TSV.
///
/// Loads the FM-index from `ref_file` and paired reads from `r1_file`/`r2_file`,
/// then calls [`process_read_pairs`] with no absolute likelihood cutoff (ε₁ = −∞)
/// and a relative log-probability cutoff of `eps_2`. Returns a TSV string of
/// (ReadID, RefID, Probability) rows — unaligned reads appear as "unclassified" —
/// and a log string summarising run parameters.
fn run_alignment(
    ref_file: &str,
    r1_file: &str,
    r2_file: &str,
    percent_mismatch: EMProb,
    eps_2: EMProb,
    threads: usize,
) -> Result<(String, String)> {
    let num_threads = if threads == 0 {
        thread::available_parallelism()?.get()
    } else {
        threads
    };
    let _ = rayon::ThreadPoolBuilder::new().num_threads(num_threads).build_global();

    let forward_fastq_records = load_fastq_forward(r1_file)?;
    let reverse_fastq_records = load_fastq_reverse(r2_file)?;

    let all_reads: Vec<ReadPair> = reverse_fastq_records
        .into_iter()
        .filter(|(read_id, _)| forward_fastq_records.contains_key(read_id))
        .map(|(read_id, rev_rec)| {
            let fw_rec = forward_fastq_records.get(&read_id).unwrap();
            ReadPair { read_id, r1: fw_rec.clone(), r2: rev_rec }
        })
        .collect();

    let all_read_ids: BTreeSet<ReadID> = forward_fastq_records.into_keys().collect();
    let num_reads = all_read_ids.len();

    let iofmidx: IOFMIndex = load_file(ref_file, 0)?;
    let fmidx = iofmidx.fmidx;
    let ref_ids_rev = iofmidx.idx_to_id;
    let refs = iofmidx.idx_to_seq;

    let log_str = format!("Timestamp: {}\nNum Threads: {}\nPercent Mismatch: {}\nEps_2: {:e}",
                Local::now().format("%Y-%m-%d %H:%M:%S"),
                num_threads,
                percent_mismatch,
                eps_2.exp(),
            );
    println!("{}", log_str);

    let (out_alignments, _) = process_read_pairs(
        &fmidx, &refs, &all_reads, &percent_mismatch, LogProb(f64::NEG_INFINITY), LogProb(eps_2.ln()),
    )?;

    let mut read_ids: HashMap<&ReadID, ReadIdx> = HashMap::new();
    let mut read_ids_rev: HashMap<ReadIdx, &ReadID> = HashMap::new();
    for (n, entry) in out_alignments.iter().enumerate() {
        read_ids.insert(*entry.key(), ReadIdx(n));
        read_ids_rev.insert(ReadIdx(n), entry.key());
    }


    let mut out = String::from("ReadID\tRefID\tProbability\n");
    for read_id in &all_read_ids {
        if !out_alignments.contains_key(read_id) {
            out.push_str(&format!("{}\tunclassified\t-\n", **read_id));
            continue;
        }
        let read_idx = *read_ids.get(read_id).unwrap();
        let read_aligns = out_alignments.get(read_id).unwrap();
        for (ref_idx, aligns) in read_aligns.iter() {
            let ref_id = ref_ids_rev.get(ref_idx).cloned().unwrap();
            for likelihood in aligns {
                out.push_str(&format!(
                    "{}\t{}\t{:.5e}\n",
                    ***read_ids_rev.get(&read_idx).unwrap(),
                    *ref_id,
                    likelihood.exp(),
                ));
            }
        }
    }
    println!(
        "{} of {} reads could not be classified.",
        num_reads - out_alignments.len(),
        num_reads
    );
    Ok((out, log_str))
}

/// Full PREMISE query pipeline: align reads then run EM classification.
///
/// Loads the FM-index and reads, aligns all pairs via [`process_read_pairs`]
/// using the absolute cutoff `eps_1` and relative cutoff `eps_2`, then runs
/// either the penalized (`rho`, `omega`) or unpenalized EM for up to `num_iter`
/// iterations with convergence criterion `em_threshold`.
///
/// Returns six values:
/// - `matches_tsv`: per-read MAP assignment table (TSV)
/// - `posteriors_tsv`: full posterior probability matrix (TSV)
/// - `props_tsv`: estimated reference abundance proportions (TSV)
/// - `aligns_tsv`: raw per-read alignment likelihoods (TSV, same format as `run_alignment`)
/// - `em_data_likelihoods`: data log-likelihood at each EM iteration
/// - `log_str`: human-readable summary of run parameters
fn run_query(
    ref_file: &str,
    r1_file: &str,
    r2_file: &str,
    percent_mismatch: EMProb,
    eps_1: EMProb,
    eps_2: EMProb,
    num_iter: usize,
    rho: EMProb,
    omega: EMProb,
    em_threshold: EMProb,
    use_penalty: bool,
    threads: usize,
) -> Result<(String, String, String, String, Vec<EMProb>, String)> {
    let num_threads = if threads == 0 {
        thread::available_parallelism()?.get()
    } else {
        threads
    };
    let _ = rayon::ThreadPoolBuilder::new().num_threads(num_threads).build_global();

    let forward_fastq_records = load_fastq_forward(r1_file)?;
    let reverse_fastq_records = load_fastq_reverse(r2_file)?;

    let all_reads: Vec<ReadPair> = reverse_fastq_records
        .into_iter()
        .filter(|(read_id, _)| forward_fastq_records.contains_key(read_id))
        .map(|(read_id, rev_rec)| {
            let fw_rec = forward_fastq_records.get(&read_id).unwrap();
            ReadPair { read_id, r1: fw_rec.clone(), r2: rev_rec }
        })
        .collect();

    let all_read_ids: BTreeSet<ReadID> = forward_fastq_records.into_keys().collect();

    let iofmidx: IOFMIndex = load_file(ref_file, 0)?;
    let fmidx = iofmidx.fmidx;
    let ref_ids_rev = iofmidx.idx_to_id;
    let refs = iofmidx.idx_to_seq;

    let log_str = format!("Timestamp: {}\nNum Threads: {}\nPercent Mismatch: {}\nEps_1: {:e} ({:.2})\nEps_2: {:e} ({:.2})\nEM Iterations: {}\nEM Threshold: {:e}",
        Local::now().format("%Y-%m-%d %H:%M:%S"),
        num_threads,
        percent_mismatch,
        eps_1,
        eps_1.ln(),
        eps_2,
        eps_2.ln(),
        num_iter,
        em_threshold,
    );
    println!("{}", log_str);

    let (out_aligns, _all_refs) = process_read_pairs(
        &fmidx, &refs, &all_reads, &percent_mismatch, LogProb(eps_1.ln()), LogProb(eps_2.ln()),
    )?;

    let mut read_ids: HashMap<&ReadID, ReadIdx> = HashMap::new();
    let mut read_ids_rev: HashMap<ReadIdx, &ReadID> = HashMap::new();
    for (n, val) in out_aligns.iter().enumerate() {
        read_ids.insert(*val.key(), ReadIdx(n));
        read_ids_rev.insert(ReadIdx(n), *val.key());
    }

    let mut ll_array: SparseArray<EMProb> = SparseArray::default();
    for val in out_aligns.iter() {
        let read_id = val.key();
        let alignments = val.value();
        for (ref_idx, positions) in alignments {
            let read_idx = read_ids.get(read_id).unwrap();
            for score in positions.iter() {
                if score.exp().is_finite() && score.exp() != 0.0 {
                    ll_array.insert(*read_idx, *ref_idx, score.exp() as EMProb);
                }
            }
        }
    }

    let (read_assignments, props, posteriors, em_data_likelihoods) = if use_penalty {
        get_proportions_par_sparse_l1_reg(&ll_array, num_iter, rho, omega, em_threshold)
    } else {
        get_proportions_par_sparse(&ll_array, num_iter)
    };

    // Build props TSV (no header — just ref_id \t proportion)
    let mut props_tsv = String::new();
    for (ref_idx, prop) in &props {
        if *prop > 0.0 {
            let ref_id = ref_ids_rev.get(ref_idx).unwrap();
            props_tsv.push_str(&format!("{}\t{:.5e}\n", **ref_id, prop));
        }
    }

    // Build posteriors TSV
    let mut posteriors_tsv = String::from("ReadID\tRefID\tPosterior\n");
    for read_id in all_read_ids.iter() {
        let read_idx = read_ids.get(&read_id);
        if !out_aligns.contains_key(read_id)
            || read_idx.is_none()
            || !read_assignments.contains_key(read_idx.unwrap())
        {
            posteriors_tsv.push_str(&format!("{}\tunclassified\t-\n", **read_id));
            continue;
        }
        let read_idx = *read_ids.get(read_id).unwrap();
        let read_aligns = out_aligns.get(read_id).unwrap();
        for (ref_idx, aligns) in read_aligns.value() {
            let ref_id = ref_ids_rev.get(ref_idx).cloned().unwrap();
            for _ in aligns {
                posteriors_tsv.push_str(&format!(
                    "{}\t{}\t{:.5e}\n",
                    ***read_ids_rev.get(&read_idx).unwrap(),
                    *ref_id.clone(),
                    posteriors.get(&(read_idx, *ref_idx)),
                ));
            }
        }
    }

    // Build matches TSV
    let mut matches_tsv = String::from("ReadID\tRefID\tPosterior\n");
    for read_id in all_read_ids.iter() {
        let read_idx = read_ids.get(&read_id);
        if !out_aligns.contains_key(read_id)
            || read_idx.is_none()
            || !read_assignments.contains_key(read_idx.unwrap())
        {
            matches_tsv.push_str(&format!("{}\tunclassified\t-\n", **read_id));
            continue;
        }
        let ref_idx = read_assignments.get(read_idx.unwrap()).unwrap();
        let ref_id = ref_ids_rev.get(ref_idx).cloned().unwrap();
        let read_id = read_ids_rev.get(read_idx.unwrap()).unwrap();

        if !out_aligns.get(read_id).unwrap().contains_key(ref_idx) {
            continue;
        }

        if out_aligns.get(read_id).unwrap().get(ref_idx).unwrap().front().is_some()
            && ll_array.get(&(*read_idx.unwrap(), *ref_idx)).is_finite()
        {
            matches_tsv.push_str(&format!(
                "{}\t{}\t{:.5e}\n",
                ***read_ids_rev.get(read_idx.unwrap()).unwrap(),
                *ref_id.clone(),
                posteriors.get(&(*read_idx.unwrap(), *ref_idx)),
            ));
        }
    }

    // Build aligns TSV — same format as run_alignment output
    let mut aligns_tsv = String::from("ReadID\tRefID\tProbability\n");
    for read_id in all_read_ids.iter() {
        if !out_aligns.contains_key(read_id) {
            aligns_tsv.push_str(&format!("{}\tunclassified\t-\n", **read_id));
            continue;
        }
        let read_idx = *read_ids.get(read_id).unwrap();
        let read_aligns = out_aligns.get(read_id).unwrap();
        for (ref_idx, aligns) in read_aligns.iter() {
            let ref_id = ref_ids_rev.get(ref_idx).cloned().unwrap();
            for likelihood in aligns {
                aligns_tsv.push_str(&format!(
                    "{}\t{}\t{:.5e}\n",
                    ***read_ids_rev.get(&read_idx).unwrap(),
                    *ref_id,
                    likelihood.exp(),
                ));
            }
        }
    }

    Ok((matches_tsv, posteriors_tsv, props_tsv, aligns_tsv, em_data_likelihoods, log_str))
}

/// Handle a file-upload request for the alignment workflow.
///
/// The `part` query parameter selects which file is being uploaded:
/// `"index"` → FM-index, `"r1"` → forward reads, `"r2"` → reverse reads.
/// A new session directory is created when `session=new` (or the parameter is
/// absent); subsequent parts reuse the existing session. Returns a JSON object
/// with the session ID.
fn do_align_upload(
    request: &mut tiny_http::Request,
    sessions: &mut HashMap<String, AlignSession>,
) -> Result<String> {
    let url = request.url().to_string();
    let qs = parse_qs(&url);
    let part = qs.get("part").map(|s| s.as_str()).unwrap_or("").to_string();
    let ext  = qs.get("ext").cloned().unwrap_or_else(|| "fastq".to_string());

    let session_id = if qs.get("session").map(|s| s.as_str()) == Some("new")
        || !qs.contains_key("session")
    {
        let id = new_session_id();
        let dir = std::env::temp_dir().join(format!("premise_{}", id));
        std::fs::create_dir_all(&dir)?;
        sessions.insert(id.clone(), AlignSession { dir, r1_ext: None, r2_ext: None });
        id
    } else {
        qs.get("session").unwrap().clone()
    };

    let session = sessions
        .get_mut(&session_id)
        .ok_or_else(|| anyhow::anyhow!("Session not found: {}", session_id))?;

    let file_name = match part.as_str() {
        "index" => "index.fmidx".to_string(),
        "r1"    => { session.r1_ext = Some(ext.clone()); format!("r1.{}", ext) }
        "r2"    => { session.r2_ext = Some(ext.clone()); format!("r2.{}", ext) }
        _       => return Err(anyhow::anyhow!("Unknown upload part: {}", part)),
    };

    let file_path = session.dir.join(&file_name);
    let mut out_file = File::create(&file_path)?;
    std::io::copy(request.as_reader(), &mut out_file)?;

    Ok(format!(r#"{{"session":"{}","ok":true}}"#, session_id))
}

/// HTTP handler for `POST /api/align/upload`. Delegates to [`do_align_upload`]
/// and responds with JSON or a 500 error.
fn handle_align_upload(mut request: tiny_http::Request, sessions: &mut HashMap<String, AlignSession>) {
    match do_align_upload(&mut request, sessions) {
        Ok(json) => {
            let ct = Header::from_bytes(b"Content-Type", b"application/json").unwrap();
            let _ = request.respond(Response::from_string(json).with_header(ct));
        }
        Err(e) => {
            let _ = request.respond(
                Response::from_string(e.to_string()).with_status_code(StatusCode(500)),
            );
        }
    }
}

/// Parse query parameters and run pairwise alignment for an uploaded session.
///
/// Reads `mismatch`, `eps_2`, and `threads` from the URL query string, resolves
/// the session's uploaded file paths, calls [`run_alignment`], and returns the
/// result as a JSON object with `tsv` and `log` fields.
fn do_align_run(
    request: &tiny_http::Request,
    sessions: &HashMap<String, AlignSession>,
) -> Result<String> {
    let url = request.url().to_string();
    let qs = parse_qs(&url);
    let session_id = qs.get("session").ok_or_else(|| anyhow::anyhow!("Missing session"))?;
    let mismatch: EMProb = qs.get("mismatch").and_then(|s| s.parse().ok()).unwrap_or(5.0);
    let eps_2: EMProb    = qs.get("eps_2").and_then(|s| s.parse().ok()).unwrap_or(1e-18);
    let threads: usize   = qs.get("threads").and_then(|s| s.parse().ok()).unwrap_or(0);

    let session = sessions
        .get(session_id)
        .ok_or_else(|| anyhow::anyhow!("Session not found: {}", session_id))?;

    let r1_ext = session.r1_ext.as_deref().unwrap_or("fastq");
    let r2_ext = session.r2_ext.as_deref().unwrap_or("fastq");

    let ref_path = session.dir.join("index.fmidx");
    let r1_path  = session.dir.join(format!("r1.{}", r1_ext));
    let r2_path  = session.dir.join(format!("r2.{}", r2_ext));

    let (tsv, log_str) = run_alignment(
        ref_path.to_str().unwrap(),
        r1_path.to_str().unwrap(),
        r2_path.to_str().unwrap(),
        mismatch,
        eps_2,
        threads,
    )?;
    Ok(format!(
        r#"{{"ok":true,"tsv":"{}","log":"{}"}}"#,
        json_escape(&tsv),
        json_escape(&log_str),
    ))
}

/// HTTP handler for `POST /api/align/run`. Delegates to [`do_align_run`]
/// and responds with JSON or a 500 error.
fn handle_align_run(request: tiny_http::Request, sessions: &HashMap<String, AlignSession>) {
    match do_align_run(&request, sessions) {
        Ok(json) => {
            let ct = Header::from_bytes(b"Content-Type", b"application/json").unwrap();
            let _ = request.respond(Response::from_string(json).with_header(ct));
        }
        Err(e) => {
            let _ = request.respond(
                Response::from_string(e.to_string()).with_status_code(StatusCode(500)),
            );
        }
    }
}

/// Escape a string for embedding inside a JSON double-quoted value.
///
/// Replaces `"`, `\`, newline, carriage return, and tab with their JSON escape
/// sequences. Other characters are passed through unchanged.
fn json_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '"'  => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c    => out.push(c),
        }
    }
    out
}

/// Parse query parameters and run the full PREMISE classification pipeline for an uploaded session.
///
/// Reads all EM parameters (`mismatch`, `eps_1`, `eps_2`, `iter`, `rho`, `omega`,
/// `em_threshold`, `no_penalty`, `threads`) from the URL query string, resolves
/// session file paths, calls [`run_query`], and returns a JSON object containing
/// `matches`, `posteriors`, `props`, `convergence` (array of per-iteration data
/// log-likelihoods), and `log` fields.
fn do_query_run(
    request: &tiny_http::Request,
    sessions: &HashMap<String, AlignSession>,
) -> Result<String> {
    let url = request.url().to_string();
    let qs = parse_qs(&url);
    let session_id = qs.get("session").ok_or_else(|| anyhow::anyhow!("Missing session"))?;
    let mismatch: EMProb  = qs.get("mismatch").and_then(|s| s.parse().ok()).unwrap_or(5.0);
    let eps_1: EMProb     = qs.get("eps_1").and_then(|s| s.parse().ok()).unwrap_or(1e-32);
    let eps_2: EMProb     = qs.get("eps_2").and_then(|s| s.parse().ok()).unwrap_or(1e-18);
    let num_iter: usize   = qs.get("iter").and_then(|s| s.parse().ok()).unwrap_or(100);
    let rho: EMProb       = qs.get("rho").and_then(|s| s.parse().ok()).unwrap_or(20.0);
    let omega: EMProb     = qs.get("omega").and_then(|s| s.parse().ok()).unwrap_or(1e-20);
    let em_threshold: EMProb = qs.get("em_threshold").and_then(|s| s.parse().ok()).unwrap_or(1e-6);
    let no_penalty: bool  = qs.get("no_penalty").map(|s| s == "true").unwrap_or(false);
    let use_penalty       = !no_penalty;
    let threads: usize    = qs.get("threads").and_then(|s| s.parse().ok()).unwrap_or(0);

    let session = sessions
        .get(session_id)
        .ok_or_else(|| anyhow::anyhow!("Session not found: {}", session_id))?;

    let r1_ext = session.r1_ext.as_deref().unwrap_or("fastq");
    let r2_ext = session.r2_ext.as_deref().unwrap_or("fastq");

    let ref_path = session.dir.join("index.fmidx");
    let r1_path  = session.dir.join(format!("r1.{}", r1_ext));
    let r2_path  = session.dir.join(format!("r2.{}", r2_ext));

    let (matches_tsv, posteriors_tsv, props_tsv, aligns_tsv, em_likelihoods, log_str) = run_query(
        ref_path.to_str().unwrap(),
        r1_path.to_str().unwrap(),
        r2_path.to_str().unwrap(),
        mismatch, eps_1, eps_2, num_iter, rho, omega, em_threshold, use_penalty, threads,
    )?;

    let convergence_json = format!("[{}]",
        em_likelihoods.iter()
            .map(|v| if v.is_finite() { format!("{:.10e}", v) } else { "null".to_string() })
            .collect::<Vec<_>>()
            .join(",")
    );

    Ok(format!(
        r#"{{"ok":true,"matches":"{}","posteriors":"{}","props":"{}","aligns":"{}","convergence":{},"log":"{}"}}"#,
        json_escape(&matches_tsv),
        json_escape(&posteriors_tsv),
        json_escape(&props_tsv),
        json_escape(&aligns_tsv),
        convergence_json,
        json_escape(&log_str),
    ))
}

/// HTTP handler for `POST /api/query/run`. Delegates to [`do_query_run`]
/// and responds with JSON or a 500 error.
fn handle_query_run(request: tiny_http::Request, sessions: &HashMap<String, AlignSession>) {
    match do_query_run(&request, sessions) {
        Ok(json) => {
            let ct = Header::from_bytes(b"Content-Type", b"application/json").unwrap();
            let _ = request.respond(Response::from_string(json).with_header(ct));
        }
        Err(e) => {
            let _ = request.respond(
                Response::from_string(e.to_string()).with_status_code(StatusCode(500)),
            );
        }
    }
}

fn main() -> Result<()>{
    let matches = Command::new("Maximum Likelihood Metagenomic Classification")
        .version("1.0")
        .author("Sriram Vijendran <vijendran.sriram@gmail.com>")
        .subcommand(
            Command::new("build")
                .about("Build FM-index from reference fasta file")
                .arg(arg!(-s --source <SRC_FILE> "Source file with sequences(fasta)")
                    .required(true)
                )
                .arg(arg!(-o --out <OUTFILE> "Output index file name")
                    .default_value("")
                    .value_parser(clap::value_parser!(String))
                )

        )
        .subcommand(
            Command::new("fasta")
               .about("Generate Fasta file containing sequences present in index")
               .arg(arg!(-i --index <INDEX_FILE> "Source index file of reference sequences(.fmidx)")
                    .required(true)
                    .value_parser(clap::value_parser!(String))
                )
        )
        .subcommand(
            Command::new("inspect")
               .about("Inspect a pre-built index")
               .arg(arg!(-i --index <INDEX_FILE> "Source index file of reference sequences(.fmidx)")
                    .required(true)
                    .value_parser(clap::value_parser!(String))
                )
        )
        .subcommand(
            Command::new("align")
                .arg(arg!(-s --source <SRC_FILE> "Source index file with reference sequences(.fmidx)")
                    .required(true)
                    .value_parser(clap::value_parser!(String))
                    )
                .arg(arg!(-p --percent_mismatch <PERCENT_MISMATCH>"Percent mismatch to reference sequences")
                    .required(true)
                    .value_parser(clap::value_parser!(EMProb))
                    )
                .arg(arg!(-'1' --r1 <READS1>"Source file with forward read sequences(fastq or fastq.gz)")
                    .required(true)
                    .value_parser(clap::value_parser!(String))
                    )
                .arg(arg!(-'2' --r2 <READS2>"Source file with reverse read sequences(fastq or fastq.gz)")
                    .required(true)
                    .value_parser(clap::value_parser!(String))
                    )
                .arg(arg!(-o --out <OUT_FILE>"Output file")
                    .default_value("out.aligns")
                    .value_parser(clap::value_parser!(String))
                    )
                .arg(arg!(--eps_2 <EPS_2>"Minimum match log-probability threshold")
                    .default_value("1e-18")
                    .value_parser(clap::value_parser!(EMProb))
                    )
                .arg(arg!(-t --threads <THREADS>"Number of threads (defaults to 2; 0 uses maximum number of threads)")
                    .default_value("2")
                    .value_parser(clap::value_parser!(usize))
                    )
        )
        .subcommand(
            Command::new("query")
                .arg(arg!(-s --source <SRC_FILE> "Source index file with reference sequences(.fmidx)")
                    .required(true)
                    .value_parser(clap::value_parser!(String))
                    )
                .arg(arg!(-p --percent_mismatch <PERCENT_MISMATCH>"Percent mismatch to reference sequences")
                    .required(true)
                    .value_parser(clap::value_parser!(EMProb))
                    )
                .arg(arg!(--eps_1 <EPS_1>"Cutoff likelihood for dropping alignments")
                    .default_value("1e-32")
                    .value_parser(clap::value_parser!(EMProb))
                    )
                .arg(arg!(-'1' --r1 <READS1>"Source file with forward read sequences(fastq or fastq.gz)")
                    .required(true)
                    .value_parser(clap::value_parser!(String))
                    )
                .arg(arg!(-'2' --r2 <READS2>"Source file with reverse read sequences(fastq or fastq.gz)")
                    .required(true)
                    .value_parser(clap::value_parser!(String))
                    )
                .arg(arg!(-i --iter <ITER>"Number of iterations for EM")
                    .default_value("100")
                    .value_parser(clap::value_parser!(usize))
                    )
                .arg(arg!(--omega <OMEGA>"penalty weight")
                    .default_value("1e-20")
                    .value_parser(clap::value_parser!(EMProb))
                    )
                .arg(arg!(--rho <RHO>"penalty weight")
                    .default_value("20")
                    .value_parser(clap::value_parser!(EMProb))
                    )
                .arg(arg!(-o --out <OUT_FILE>"Output file")
                    .default_value("out")
                    .value_parser(clap::value_parser!(String))
                    )
                .arg(Arg::new("no-penalty")
                    .long("no-penalty")
                    .help("Disable penalty")
                    .required(false)
                    .num_args(0)
                    .action(ArgAction::SetFalse))
                .arg(arg!(--eps_2 <EPS_2>"Minimum match log-probability threshold")
                    .default_value("1e-18")
                    .value_parser(clap::value_parser!(EMProb))
                    )
                .arg(arg!(--em_threshold <EM_THRESHOLD>"EM convergence threshold")
                    .default_value("1e-6")
                    .value_parser(clap::value_parser!(EMProb))
                    )
                .arg(arg!(-t --threads <THREADS>"Number of threads (defaults to 2; 0 uses maximum number of threads)")
                    .default_value("2")
                    .value_parser(clap::value_parser!(usize))
                    )
        )
        .subcommand(
            Command::new("server")
                .about("Start a local HTTP server serving the web interface")
                .arg(Arg::new("port")
                    .long("port")
                    .short('p')
                    .help("Port to listen on")
                    .default_value("8080")
                    .value_parser(clap::value_parser!(u16))
                )
                .arg(Arg::new("ip")
                    .long("ip")
                    .help("IP address to bind to")
                    .default_value("127.0.0.1")
                    .value_parser(clap::value_parser!(String))
                )
        )
        .about("Maximum Likelihood Metagenomic classifier using Suffix trees")
        .get_matches();
  
    match matches.subcommand(){
        Some(("build",  sub_m)) => {
            let src_file = sub_m.get_one::<String>("source").expect("required").as_str();
            let outfile = sub_m.get_one::<String>("out").unwrap().as_str();

            let fasta_data = std::fs::read(src_file)?;
            let (idx_bytes, _) = build_index_from_bytes(&fasta_data)?;

            let out_path = match outfile {
                "" => format!("{}.fmidx", src_file),
                p  => p.to_string(),
            };
            std::fs::write(&out_path, &idx_bytes)?;
            println!("Index written to {}", out_path);

        },
        Some(("align",  sub_m)) => {
            let ref_file = sub_m.get_one::<String>("source").expect("required").as_str();
            let r1_file = sub_m.get_one::<String>("r1").expect("required").as_str();
            let r2_file = sub_m.get_one::<String>("r2").expect("required").as_str();
            let percent_mismatch = *sub_m.get_one::<EMProb>("percent_mismatch").expect("required");
            let eps_2 = *sub_m.get_one::<EMProb>("eps_2").expect("required");
            let outfile = sub_m.get_one::<String>("out").unwrap().as_str();
            let threads = *sub_m.get_one::<usize>("threads").expect("required");

            let now = Instant::now();
            let (tsv, _) = run_alignment(ref_file, r1_file, r2_file, percent_mismatch, eps_2, threads)?;
            std::fs::write(outfile, tsv.as_bytes())?;
            println!("Alignment written to {} ({:.2?})", outfile, now.elapsed());
        },
        Some(("query",  sub_m)) => {
            let ref_file = sub_m.get_one::<String>("source").expect("required").as_str();
            let r1_file = sub_m.get_one::<String>("r1").expect("required").as_str();
            let r2_file = sub_m.get_one::<String>("r2").expect("required").as_str();
            let num_iter = *sub_m.get_one::<usize>("iter").expect("required");
            let percent_mismatch = *sub_m.get_one::<EMProb>("percent_mismatch").expect("required");
            let eps_1 = *sub_m.get_one::<EMProb>("eps_1").expect("required");
            let eps_2 = *sub_m.get_one::<EMProb>("eps_2").expect("required");
            let omega = *sub_m.get_one::<EMProb>("omega").expect("required");
            let rho = *sub_m.get_one::<EMProb>("rho").expect("required");
            let outfile = sub_m.get_one::<String>("out").unwrap().as_str();
            let use_penalty = *sub_m.get_one::<bool>("no-penalty").unwrap();
            let em_threshold = *sub_m.get_one::<EMProb>("em_threshold").expect("required");
            let threads = *sub_m.get_one::<usize>("threads").expect("required");

            let now = Instant::now();
            let (matches_tsv, posteriors_tsv, props_tsv, aligns_tsv, _, _) = run_query(
                ref_file, r1_file, r2_file,
                percent_mismatch, eps_1, eps_2, num_iter, rho, omega, em_threshold, use_penalty, threads,
            )?;
            std::fs::write(format!("{}.matches", outfile), matches_tsv.as_bytes())?;
            std::fs::write(format!("{}.posteriors", outfile), posteriors_tsv.as_bytes())?;
            std::fs::write(format!("{}.props", outfile), props_tsv.as_bytes())?;
            std::fs::write(format!("{}.aligns", outfile), aligns_tsv.as_bytes())?;
            println!("Query written to {}.{{matches,posteriors,props,aligns}} ({:.2?})", outfile, now.elapsed());
        },
        Some(("server", sub_m)) => {
            let port = *sub_m.get_one::<u16>("port").unwrap();
            let ip = sub_m.get_one::<String>("ip").unwrap();
            let addr = format!("{}:{}", ip, port);

            let html = include_str!("templates/index.html")
                .replace("{{CSS}}", include_str!("templates/styles.css"))
                .replace("{{JS}}",  include_str!("templates/app.js"));

            let server = Server::http(&addr)
                .map_err(|e| anyhow::anyhow!("Failed to bind to {}: {}", addr, e))?;

            println!("Premise web interface running at http://{}", addr);
            println!("Press Ctrl+C to stop.");

            let mut sessions: HashMap<String, AlignSession> = HashMap::new();
            let mut query_sessions: HashMap<String, AlignSession> = HashMap::new();

            for request in server.incoming_requests() {
                let url = request.url().to_string();
                let path = url.split('?').next().unwrap_or("/").to_string();
                let method = request.method().clone();

                match (method, path.as_str()) {
                    (tiny_http::Method::Get, "/") => {
                        let ct = Header::from_bytes(b"Content-Type", b"text/html; charset=utf-8").unwrap();
                        let _ = request.respond(Response::from_string(&html).with_header(ct));
                    }
                    (tiny_http::Method::Post, "/api/build") => {
                        handle_api_build(request);
                    }
                    (tiny_http::Method::Post, "/api/align/upload") => {
                        handle_align_upload(request, &mut sessions);
                    }
                    (tiny_http::Method::Post, "/api/align/run") => {
                        handle_align_run(request, &sessions);
                    }
                    (tiny_http::Method::Post, "/api/query/upload") => {
                        handle_align_upload(request, &mut query_sessions);
                    }
                    (tiny_http::Method::Post, "/api/query/run") => {
                        handle_query_run(request, &query_sessions);
                    }
                    _ => {
                        let _ = request.respond(
                            Response::from_string("Not Found").with_status_code(StatusCode(404)),
                        );
                    }
                }
            }
        },
        _ => {
            println!("No option selected! Refer help page (-h flag)");
        }
    }

    Ok(())
}
