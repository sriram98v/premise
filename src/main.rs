extern crate clap;
pub mod utils;
use bio::alphabets;
use bio::alphabets::Alphabet;
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

#[derive(Debug, Clone, Hash, PartialEq, Eq, Savefile, PartialOrd, Ord)]
pub struct ReadID(String);

impl std::ops::Deref for ReadID {
    type Target = String;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}


#[derive(Debug, Clone, Hash, PartialEq, Eq, Savefile)]
pub struct RefID(String);

impl std::ops::Deref for RefID {
    type Target = String;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[derive(Debug, Clone, Hash, PartialEq, Eq, Copy, Savefile)]
pub struct ReadIdx(usize);

impl std::ops::Deref for ReadIdx {
    type Target = usize;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}


#[derive(Debug, Clone, Hash, PartialEq, Eq, Copy, Savefile)]
pub struct RefIdx(usize);

impl std::ops::Deref for RefIdx {
    type Target = usize;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[derive(Savefile)]
pub struct IOFMIndex{
    fmidx: FmIndexFlat64<i64>,
    idx_to_id: HashMap<RefIdx, RefID>,
    id_to_idx: HashMap<RefID, RefIdx>,
    idx_to_seq: HashMap<RefIdx, Vec<u8>>
}

type EMProb = f64;

/// type to store all alignments of a read to references
// type Alignments = HashMap<RefIdx, LogProb>;
/// type to store the best alignments of reads to references
type MatchLikelihoods = HashMap<RefIdx, LogProb>;
/// type to store all alignments of all reads to references
type ReadAlignments<'a> = DashMap<&'a ReadID, HashMap<RefIdx, VecDeque<LogProb>>>;

#[derive(Debug, Clone, Copy, Hash, PartialEq, Eq, Savefile)]
pub struct MEMPos{
    pub ref_start: usize,
    pub ref_end: usize,
    pub read_start: usize, 
    pub read_end: usize,
}

impl MEMPos{
    fn is_consecutive(&self, kmer: &MEMPos)->bool{
        if kmer.read_start>=self.read_start && kmer.ref_start>=self.ref_start{
            return true;
        }
        false
    }
    fn is_overlapping(&self, kmer: &MEMPos)->bool{
        if kmer.read_start-self.read_start == kmer.ref_start-self.ref_start && kmer.read_end-self.read_end == kmer.ref_end-self.ref_end{
            return true;
        }
        false
    }
    fn merge(&mut self, kmer: &MEMPos){
        self.read_start = cmp::min(kmer.read_start, self.read_start);
        self.read_end = cmp::max(kmer.read_end, self.read_end);
        self.ref_start = cmp::min(kmer.ref_start, self.ref_start);
        self.ref_end = cmp::max(kmer.ref_end, self.ref_end);
    }
}

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
    fn insert(&mut self, read_idx: ReadIdx, ref_idx: RefIdx, val: T){
        self.values.insert((read_idx, ref_idx), val);
        self.read_idxs_ref_map.entry(read_idx).or_default().insert(ref_idx);
        self.ref_idxs_read_map.entry(ref_idx).or_default().insert(read_idx);
    }

    fn get(&self, index: &(ReadIdx, RefIdx)) -> T {
        if self.values.contains_key(index){
            *self.values.get(index).unwrap()
        }
        else {
            T::zero()
        }
    }

    fn get_read_idxs(&self)->HashSet<ReadIdx>{
        self.read_idxs_ref_map.iter().map(|x| *x.key()).collect()
    }

    fn get_ref_idxs(&self)->HashSet<RefIdx>{
        self.ref_idxs_read_map.iter().map(|x| *x.key()).collect()
    }

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

    fn num_reads(&self)->usize{
        self.get_read_idxs().len()
    }
}

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
fn query_read(fmidx: &FmIndexFlat64<i64>, refs: &HashMap<RefIdx, Vec<u8>>, record: &fastq::Record, percent_mismatch: &EMProb, complement: bool, tau: &EMProb)->Result<(MatchLikelihoods, HashMap<RefIdx, HashSet<MEMPos>>)>{

    let read_len = record.seq().len();
    let read_seq = match complement {
        true => bio::alphabets::dna::revcomp(record.seq()),
        false => record.seq().to_vec(),
    };
    let read_qual = match complement {
        true => record.qual().iter().rev().cloned().collect_vec(),
        false => record.qual().to_vec(),
    };

    let mut match_likelihood: HashMap<RefIdx, LogProb> = HashMap::new();
    let mut softclipped: HashMap<RefIdx, HashSet<MEMPos>> = HashMap::new();

    let (_other_matches, mems) = clean_kmer_matches(fmidx, refs, record, percent_mismatch, complement);
    
    mems.into_iter().for_each(|hit| {
        let ref_id = hit.0;
        let ref_seq = refs.get(&ref_id).unwrap();
        let ref_len = ref_seq.len();

        
        for mem in hit.1.iter(){
            let read_pos = mem.read_start;

            if mem.ref_start<read_pos{
                softclipped.entry(ref_id)
                    .or_default()
                    .insert(*mem);
                continue;
            }
            let ref_pos = mem.ref_start-mem.read_start;
            if ref_pos+read_len<ref_len{
                let ref_match_seg = &ref_seq[ref_pos..ref_pos+read_len];
                    let match_log_prob = compute_match_log_prob(&read_seq, &read_qual, ref_match_seg);

                    if match_log_prob.exp()!=0.0 && match_log_prob>=LogProb(tau.ln()){
                        // update match score
                        match_likelihood.entry(ref_id)
                            .and_modify(|e| {
                                *e += match_log_prob;
                            })
                            .or_insert(match_log_prob);
                    }
            }
        }
    });

    Ok((match_likelihood, softclipped))
}

fn process_read_pairs<'a>(fmidx: &FmIndexFlat64<i64>,
                    refs: &HashMap<RefIdx, Vec<u8>>,
                    read_pairs: &'a[ReadPair],
                    percent_mismatch: &EMProb,
                    cutoff: LogProb,
                    tau: &EMProb)-> Result<(ReadAlignments<'a>, DashSet<RefIdx>)>
{
    let out_aligns: ReadAlignments = DashMap::new();

    let pb = ProgressBar::with_draw_target(Some(read_pairs.len() as u64), ProgressDrawTarget::stderr());
    pb.set_style(ProgressStyle::with_template("Finding pairwise alignments: {spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent}% ({eta})").unwrap());

    let all_ref_ids: DashSet<RefIdx> = DashSet::new();

    let half_log_prob = LogProb::from(Prob(0.5));

    read_pairs
        .par_iter()
        .map(|read_pair| {
            let r1_rec = &read_pair.r1;
            let r2_rec = &read_pair.r2;

            let mut match_likelihoods: HashMap<RefIdx, LogProb> = HashMap::new();

            let (r1_match_likelihoods, _) = query_read(fmidx, refs, r1_rec, percent_mismatch, false, tau).unwrap();
            let (r1_rc_match_likelihoods, _) = query_read(fmidx, refs, r1_rec, percent_mismatch, true, tau).unwrap();

            let (r2_match_likelihoods, _) = query_read(fmidx, refs, r2_rec, percent_mismatch, false, tau).unwrap();
            let (r2_rc_match_likelihoods, _) = query_read(fmidx, refs, r2_rec, percent_mismatch, true, tau).unwrap();

            r1_match_likelihoods.keys().filter(|k| r2_rc_match_likelihoods.contains_key(k)).for_each(|x| {
                let match_prob = half_log_prob + r1_match_likelihoods.get(x).unwrap() + r2_rc_match_likelihoods.get(x).unwrap();
                match_likelihoods.insert(*x, match_prob);
            });

            r1_rc_match_likelihoods.keys().filter(|k| r2_match_likelihoods.contains_key(k)).for_each(|x| {
                let match_prob = half_log_prob + r1_rc_match_likelihoods.get(x).unwrap() + r2_match_likelihoods.get(x).unwrap();
                match_likelihoods.entry(*x).and_modify(|v| *v = LogProb::from(Prob(v.exp() + match_prob.exp()))).or_insert(match_prob);
            });

            // let mut sorted_alignments = match_likelihoods.values().map(|x| x.exp()).collect_vec();
            // sorted_alignments.sort_by(|a, b| a.total_cmp(&b));
            let all_keys = match_likelihoods.keys().cloned().collect_vec();
            if match_likelihoods.len()>=1{
                let max_likelihood = match_likelihoods.values().max_by(|a, b| a.total_cmp(&b)).unwrap().clone();
                for k in all_keys.iter(){
                    if match_likelihoods.get(k).unwrap()-max_likelihood<=cutoff{
                        // dbg!(k, match_likelihoods.get(k).unwrap(), max_likelihood, LogProb::from(Prob::from(1e-6)));
                        match_likelihoods.remove(k);
                    }
                }
            }

            // if all_keys.len()>match_likelihoods.len(){
            //     dbg!(all_keys.len()-match_likelihoods.len());
            // }

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

fn _update_pi(rho: EMProb, omega: EMProb, ej: EMProb, lambda: EMProb)-> EMProb{
    let phi = _compute_phi(lambda, omega, rho, ej);
    (-phi + (phi*phi - 4.0*lambda*ej*omega).sqrt())/(2.0*lambda)
}

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

fn _compute_f(lambda: EMProb, omega: EMProb, rho: EMProb, ejs: &HashMap<RefIdx, EMProb>)-> EMProb{
    ejs.keys()
        .map(|ref_idx| {
            let ej = *ejs.get(ref_idx).unwrap();
            let phi = _compute_phi(lambda, omega, rho, ej);
            -phi + (phi*phi - 4.0*lambda*ej*omega).sqrt()
        }).sum::<EMProb>() - 2.0*lambda
}

fn _compute_deriv_f(lambda: EMProb, omega: EMProb, rho: EMProb, ejs: &HashMap<RefIdx, EMProb>)-> EMProb{
    ejs.keys()
        .map(|ref_idx| {
            let ej = *ejs.get(ref_idx).unwrap();
            let phi = _compute_phi(lambda, omega, rho, ej);
            omega + 0.5*(1.0/(phi*phi - 4.0*lambda*ej*omega).sqrt())*(2.0*phi*omega-4.0*ej*omega)
        }).sum::<EMProb>() - 2.0
}

fn _compute_phi(lambda: EMProb, omega: EMProb, rho: EMProb, ej: EMProb)-> EMProb{
    return lambda*omega+rho-ej
}

fn _l1_norm(props: &HashMap<RefIdx, EMProb>)->EMProb{
    props.values().sum()
}

fn _l2_norm_sq(props: &HashMap<RefIdx, EMProb>)->EMProb{
    props.values().map(|x| x*x).sum()
}

fn _penalty(props: &HashMap<RefIdx, EMProb>, omega: EMProb)->EMProb{
    props.values().map(|x| (1.0+(x/omega)).ln()).sum()
}

fn get_proportions_par_sparse_l1_reg(ll_array: &SparseArray<EMProb>, num_iter: usize, rho: EMProb, omega: EMProb)->(HashMap<ReadIdx, RefIdx>, HashMap<RefIdx, EMProb>, SparseArray<EMProb>, Vec<EMProb>){
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

        if i>20 && data_loglikelihood_diff>0.0 && data_loglikelihood_diff<=1e-6{
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

struct AlignSession {
    dir: std::path::PathBuf,
    r1_ext: Option<String>,
    r2_ext: Option<String>,
}

fn new_session_id() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    format!("{:x}", SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos())
}

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

fn run_alignment(
    ref_file: &str,
    r1_file: &str,
    r2_file: &str,
    percent_mismatch: EMProb,
    tau: EMProb,
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

    let log_str = format!("Timestamp: {}\nNum Threads: {}\nPercent Mismatch: {}\nTau: {:e}",
                Local::now().format("%Y-%m-%d %H:%M:%S"),
                num_threads,
                percent_mismatch,
                tau,
            );
    println!("{}", log_str);

    let (out_alignments, _) = process_read_pairs(
        &fmidx, &refs, &all_reads, &percent_mismatch, LogProb(f64::NEG_INFINITY), &tau,
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

fn run_query(
    ref_file: &str,
    r1_file: &str,
    r2_file: &str,
    percent_mismatch: EMProb,
    cutoff: EMProb,
    tau: EMProb,
    num_iter: usize,
    lambda: EMProb,
    gamma: EMProb,
    use_penalty: bool,
    threads: usize,
) -> Result<(String, String, String, Vec<EMProb>, String)> {
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

    let log_str = format!("Timestamp: {}\nNum Threads: {}\nPercent Mismatch: {}\nCutoff likelihood: {:e} ({:.2})\nTau: {:e}\nEM Iterations: {}",
        Local::now().format("%Y-%m-%d %H:%M:%S"),
        num_threads,
        percent_mismatch,
        cutoff,
        cutoff.ln(),
        tau,
        num_iter,
    );
    println!("{}", log_str);

    let (out_aligns, _all_refs) = process_read_pairs(
        &fmidx, &refs, &all_reads, &percent_mismatch, LogProb(cutoff.ln()), &tau,
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
        get_proportions_par_sparse_l1_reg(&ll_array, num_iter, lambda, gamma)
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

    Ok((matches_tsv, posteriors_tsv, props_tsv, em_data_likelihoods, log_str))
}

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

fn do_align_run(
    request: &tiny_http::Request,
    sessions: &HashMap<String, AlignSession>,
) -> Result<String> {
    let url = request.url().to_string();
    let qs = parse_qs(&url);
    let session_id = qs.get("session").ok_or_else(|| anyhow::anyhow!("Missing session"))?;
    let mismatch: EMProb = qs.get("mismatch").and_then(|s| s.parse().ok()).unwrap_or(5.0);
    let tau: EMProb      = qs.get("tau").and_then(|s| s.parse().ok()).unwrap_or(1e-18);
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
        tau,
        threads,
    )?;
    Ok(format!(
        r#"{{"ok":true,"tsv":"{}","log":"{}"}}"#,
        json_escape(&tsv),
        json_escape(&log_str),
    ))
}

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

fn do_query_run(
    request: &tiny_http::Request,
    sessions: &HashMap<String, AlignSession>,
) -> Result<String> {
    let url = request.url().to_string();
    let qs = parse_qs(&url);
    let session_id = qs.get("session").ok_or_else(|| anyhow::anyhow!("Missing session"))?;
    let mismatch: EMProb  = qs.get("mismatch").and_then(|s| s.parse().ok()).unwrap_or(5.0);
    let cutoff: EMProb    = qs.get("cutoff").and_then(|s| s.parse().ok()).unwrap_or(1e-4);
    let tau: EMProb       = qs.get("tau").and_then(|s| s.parse().ok()).unwrap_or(1e-18);
    let num_iter: usize   = qs.get("iter").and_then(|s| s.parse().ok()).unwrap_or(100);
    let lambda: EMProb    = qs.get("lambda").and_then(|s| s.parse().ok()).unwrap_or(20.0);
    let gamma: EMProb     = qs.get("gamma").and_then(|s| s.parse().ok()).unwrap_or(1e-20);
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

    let (matches_tsv, posteriors_tsv, props_tsv, em_likelihoods, log_str) = run_query(
        ref_path.to_str().unwrap(),
        r1_path.to_str().unwrap(),
        r2_path.to_str().unwrap(),
        mismatch, cutoff, tau, num_iter, lambda, gamma, use_penalty, threads,
    )?;

    let convergence_json = format!("[{}]",
        em_likelihoods.iter()
            .map(|v| if v.is_finite() { format!("{:.10e}", v) } else { "null".to_string() })
            .collect::<Vec<_>>()
            .join(",")
    );

    Ok(format!(
        r#"{{"ok":true,"matches":"{}","posteriors":"{}","props":"{}","convergence":{},"log":"{}"}}"#,
        json_escape(&matches_tsv),
        json_escape(&posteriors_tsv),
        json_escape(&props_tsv),
        convergence_json,
        json_escape(&log_str),
    ))
}

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
                .arg(arg!(--tau <TAU>"Minimum match log-probability threshold")
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
                .arg(arg!(-c --cutoff <CUTOFF>"Cutoff likelihood for dropping alignments")
                    .default_value("1e-4")
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
                .arg(arg!(-g --gamma <GAMMA>"penalty weight")
                    .default_value("1e-20")
                    .value_parser(clap::value_parser!(EMProb))
                    )
                .arg(arg!(-l --lambda <LAMBDA>"penalty weight") 
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
                .arg(arg!(--tau <TAU>"Minimum match log-probability threshold")
                    .default_value("1e-18")
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
            let tau = *sub_m.get_one::<EMProb>("tau").expect("required");
            let outfile = sub_m.get_one::<String>("out").unwrap().as_str();
            let threads = *sub_m.get_one::<usize>("threads").expect("required");

            let now = Instant::now();
            let (tsv, _) = run_alignment(ref_file, r1_file, r2_file, percent_mismatch, tau, threads)?;
            std::fs::write(outfile, tsv.as_bytes())?;
            println!("Alignment written to {} ({:.2?})", outfile, now.elapsed());
        },
        Some(("query",  sub_m)) => {
            let ref_file = sub_m.get_one::<String>("source").expect("required").as_str();
            let r1_file = sub_m.get_one::<String>("r1").expect("required").as_str();
            let r2_file = sub_m.get_one::<String>("r2").expect("required").as_str();
            let num_iter = *sub_m.get_one::<usize>("iter").expect("required");
            let percent_mismatch = *sub_m.get_one::<EMProb>("percent_mismatch").expect("required");
            let cutoff = *sub_m.get_one::<EMProb>("cutoff").expect("required");
            let tau = *sub_m.get_one::<EMProb>("tau").expect("required");
            let gamma = *sub_m.get_one::<EMProb>("gamma").expect("required");
            let lambda = *sub_m.get_one::<EMProb>("lambda").expect("required");
            let outfile = sub_m.get_one::<String>("out").unwrap().as_str();
            let use_penalty = *sub_m.get_one::<bool>("no-penalty").unwrap();
            let threads = *sub_m.get_one::<usize>("threads").expect("required");

            let now = Instant::now();
            let (matches_tsv, posteriors_tsv, props_tsv, _, _) = run_query(
                ref_file, r1_file, r2_file,
                percent_mismatch, cutoff, tau, num_iter, lambda, gamma, use_penalty, threads,
            )?;
            std::fs::write(format!("{}.matches", outfile), matches_tsv.as_bytes())?;
            std::fs::write(format!("{}.posteriors", outfile), posteriors_tsv.as_bytes())?;
            std::fs::write(format!("{}.props", outfile), props_tsv.as_bytes())?;
            println!("Query written to {}.{{matches,posteriors,props}} ({:.2?})", outfile, now.elapsed());
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
