extern crate clap;
pub mod utils;
use bio::alphabets;
use bio::alphabets::Alphabet;
use bio::io::fasta;
use bio::stats::LogProb;
use clap::{arg, Command};
use dashmap::DashMap;
use indicatif::{ProgressBar, ProgressDrawTarget};
use itertools::Itertools;
use utils::*;
use std::collections::{HashSet, VecDeque};
use std::fmt::Debug;
use std::ops::IndexMut;
use std::thread;
use std::{collections::HashMap, fs::File, io::{BufReader, Write}, sync::Mutex};
use bio::io::fastq;
use indicatif::ProgressStyle;
use rayon::prelude::*;
use std::path::Path;
use anyhow::Result;
use genedex::{FmIndexConfig, alphabet, FmIndexFlat64};
use flate2::read::GzDecoder;
use chrono::Local;
use savefile::prelude::*;
use std::time::Instant;
use genedex::text_with_rank_support::{FlatTextWithRankSupport, Block64};
use std::cmp;
use std::ops::Index;

#[derive(Debug, Clone, Hash, PartialEq, Eq, Savefile)]
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
type Alignments = HashMap<RefIdx, (usize, LogProb)>;
/// type to store the best alignments of reads to references
type MatchLikelihoods = HashMap<RefIdx, LogProb>;
/// type to store all alignments of all reads to references
type ReadAlignments = HashMap<ReadID, HashMap<RefIdx, VecDeque<(usize, LogProb)>>>;

#[derive(Debug, Clone, Copy, Hash, PartialEq, Eq, Savefile)]
pub struct MEMPos{
    pub ref_start: usize,
    pub ref_end: usize,
    pub read_start: usize, 
    pub read_end: usize,
}

impl MEMPos{
    fn is_consecutive(&self, kmer: &MEMPos)->bool{
        if kmer.read_start>=self.read_start{
            if kmer.ref_start>=self.ref_start{
                return true;
            }
        }
        false
    }
    fn is_overlapping(&self, kmer: &MEMPos)->bool{
        if kmer.read_start-self.read_start == kmer.ref_start-self.ref_start{
            if kmer.read_end-self.read_end == kmer.ref_end-self.ref_end{
                return true;
            }
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

/// Dictionary of Keys (DoK) format for sparse array to store likelihoods
#[derive(Debug, Clone, Savefile, Default)]
pub struct SparseArray<T: Copy + Send + Sync + Debug>{
    values: HashMap<(ReadIdx,RefIdx), T>,
    read_idxs_ref_map: HashMap<ReadIdx, Vec<RefIdx>>,
    ref_idxs_read_map: HashMap<RefIdx, Vec<ReadIdx>>,
}

impl<T: Copy + Send + Sync + Debug> SparseArray<T>{
    fn insert(&mut self, read_idx: ReadIdx, ref_idx: RefIdx, val: T){
        self.values.insert((read_idx, ref_idx), val);
        self.read_idxs_ref_map.entry(read_idx).or_insert(vec![]).push(ref_idx);
        self.ref_idxs_read_map.entry(ref_idx).or_insert(vec![]).push(read_idx);
    }

    fn get(&self, index: &(ReadIdx, RefIdx)) -> Option<&T> {
        self.values.get(&index)
    }

    fn contains_key(&self, index: &(ReadIdx, RefIdx)) -> bool {
        self.values.contains_key(&index)
    }

    fn contains_ref(&self, index: &RefIdx) -> bool {
        self.ref_idxs_read_map.contains_key(&index)
    }

    fn get_read_idxs(&self)->HashSet<ReadIdx>{
        self.read_idxs_ref_map.keys().cloned().collect()
    }

    fn get_ref_idxs(&self)->HashSet<RefIdx>{
        self.ref_idxs_read_map.keys().cloned().collect()
    }

    fn _get_all_read_hits(&self, read_idx: &ReadIdx)->Vec<T>{
        let ref_idxs: &Vec<RefIdx> = self.read_idxs_ref_map.get(read_idx).unwrap();
        let mut out: Vec<T> = Vec::with_capacity(ref_idxs.len());

        for ref_idx in ref_idxs{
            if let Some(x) = self.values.get(&(*read_idx, *ref_idx)) {
                out.push(*x);
            }
        }
        out
    }

    fn _get_all_read_hits_par(&self, read_idx: &ReadIdx)->Vec<T>{
        let ref_idxs: &Vec<RefIdx> = self.read_idxs_ref_map.get(read_idx).unwrap();
        let out: Mutex<Vec<T>> = Mutex::new(Vec::with_capacity(ref_idxs.len()));

        ref_idxs.par_iter().for_each(|ref_idx | {
            if let Some(x) = self.values.get(&(*read_idx, *ref_idx)) {
                out.lock().unwrap().push(*x);
            }
        });
        out.into_inner().unwrap()
    }

    fn get_all_read_hits_idx(&self, read_idx: &ReadIdx)->Vec<(RefIdx, T)>{
        let ref_idxs: &Vec<RefIdx> = self.read_idxs_ref_map.get(read_idx).unwrap();
        let mut out: Vec<(RefIdx, T)> = Vec::with_capacity(ref_idxs.len());

        for ref_idx in ref_idxs{
            if let Some(x) = self.values.get(&(*read_idx, *ref_idx)) {
                out.push((*ref_idx, *x));
            }
        }
        out
    }

    fn get_all_ref_hits(&self, ref_idx: &RefIdx)->Vec<T>{
        let read_idxs: &Vec<ReadIdx> = self.ref_idxs_read_map.get(ref_idx).unwrap();
        let mut out: Vec<T> = Vec::with_capacity(read_idxs.len());

        for read_idx in read_idxs{
            if let Some(x) = self.values.get(&(*read_idx, *ref_idx)) {
                out.push(*x);
            }
        }
        out
    }

    fn _get_all_ref_hits_par(&self, ref_idx: &RefIdx)->Vec<T>{
        let read_idxs: &Vec<ReadIdx> = self.ref_idxs_read_map.get(ref_idx).unwrap();
        let out: Mutex<Vec<T>> = Mutex::new(Vec::with_capacity(read_idxs.len()));

        read_idxs.par_iter().for_each(|read_idx | {
            if let Some(x) = self.values.get(&(*read_idx, *ref_idx)) {
                out.lock().unwrap().push(*x);
            }
        });
        out.into_inner().unwrap()
    }

    fn num_reads(&self)->usize{
        self.get_read_idxs().len()
    }

    fn _num_refs(&self)->usize{
        self.get_ref_idxs().len()
    }

    fn triples(&self)->impl Iterator<Item=(&(ReadIdx, RefIdx), &T)>{
        self.values.iter()
    }

}

impl<T: Copy + Send + Sync + Debug> Index<(ReadIdx, RefIdx)> for SparseArray<T>{
    type Output = T;

    fn index(&self, index: (ReadIdx, RefIdx)) -> &Self::Output {
        self.values.get(&index).expect("Invalid key!")
    }
}

impl<T: Copy + Send + Sync + Debug> IndexMut<(ReadIdx, RefIdx)> for SparseArray<T>{
    fn index_mut(&mut self, index: (ReadIdx, RefIdx)) -> &mut Self::Output {
        self.values.get_mut(&index).expect("Invalid key!")
    }
}


fn merge_kmer_matches(kmer_matches: &VecDeque<MEMPos>)->Vec<MEMPos>{
    let mut out_vec = Vec::with_capacity(kmer_matches.len());
    out_vec.push(kmer_matches[0]);
    for kmer_match in kmer_matches{
        for running_mem in out_vec.iter_mut(){
            if running_mem.is_consecutive(kmer_match){
                if running_mem.is_overlapping(kmer_match){
                    running_mem.merge(kmer_match);
                }
            }
        }
    }

    out_vec
}

/// Filters the matches found for different kmers and removes repeated alignments.
fn clean_kmer_matches(fmidx: &FmIndexFlat64<i64>, refs: &HashMap<RefIdx, Vec<u8>>, record: &fastq::Record, percent_mismatch: &EMProb)->(HashMap<RefIdx, HashSet<usize>>, HashMap<RefIdx, Vec<MEMPos>>){
    let mut match_positions: HashMap<RefIdx, HashSet<usize>> = HashMap::new();
    let mut mems: HashMap<RefIdx, VecDeque<MEMPos>> = HashMap::new(); // Contains all the MEMs between a read and a references denoted by 
    let read_seq = record.seq();
    let kmer_size = kmer_length(read_seq.len(), *percent_mismatch);

    fmidx.locate_many(
    record.seq().windows(kmer_size).filter(|kmer| {
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
                        .or_insert(VecDeque::new())
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
fn query_read(fmidx: &FmIndexFlat64<i64>, refs: &HashMap<RefIdx, Vec<u8>>, record: &fastq::Record, percent_mismatch: &EMProb)->Result<(Alignments, MatchLikelihoods, HashMap<RefIdx, HashSet<MEMPos>>)>{

    let read_len = record.seq().len();
    let read_seq = record.seq();
    let read_qual = record.qual();
    let max_num_mismatches: usize = (read_len as EMProb * (percent_mismatch/100 as EMProb)).floor() as usize;


    let mut best_match: HashMap<RefIdx, (usize, LogProb)> = HashMap::new();
    let mut match_likelihood: HashMap<RefIdx, LogProb> = HashMap::new();
    let mut softclipped: HashMap<RefIdx, HashSet<MEMPos>> = HashMap::new();

    let (_other_matches, mems) = clean_kmer_matches(fmidx, refs, record, percent_mismatch);
    
    mems.into_iter().for_each(|hit| {
        let ref_id = hit.0;
        let ref_seq = refs.get(&ref_id).unwrap();
        let ref_len = ref_seq.len();
        
        for mem in hit.1.iter(){
            let read_pos = mem.read_start;

            if mem.ref_start<read_pos{
                softclipped.entry(ref_id)
                    .or_insert(HashSet::new())
                    .insert(*mem);
                continue;
            }
            let ref_pos = mem.ref_start-mem.read_start;
            if ref_pos+read_len<ref_len{
                let ref_match_seg = &ref_seq[ref_pos..ref_pos+read_len];
                if num_mismatches(read_seq, ref_match_seg)<=max_num_mismatches{
                    let match_log_prob = compute_match_log_prob(read_seq, read_qual, ref_match_seg);

                    best_match.entry(ref_id)
                        .and_modify(|e| {
                            if e.1<=match_log_prob{
                                *e = (ref_pos, match_log_prob);
                            }
                        })
                        .or_insert((ref_pos, match_log_prob));

                    // update match score
                    match_likelihood.entry(ref_id)
                        .and_modify(|e| {
                            *e = *e + match_log_prob;
                        })
                        .or_insert(match_log_prob);
                }
            }
        }
    });

    Ok((best_match, match_likelihood, softclipped))
}

fn process_fastq_file(fmidx: &FmIndexFlat64<i64>, 
                    refs: &HashMap<RefIdx, Vec<u8>>, 
                    fastq_file: &Path, 
                    percent_mismatch: &EMProb)-> Result<(ReadAlignments, HashSet<RefIdx>)>
{
    let fastq_records = match get_extension_from_filename(fastq_file.to_str().expect("Invalid reads file!")){
       Some("gz") => {
            let f = File::open(fastq_file)?;
            let decoder = GzDecoder::new(f);

            fastq::Reader::from_bufread(BufReader::new(decoder)).records().filter_map(|x| x.ok()).collect_vec()
            // fastq::Reader::from_bufread(GzDecoder::new(reader)).records()
       },
       Some("fastq")|Some("fq") => {
            let f = File::open(fastq_file)?;
            let reader = BufReader::new(f);
            fastq::Reader::from_bufread(reader).records().filter_map(|x| x.ok()).collect_vec()
        },
       _ => {panic!("Invalid file type for reads!")}
    };

    let out_aligns: Mutex<HashMap<ReadID, HashMap<RefIdx, VecDeque<(usize, LogProb)>>>> = Mutex::new(HashMap::new());

    let pb = ProgressBar::with_draw_target(Some(fastq_records.len() as u64), ProgressDrawTarget::stderr());
    pb.set_style(ProgressStyle::with_template("Finding pairwise alignments: {spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent}% ({eta})").unwrap());

    let all_ref_ids: Mutex<HashSet<RefIdx>> = Mutex::new(HashSet::new());

    let softclipped_alignments: Mutex<HashMap<ReadID, HashMap<RefIdx, HashSet<MEMPos>>>> = Mutex::new(HashMap::new());
    
    fastq_records
        .par_iter()
        .map(|record| {
            let (best_hits, match_likelihoods, softclipped) = query_read(fmidx, refs, record, percent_mismatch).unwrap();

            if softclipped.len()>0 && match_likelihoods.len()==0{
                softclipped_alignments.lock().unwrap()
                    .entry(ReadID(record.id().to_string()))
                    .or_insert(HashMap::new())
                    .extend(softclipped);
            }

            (record.id().to_string(), best_hits, match_likelihoods)
        })
        .for_each(|(read_id, best_hits, match_likelihoods)| {
            match_likelihoods.iter()
                .for_each(|x| {
                    let best_align = best_hits.get(x.0).unwrap();
                    all_ref_ids.lock().unwrap().insert(*x.0);

                    out_aligns.lock().unwrap().entry(ReadID(read_id.clone()))
                        .or_default()
                        .entry(*x.0).or_default().push_back((best_align.0, *x.1));

                });
            pb.inc(1);
        });

    pb.finish_with_message("");
    
    Ok((out_aligns.into_inner().unwrap(), all_ref_ids.into_inner().unwrap()))
}

fn get_proportions_par_sparse(ll_array: &SparseArray<EMProb>, num_iter: usize, _penalty_weight: EMProb, _penalty_gamma: EMProb)->(HashMap<ReadIdx, RefIdx>, HashMap<RefIdx, EMProb>, SparseArray<EMProb>){
    let num_reads = ll_array.num_reads();

    let read_idxs = ll_array.get_read_idxs();
    // let ref_idxs = ll_array.get_ref_idxs();

    let mut props: Vec<HashMap<RefIdx, EMProb>> = vec![HashMap::new();num_iter+1];
    let mut w: SparseArray<EMProb> = SparseArray::default();

    let initial_props: DashMap<RefIdx, usize> = DashMap::new();

    read_idxs.par_iter().for_each(|read_idx| {
        let row_argmax = ll_array.get_all_read_hits_idx(&read_idx).iter().max_by(|&(_, f1), &(_, f2)| {
            EMProb::total_cmp(f1, f2)
        }).unwrap().0;
        initial_props.entry(row_argmax).and_modify(|e| { *e += 1 }).or_insert(1);
    });

    let total = initial_props.iter().map(|x| *x.value()).sum::<usize>();

    for (ref_idx, count) in initial_props.into_iter(){
        props[0].insert(ref_idx, (count as EMProb)/(total as EMProb));
    }

    let pb = ProgressBar::with_draw_target(Some(num_iter as u64), ProgressDrawTarget::stderr());
    pb.set_style(ProgressStyle::with_template("Running EM: {spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent}% ({eta})").unwrap());


    for i in 0..num_iter {

        // dbg!(i, "update w");
        ll_array.triples().for_each(|(x, val)| {
            if props[i].contains_key(&x.1){
                w.insert(x.0, x.1, val*props[i].get(&x.1).unwrap_or(&(0 as EMProb)));
            }
        });

        // dbg!(i, "compute clik");
        let clik: DashMap<ReadIdx, EMProb> = DashMap::new();
        
        read_idxs
            .iter()
            .par_bridge()
            .for_each(|read_idx| {
                let pr_rspr_s = ll_array.get_all_read_hits_idx(read_idx).into_iter().map(|(ref_idx, val)| val*props[i].get(&ref_idx).unwrap_or(&(0 as EMProb))).sum::<EMProb>();
                clik.entry(*read_idx).and_modify(|e| *e += pr_rspr_s).or_insert(pr_rspr_s);
            });


        // dbg!(i, "compute posterior");
        ll_array.triples().for_each(|(x, _) | {
            if w.contains_key(x){
                w[*x] /= *clik.get(&x.0).unwrap();
            }
        });

        // dbg!(i, "update props");
        props[i+1] = ll_array.get_ref_idxs()
            .par_iter()
            .filter(|x| w.contains_ref(x))
            .map(|ref_idx| {
                let vals = w.get_all_ref_hits(ref_idx);
                (*ref_idx, vals.iter().sum::<EMProb>()/(num_reads as EMProb))
            })
            .collect();

        pb.inc(1);

    }
    pb.finish_with_message("");

    let pb = ProgressBar::with_draw_target(Some(num_reads as u64), ProgressDrawTarget::stderr());
    pb.set_style(ProgressStyle::with_template("Finding optimal assignments: {spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent}% ({eta})").unwrap());

    let ref_idxs = w.get_read_idxs();

    let results: HashMap<ReadIdx, RefIdx> = ref_idxs.iter().par_bridge().map(|read_idx| {
        let row_argmax = ll_array.get_all_read_hits_idx(read_idx).iter().max_by(|&(_, f1), &(_, f2)| {
            EMProb::total_cmp(f1, f2)
        }).unwrap().0;
        pb.inc(1);
        (*read_idx, row_argmax)
    }).collect();

    pb.finish_with_message("");

    (results, props[num_iter].clone(), w)
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
                .arg(arg!(-r --reads <READS>"Source file with read sequences(fastq or fastq.gz)")
                    .required(true)
                    .value_parser(clap::value_parser!(String))
                    )
                .arg(arg!(-o --out <OUT_FILE>"Output file")
                    .default_value("out.aligns")
                    .value_parser(clap::value_parser!(String))
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
                .arg(arg!(-r --reads <READS>"Source file with read sequences(fastq or fastq.gz)")
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
                    .default_value("0")
                    .value_parser(clap::value_parser!(EMProb))
                    )
                .arg(arg!(-o --out <OUT_FILE>"Output file")
                    .default_value("out.matches")
                    .value_parser(clap::value_parser!(String))
                    )
                .arg(arg!(-t --threads <THREADS>"Number of threads (defaults to 2; 0 uses maximum number of threads)")
                    .default_value("2")
                    .value_parser(clap::value_parser!(usize))
                    )
        )
        .about("Maximum Likelihood Metagenomic classifier using Suffix trees")
        .get_matches();
  
    match matches.subcommand(){
        Some(("build",  sub_m)) => {
            let src_file = sub_m.get_one::<String>("source").expect("required").as_str();
            let outfile = sub_m.get_one::<String>("out").unwrap().as_str();

            let f = File::open(src_file)?;
            let reader = BufReader::new(f);

            let records = fasta::Reader::new(reader).records();

            let mut nb_reads = 0;
            let mut nb_bases = 0;

            let mut ref_ids: HashMap<RefID, RefIdx> = HashMap::new();
            let mut ref_ids_rev: HashMap<RefIdx, RefID> = HashMap::new();
            let mut refs: HashMap<RefIdx, Vec<u8>> = HashMap::new();
            let mut refs_texts = vec![];

            for (idx, result) in records.enumerate() {
                let record = result.expect("Error during reference parsing");

                ref_ids.insert(RefID(record.id().to_string()), RefIdx(idx));
                ref_ids_rev.insert(RefIdx(idx), RefID(record.id().to_string()));
                refs.insert(RefIdx(idx), record.seq().to_vec());
                refs_texts.push(record.seq().to_vec());

                nb_reads += 1;
                nb_bases += record.seq().len();

            }

            let log_str = format!("Timestamp: {}\nReference Index: {}\nNum references: {}\n Num bases: {}\nOutput File: {}",
                Local::now().format("%Y-%m-%d %H:%M:%S"),
                src_file,
                nb_reads,
                nb_bases,
                outfile,
            );

            println!("{}", log_str);

            let dna_alphabet = alphabet::ascii_dna_iupac_as_dna_with_n();
            let fmidx: FmIndexFlat64<_> = FmIndexConfig::<i64, FlatTextWithRankSupport<i64, Block64>>::new()
                .suffix_array_sampling_rate(1)
                .lookup_table_depth(13)
                .construct_index(refs_texts.iter().map(|x| str::from_utf8(x).unwrap()).collect_vec(), dna_alphabet);

            let io_struct = IOFMIndex{
                fmidx,
                idx_to_id: ref_ids_rev,
                id_to_idx: ref_ids,
                idx_to_seq: refs,
            };

            
            match outfile{
                "" => {
                    save_file(format!("{}.fmidx", src_file), 0, &io_struct)?;
                },
                _ => {
                    save_file(outfile, 0, &io_struct)?;
                }
            };

        },
        Some(("align",  sub_m)) => {
            let ref_file = sub_m.get_one::<String>("source").expect("required").as_str();
            let reads_file = sub_m.get_one::<String>("reads").expect("required").as_str();
            let percent_mismatch = sub_m.get_one::<EMProb>("percent_mismatch").expect("required");
            let outfile = sub_m.get_one::<String>("out").unwrap().as_str();

            let outpath = Some(Mutex::new(File::create(outfile).unwrap()));
            let outstr = "ReadID\tRefID\tStart_pos\tLikelihood\n".to_string();
            match outpath.as_ref(){
                Some(file) => {file.lock().unwrap().write_all(outstr.as_bytes()).unwrap();},
                None => {print!("{}", outstr)},
            };


            let num_threads = match sub_m.get_one::<usize>("threads").expect("required"){
                0 => thread::available_parallelism()?.get(),
                _ => *sub_m.get_one::<usize>("threads").expect("required"),
            };
            
            rayon::ThreadPoolBuilder::new().num_threads(num_threads).build_global()?;

            let fastq_records = match get_extension_from_filename(reads_file){
                Some("gz") => {
                    let f = File::open(reads_file)?;
                    let decoder = GzDecoder::new(f);
                    fastq::Reader::from_bufread(BufReader::new(decoder)).records().filter_map(|x| x.ok()).collect_vec()
                },
                Some("fastq")|Some("fq") => {
                    let f = File::open(reads_file)?;
                    let reader = BufReader::new(f);
                    fastq::Reader::from_bufread(reader).records().filter_map(|x| x.ok()).collect_vec()
                },
                _ => panic!("Invalid file type for reads!")
            };

            let read_len = fastq_records.iter().map(|x| x.seq().len()).sum::<usize>() as f32/(fastq_records.len() as f32);
            let all_read_ids: HashSet<ReadID> = fastq_records.iter().map(|x| ReadID(x.id().to_string())).collect();

            let log_str = format!("Timestamp: {}\nReads File: {}\nReference Index: {}\nNum Threads: {}\nPercent Mismatch: {}\nOutput File: {}",
                Local::now().format("%Y-%m-%d %H:%M:%S"),
                reads_file,
                ref_file,
                num_threads,
                percent_mismatch,
                outfile,
            );

            let num_reads = all_read_ids.len();

            println!("{}", log_str);
            println!("Num reads: {}", num_reads);
            println!("Average read len: {:.3} bp", read_len);

            let iofmidx: IOFMIndex = load_file(ref_file, 0)?;

            let fmidx: FmIndexFlat64<i64> = iofmidx.fmidx;

            let _ref_ids = iofmidx.id_to_idx;
            let ref_ids_rev = iofmidx.idx_to_id;
            let refs = iofmidx.idx_to_seq;

            let now = Instant::now();

            let (out_alignments, _all_refs) = process_fastq_file(&fmidx, &refs, Path::new(reads_file), percent_mismatch)?;

            let mut read_ids: HashMap<ReadID, ReadIdx> = HashMap::new();
            let mut read_ids_rev: HashMap<ReadIdx, ReadID> = HashMap::new();
            for (n, read_id) in out_alignments.keys().enumerate() {
                read_ids.insert(read_id.clone(), ReadIdx(n));
                read_ids_rev.insert(ReadIdx(n), read_id.clone());
            }

            println!("Number of reads aligned: {} ({:.2}%)", out_alignments.len(), (out_alignments.len() as f64/num_reads as f64)*100_f64);

            let pb = ProgressBar::new(out_alignments.len() as u64);
            pb.set_style(ProgressStyle::with_template("Writing output: {spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent}% ({eta})").unwrap());


            for read_id in all_read_ids{
                if !out_alignments.contains_key(&read_id){
                    let outstr = format!("{}\t{}\t{}\t{:.5}\n", 
                        *read_id,
                        "unclassified".to_string(), 
                        "-".to_string(), 
                        "-".to_string()
                    );
                    match outpath.as_ref(){
                        Some(file) => {file.lock().unwrap().write_all(outstr.as_bytes()).unwrap();},
                        None => {print!("{}", outstr)},
                    };
                    continue;
                }

                let read_idx = *read_ids.get(&read_id).unwrap();

                let read_aligns = out_alignments.get(&read_id).unwrap();
                for (ref_idx,aligns) in read_aligns{

                    let ref_id = ref_ids_rev.get(ref_idx).cloned().unwrap();

                    for (pos,likelihood) in aligns{
                        let outstr = format!("{}\t{}\t{}\t{:.5}\n", 
                            *read_ids_rev.get(&read_idx).unwrap().clone(),
                            *ref_id.clone(), 
                            pos, 
                            likelihood.exp(),
                        );
                        match outpath.as_ref(){
                            Some(file) => {file.lock().unwrap().write_all(outstr.as_bytes()).unwrap();},
                            None => {print!("{}", outstr)},
                        };
                    }

                }

            }

            println!("{} reads ({:.3}%) could not be classified!", num_reads-out_alignments.len(), (((num_reads-out_alignments.len()) as f64)/num_reads as f64)*100_f64);
            let elapsed_time = now.elapsed();
            println!("Total runtime: {:.2?}", elapsed_time);

            pb.finish_with_message("");
        },
        Some(("query",  sub_m)) => {
            let ref_file = sub_m.get_one::<String>("source").expect("required").as_str();
            let reads_file = sub_m.get_one::<String>("reads").expect("required").as_str();
            let num_iter = sub_m.get_one::<usize>("iter").expect("required");
            let percent_mismatch = sub_m.get_one::<EMProb>("percent_mismatch").expect("required");
            let cutoff = *sub_m.get_one::<EMProb>("cutoff").expect("required");
            let gamma = sub_m.get_one::<EMProb>("gamma").expect("required");
            let lambda = sub_m.get_one::<EMProb>("lambda").expect("required");
            let outfile = sub_m.get_one::<String>("out").unwrap().as_str();

            let outpath = Some(Mutex::new(File::create(outfile).unwrap()));
            let outstr = "ReadID\tRefID\tStart_pos\tPosterior\n".to_string();
            match outpath.as_ref(){
                Some(file) => {file.lock().unwrap().write_all(outstr.as_bytes()).unwrap();},
                None => {print!("{}", outstr)},
            };


            let num_threads = match sub_m.get_one::<usize>("threads").expect("required"){
                0 => thread::available_parallelism()?.get(),
                _ => *sub_m.get_one::<usize>("threads").expect("required"),
            };
            
            rayon::ThreadPoolBuilder::new().num_threads(num_threads).build_global()?;

            let fastq_records = match get_extension_from_filename(reads_file){
                Some("gz") => {
                    let f = File::open(reads_file)?;
                    let decoder = GzDecoder::new(f);
                    fastq::Reader::from_bufread(BufReader::new(decoder)).records().filter_map(|x| x.ok()).collect_vec()
                },
                Some("fastq")|Some("fq") => {
                    let f = File::open(reads_file)?;
                    let reader = BufReader::new(f);
                    fastq::Reader::from_bufread(reader).records().filter_map(|x| x.ok()).collect_vec()
                },
                _ => panic!("Invalid file type for reads!")
            };

            let read_len = fastq_records.iter().map(|x| x.seq().len()).sum::<usize>() as f32/(fastq_records.len() as f32);
            let all_read_ids: HashSet<ReadID> = fastq_records.iter().map(|x| ReadID(x.id().to_string())).collect();

            let log_str = format!("Timestamp: {}\nReads File: {}\nReference Index: {}\nNum Threads: {}\nPercent Mismatch: {}\nCutoff likelihood: {:e}\nEM Iterations: {}\nOutput File: {}",
                Local::now().format("%Y-%m-%d %H:%M:%S"),
                reads_file,
                ref_file,
                num_threads,
                percent_mismatch,
                cutoff,
                num_iter,
                outfile,
            );

            let num_reads = all_read_ids.len();

            println!("{}", log_str);
            println!("Num reads: {}", num_reads);
            println!("Average read len: {:.3} bp", read_len);

            let iofmidx: IOFMIndex = load_file(ref_file, 0)?;

            let fmidx: FmIndexFlat64<i64> = iofmidx.fmidx;

            let _ref_ids = iofmidx.id_to_idx;
            let ref_ids_rev = iofmidx.idx_to_id;
            let refs = iofmidx.idx_to_seq;

            let now = Instant::now();

            let (out_aligns, _all_refs) = process_fastq_file(&fmidx, &refs, Path::new(reads_file), percent_mismatch)?;

            let pb = ProgressBar::with_draw_target(Some(out_aligns.len() as u64), ProgressDrawTarget::stderr());
            pb.set_style(ProgressStyle::with_template("Filtering alignments: {spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent}% ({eta})").unwrap());

            let all_refs_filtered: Mutex<HashSet<RefIdx>> = Mutex::new(HashSet::new());
            // let out_alignments: ReadAlignments = 
            out_aligns.iter().for_each(|(_, aligns)| {
                let mut alignment_dist = aligns.iter()
                    .map(|(ref_idx, hits)| (*ref_idx, hits.iter().map(|(pos, ll)| (*pos, ll.0 as EMProb)).max_by(|a, b| a.1.total_cmp(&b.1)).unwrap()))
                    .collect_vec();
                alignment_dist.sort_by(|a, b| a.1.1.total_cmp(&b.1.1));
                alignment_dist.reverse();

                let num_alignments = alignment_dist.len();

                let read_alignment_max_likelihood = alignment_dist[0].1.1.exp();
                
                let mut n = 0;
                while n<num_alignments && alignment_dist[n].1.1.exp()/read_alignment_max_likelihood>= cutoff{
                    n+=1;
                }
                let usable_alignments = n;

                let best_aligns: HashMap<RefIdx, VecDeque<(usize, LogProb)>> = alignment_dist.iter()
                    .take(usable_alignments)
                    .map(|x| (x.0, VecDeque::from([(x.1.0, LogProb(x.1.1.into()))])))
                    .collect();
                for ref_idx in best_aligns.keys(){
                    all_refs_filtered.lock().unwrap().insert(*ref_idx);
                }
                pb.inc(1);
            });

            pb.finish_with_message("");

            let out_alignments: ReadAlignments = out_aligns.into_iter()
                .map(|(read_id, mut aligns)| {
                    aligns.retain(|&k, _| all_refs_filtered.lock().unwrap().contains(&k));
                    (read_id, aligns)
                })
                .collect();


            let mut em_ref_ids: HashMap<usize, RefIdx> = HashMap::new();
            // let mut em_ref_ids_rev: HashMap<RefIdx, usize> = HashMap::new();
            for (n,ref_idx) in all_refs_filtered.into_inner().unwrap().into_iter().enumerate(){
                em_ref_ids.insert(n, ref_idx);
                // em_ref_ids_rev.insert(ref_idx, n);
            }
            let num_em_refs = em_ref_ids.len();


            let mut read_ids: HashMap<ReadID, ReadIdx> = HashMap::new();
            let mut read_ids_rev: HashMap<ReadIdx, ReadID> = HashMap::new();
            for (n, read_id) in out_alignments.keys().enumerate() {
                read_ids.insert(read_id.clone(), ReadIdx(n));
                read_ids_rev.insert(ReadIdx(n), read_id.clone());
            }

            println!("Number of reads aligned: {} ({:.2}%)", out_alignments.len(), (out_alignments.len() as f64/num_reads as f64)*100_f64);
            println!("Number of potential sources: {}", num_em_refs);

            let mut ll_array: SparseArray<EMProb> = SparseArray::default();

            for (read_id, alignments) in out_alignments.iter(){
                for (ref_idx, positions) in alignments{
                    let read_idx = read_ids.get(read_id).unwrap();

                    for (_, score) in positions.iter(){
                        ll_array.insert(*read_idx, *ref_idx, score.exp() as EMProb);
                    }
                }
            }


            let (read_assignments, _props, posteriors) = get_proportions_par_sparse(&ll_array, *num_iter, *lambda, *gamma);

            let pb = ProgressBar::new(read_assignments.len() as u64);
            pb.set_style(ProgressStyle::with_template("Writing output: {spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent}% ({eta})").unwrap());


            for read_id in all_read_ids{
                if !out_alignments.contains_key(&read_id){
                    let outstr = format!("{}\t{}\t{}\t{:.5}\n", 
                        *read_id,
                        "unclassified".to_string(), 
                        "-".to_string(), 
                        "-".to_string()
                    );
                    match outpath.as_ref(){
                        Some(file) => {file.lock().unwrap().write_all(outstr.as_bytes()).unwrap();},
                        None => {print!("{}", outstr)},
                    };
                    continue;
                }

                let read_idx = *read_ids.get(&read_id).unwrap();
                let ref_idx = read_assignments.get(&read_idx).unwrap();
                let ref_id = ref_ids_rev.get(ref_idx).cloned().unwrap();
                let read_id = read_ids_rev.get(&read_idx).unwrap();

                if !out_alignments.get(read_id).unwrap().contains_key(ref_idx){
                    continue;
                }

                if out_alignments.get(read_id).unwrap().get(ref_idx).unwrap().front().is_some() && ll_array[(read_idx, *ref_idx)].is_finite(){
                    let outstr = format!("{}\t{}\t{}\t{:.5e}\n", 
                        *read_ids_rev.get(&read_idx).unwrap().clone(),
                        *ref_id.clone(), 
                        out_alignments.get(read_id).unwrap().get(ref_idx).unwrap().front().unwrap().0, 
                        posteriors.get(&(read_idx, *ref_idx)).expect("Read-reference alignment was never found!"),
                    );
                    match outpath.as_ref(){
                        Some(file) => {file.lock().unwrap().write_all(outstr.as_bytes()).unwrap();},
                        None => {print!("{}", outstr)},
                    };
                }
            }

            println!("{} reads ({:.3}%) could not be classified!", num_reads-out_alignments.len(), (((num_reads-out_alignments.len()) as f64)/num_reads as f64)*100_f64);
            let elapsed_time = now.elapsed();
            println!("Total runtime: {:.2?}", elapsed_time);

            pb.finish_with_message("");
        },
        _ => {
            println!("No option selected! Refer help page (-h flag)");
        }
    }
    
    Ok(())
}
