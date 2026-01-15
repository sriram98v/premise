extern crate blas_src;
extern crate clap;
extern crate shrust;
pub mod utils;
use bio::alphabets;
use bio::alphabets::Alphabet;
use bio::io::fasta;
use bio::stats::Prob;
use clap::{arg, Command};
use indicatif::{ProgressBar, ProgressDrawTarget};
use itertools::Itertools;
use ndarray::Array2;
use ndarray::prelude::*;
use utils::*;
use std::collections::{HashSet, VecDeque};
use std::thread;
use std::{collections::HashMap, fs::File, io::{BufReader, Write}, sync::Mutex};
use bio::io::fastq;
use indicatif::ProgressStyle;
use rayon::prelude::*;
use std::path::Path;
use anyhow::Result;
use genedex::{FmIndexConfig, alphabet, FmIndexFlat64};
use flate2::read::GzDecoder;
use ndarray_stats::QuantileExt;
use chrono::Local;
use savefile::prelude::*;
use std::time::Instant;
use genedex::text_with_rank_support::{FlatTextWithRankSupport, Block64};

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

/// type to store all alignments of a read to references
type Alignments = HashMap<RefIdx, (usize, Prob)>;
/// type to store the best alignments of reads to references
type MatchLikelihoods = HashMap<RefIdx, Prob>;
/// type to store all alignments of all reads to references
type ReadAlignments = HashMap<ReadID, HashMap<RefIdx, VecDeque<(usize, Prob)>>>;

/// Filters the matches found for different kmers and removes repeated alignments.
fn clean_kmer_matches(fmidx: &FmIndexFlat64<i64>, refs: &HashMap<RefIdx, Vec<u8>>, record: &fastq::Record, percent_mismatch: &f32)->HashMap<RefIdx, HashSet<usize>>{
    let mut match_positions: HashMap<RefIdx, HashSet<usize>> = HashMap::new();
    let read_seq = record.seq();
    let kmer_size = kmer_length(read_seq.len(), *percent_mismatch);



    fmidx.locate_many(
    record.seq().windows(kmer_size).filter(|kmer| {
                let read_alphabet = Alphabet::new(*kmer);
                let dna_alphabet = alphabets::dna::alphabet();
                // dbg!(&read_alphabet.symbols, &dna_alphabet.symbols);
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

                    if ref_start_pos_kmer>=kmer_start_read{
                        // start position of alignment in reference for read
                        let align_start = ref_start_pos_kmer-kmer_start_read;
                        match_positions.entry(seq_idx).and_modify(|align_pos| {align_pos.insert(align_start);}).or_default().insert(align_start);
                    }


                })
        });

    match_positions
}

/// Aligns a single read to each of the references
/// Returns a pair of Hashmaps. The first maps the read to its best alignment to each reference (reference_name, (alignment_start_pos, likelihood of alignment)).
/// The second returns the sum of likelihoods of all alignments to each reference.(reference_name, (sum of likelihood of all alignments)). 
fn query_read(fmidx: &FmIndexFlat64<i64>, refs: &HashMap<RefIdx, Vec<u8>>, record: &fastq::Record, percent_mismatch: &f32)->Result<(Alignments, MatchLikelihoods)>{

    let read_len = record.seq().len();
    let read_seq = record.seq();
    let read_qual = record.qual();
    let max_num_mismatches: usize = (read_len as f32 * (percent_mismatch/100_f32)).floor() as usize;

    let best_match: Mutex<HashMap<RefIdx, (usize, Prob)>> = Mutex::new(HashMap::new());
    let match_likelihood: Mutex<HashMap<RefIdx, Prob>> = Mutex::new(HashMap::new());


    // let matches = match_read_kmers(fmidx, record, percent_mismatch)?;

    let other_matches = clean_kmer_matches(fmidx, refs, record, percent_mismatch);
    
    other_matches.into_iter().for_each(|hit| {
        let ref_id = hit.0;
        let ref_seq = refs.get(&ref_id).unwrap();
        let ref_len = ref_seq.len();
        
        for ref_pos in hit.1.iter(){
            if ref_pos+read_len<ref_len{
                let ref_match_seg = &ref_seq[*ref_pos..ref_pos+read_len];
                if num_mismatches(read_seq, ref_match_seg)<=max_num_mismatches{
                    // dbg!(str::from_utf8(read_seq).unwrap(), str::from_utf8(ref_match_seg).unwrap());
                    let match_log_prob = compute_match_log_prob(read_seq, read_qual, ref_match_seg);

                    best_match.lock().unwrap().entry(ref_id)
                        .and_modify(|e| {
                            if e.1<=match_log_prob{
                                *e = (*ref_pos, match_log_prob);
                            }
                        })
                        .or_insert((*ref_pos, match_log_prob));

                    // update match score
                    match_likelihood.lock().unwrap().entry(ref_id)
                        .and_modify(|e| {
                            *e = *e + match_log_prob;
                        })
                        .or_insert(match_log_prob);
                }
            }
        }
    });

    Ok((best_match.into_inner().unwrap(), match_likelihood.into_inner().unwrap()))
}

fn process_fastq_file(fmidx: &FmIndexFlat64<i64>, 
                    refs: &HashMap<RefIdx, Vec<u8>>, 
                    fastq_file: &Path, 
                    percent_mismatch: &f32)-> Result<(ReadAlignments, HashSet<RefIdx>)>
{
    let fastq_records = match get_extension_from_filename(fastq_file.to_str().expect("Invalid reads file!")){
       Some("gz") => {
            let f = File::open(fastq_file)?;
            let decoder = GzDecoder::new(f);

            fastq::Reader::from_bufread(BufReader::new(decoder)).records().map(|x| x.unwrap()).collect_vec()
            // fastq::Reader::from_bufread(GzDecoder::new(reader)).records()
       },
       Some("fastq")|Some("fq") => {
            let f = File::open(fastq_file)?;
            let reader = BufReader::new(f);
            fastq::Reader::from_bufread(reader).records().map(|x| x.unwrap()).collect_vec()
        },
       _ => {panic!("Invalid file type for reads!")}
    };

    let out_aligns: Mutex<HashMap<ReadID, HashMap<RefIdx, VecDeque<(usize, Prob)>>>> = Mutex::new(HashMap::new());

    let pb = ProgressBar::with_draw_target(Some(fastq_records.len() as u64), ProgressDrawTarget::stdout());
    pb.set_style(ProgressStyle::with_template("Finding pairwise alignments: {spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent}% ({eta})").unwrap());

    let all_ref_ids: Mutex<HashSet<RefIdx>> = Mutex::new(HashSet::new());
    
    fastq_records
        .par_iter()
        .map(|record| {
            let (best_hits, match_likelihoods) = query_read(fmidx, refs, record, percent_mismatch).unwrap();
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

    // dbg!(&out_aligns);
    
    Ok((out_aligns.into_inner().unwrap(), all_ref_ids.into_inner().unwrap()))
}

fn proportions_penalty(props: &Array1<f64>, gamma: f64)->f64{
    let new_props = 1_f64 + props/gamma;
    new_props.sum().ln()
}

#[allow(unused_assignments)]
fn get_proportions(ll_array: &Array2<f64>, num_iter: usize, penalty_weight: f64, penalty_gamma: f64)->(HashMap<ReadIdx, usize>, Vec<f64>){
    let num_reads = ll_array.shape()[0];
    let num_srcs = ll_array.shape()[1];

    let mut props = Array::<f64, _>::zeros((num_iter+1, num_srcs).f());
    let mut w = Array::<f64, _>::zeros((num_reads, num_srcs).f());

    let mut initial_props: HashMap<usize, usize> = HashMap::new();

    ll_array.axis_iter(Axis(0)).for_each(|row| {
        let row_argmax = row.argmax_skipnan().unwrap();
        initial_props.entry(row_argmax).and_modify(|e| { *e += 1 }).or_insert(1);
    });

    let total = initial_props.values().sum::<usize>();

    for (ref_idx, count) in initial_props.into_iter(){
        props[[0, ref_idx]] = (count as f64)/(total as f64);
    }

    let pb = ProgressBar::with_draw_target(Some(num_iter as u64), ProgressDrawTarget::stdout());
    pb.set_style(ProgressStyle::with_template("Running EM: {spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent}% ({eta})").unwrap());

    for i in 0..num_iter {
        w = ll_array*&props.slice_mut(s![i, ..]);

        let penalty = penalty_weight*proportions_penalty(&props.slice(s![i, ..num_srcs]).into_owned(), penalty_gamma);

        w -= penalty;

        let clik = w.sum_axis(Axis(1)).insert_axis(Axis(1));

        w = w/clik;
        w.mapv_inplace(|x| if x.is_nan() { 0. } else { x });

        props.slice_mut(s![i+1, ..num_srcs]).assign(&w.mean_axis(Axis(0)).unwrap());

        pb.inc(1);

    }
    pb.finish_with_message("");

    w = ll_array*&props.slice(s![num_iter, ..]);

    let clik = w.sum_axis(Axis(1)).insert_axis(Axis(1));

    w = w/clik;
    w.mapv_inplace(|x| if x.is_nan() { 0. } else { x });

    let pb = ProgressBar::with_draw_target(Some(w.shape()[0] as u64), ProgressDrawTarget::stdout());
    pb.set_style(ProgressStyle::with_template("Finding optimal assignments: {spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent}% ({eta})").unwrap());

    let results: HashMap<ReadIdx, usize> = w.axis_iter(Axis(0)).enumerate().par_bridge().map(|(row_idx, row)| {
        // dbg!(&row);
        let row_argmax = row.argmax_skipnan().unwrap();
        pb.inc(1);
        (ReadIdx(row_idx), row_argmax)
    }).collect();

    pb.finish_with_message("");

    (results, props.slice(s![num_iter, ..]).to_vec())
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
            Command::new("query")
                .arg(arg!(-s --source <SRC_FILE> "Source index file with reference sequences(.fmidx)")
                    .required(true)
                    .value_parser(clap::value_parser!(String))
                    )
                .arg(arg!(-p --percent_mismatch <PERCENT_MISMATCH>"Percent mismatch to reference sequences")
                    .required(true)
                    .value_parser(clap::value_parser!(f32))
                    )
                .arg(arg!(-c --cutoff <CUTOFF>"Cutoff likelihood for dropping alignments")
                    .default_value("1e-4")
                    .value_parser(clap::value_parser!(f64))
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
                    .value_parser(clap::value_parser!(f64))
                    )
                .arg(arg!(-l --lambda <LAMBDA>"penalty weight") 
                    .default_value("0")
                    .value_parser(clap::value_parser!(f64))
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
        Some(("query",  sub_m)) => {
            let ref_file = sub_m.get_one::<String>("source").expect("required").as_str();
            let reads_file = sub_m.get_one::<String>("reads").expect("required").as_str();
            let num_iter = sub_m.get_one::<usize>("iter").expect("required");
            let percent_mismatch = sub_m.get_one::<f32>("percent_mismatch").expect("required");
            let cutoff = *sub_m.get_one::<f64>("cutoff").expect("required");
            let gamma = sub_m.get_one::<f64>("gamma").expect("required");
            let lambda = sub_m.get_one::<f64>("lambda").expect("required");
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
                    fastq::Reader::from_bufread(BufReader::new(decoder)).records().map(|x| x.unwrap()).collect_vec()
                },
                Some("fastq")|Some("fq") => {
                    let f = File::open(reads_file)?;
                    let reader = BufReader::new(f);
                    fastq::Reader::from_bufread(reader).records().map(|x| x.unwrap()).collect_vec()
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

            println!("{}", log_str);
            println!("Num reads: {}", fastq_records.len());
            println!("Average read len: {:.3} bp", read_len);

            let iofmidx: IOFMIndex = load_file(ref_file, 0)?;

            let fmidx: FmIndexFlat64<i64> = iofmidx.fmidx;

            let _ref_ids = iofmidx.id_to_idx;
            let ref_ids_rev = iofmidx.idx_to_id;
            let refs = iofmidx.idx_to_seq;

            let now = Instant::now();

            let (out_aligns, _all_refs) = process_fastq_file(&fmidx, &refs, Path::new(reads_file), percent_mismatch)?;

            let pb = ProgressBar::with_draw_target(Some(out_aligns.len() as u64), ProgressDrawTarget::stdout());
            pb.set_style(ProgressStyle::with_template("Filtering alignments: {spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent}% ({eta})").unwrap());

            let all_refs_filtered: Mutex<HashSet<RefIdx>> = Mutex::new(HashSet::new());
            let out_alignments: ReadAlignments = out_aligns.into_iter().map(|(read_id, aligns)| {
                let mut alignment_dist = aligns.iter()
                    .map(|(ref_idx, hits)| (*ref_idx, hits.iter().map(|(pos, ll)| (*pos, ll.0)).max_by(|a, b| a.1.total_cmp(&b.1)).unwrap()))
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

                let best_aligns: HashMap<RefIdx, VecDeque<(usize, Prob)>> = alignment_dist.iter()
                    .take(usable_alignments)
                    .map(|x| (x.0, VecDeque::from([(x.1.0, Prob(x.1.1))])))
                    .collect();
                for ref_idx in best_aligns.keys(){
                    all_refs_filtered.lock().unwrap().insert(*ref_idx);
                }
                pb.inc(1);
                (read_id, best_aligns)
            }).collect();

            pb.finish_with_message("");



            let mut em_ref_ids: HashMap<usize, RefIdx> = HashMap::new();
            let mut em_ref_ids_rev: HashMap<RefIdx, usize> = HashMap::new();
            for (n,ref_idx) in all_refs_filtered.into_inner().unwrap().into_iter().enumerate(){
                em_ref_ids.insert(n, ref_idx);
                em_ref_ids_rev.insert(ref_idx, n);
            }


            let mut read_ids: HashMap<ReadID, ReadIdx> = HashMap::new();
            let mut read_ids_rev: HashMap<ReadIdx, ReadID> = HashMap::new();
            for (n, read_id) in out_alignments.keys().enumerate() {
                read_ids.insert(read_id.clone(), ReadIdx(n));
                read_ids_rev.insert(ReadIdx(n), read_id.clone());
            }

            println!("Number of reads aligned: {} ({:.2}%)", out_alignments.len(), (out_alignments.len() as f64/fastq_records.len() as f64)*100_f64);
            println!("Number of potential sources: {}", em_ref_ids.len());
            println!("Min Mem required: {:.3} Gb", em_ref_ids.len() as f64*read_ids.len() as f64*64_f64*1e-9);

            let mut ll_array: ArrayBase<ndarray::OwnedRepr<f64>, Dim<[usize; 2]>> = Array2::<f64>::zeros((read_ids.len(), em_ref_ids.len()));
            ll_array.fill(f64::NEG_INFINITY);

            for (read_id, alignments) in out_alignments.iter(){
                for (ref_idx, positions) in alignments{
                    let em_ref_idx = em_ref_ids_rev.get(ref_idx).unwrap();
                    let read_idx = read_ids.get(read_id).unwrap();

                    for (_, score) in positions.iter(){
                        ll_array[[**read_idx, *em_ref_idx]] = **score;
                    }
                }
            }


            let (read_assignments, _props) = get_proportions(&ll_array.exp(), *num_iter, *lambda, *gamma);

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
                let em_ref_idx = read_assignments.get(&read_idx).unwrap();
                let ref_idx = em_ref_ids.get(&em_ref_idx).unwrap();
                let ref_id = ref_ids_rev.get(ref_idx).cloned().unwrap();
                let read_id = read_ids_rev.get(&read_idx).unwrap();

                if !out_alignments.get(read_id).unwrap().contains_key(ref_idx){
                    continue;
                }

                if out_alignments.get(read_id).unwrap().get(ref_idx).unwrap().front().is_some() && ll_array[[*read_idx, *em_ref_idx]].is_finite(){
                    let outstr = format!("{}\t{}\t{}\t{:.5}\n", 
                        *read_ids_rev.get(&read_idx).unwrap().clone(),
                        *ref_id.clone(), 
                        out_alignments.get(read_id).unwrap().get(ref_idx).unwrap().front().unwrap().0, 
                        ll_array[[*read_idx, *em_ref_idx]]
                    );
                    match outpath.as_ref(){
                        Some(file) => {file.lock().unwrap().write_all(outstr.as_bytes()).unwrap();},
                        None => {print!("{}", outstr)},
                    };
                }
            }

            println!("{} reads ({:.3}%) could not be classified!", fastq_records.len()-out_alignments.len(), (((fastq_records.len()-out_alignments.len()) as f64)/fastq_records.len() as f64)*100_f64);
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
