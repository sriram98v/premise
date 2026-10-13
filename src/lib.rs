extern crate clap;
pub mod align;
pub mod em;
pub mod group;
pub mod io;
pub mod server;
pub mod tree;
pub mod utils;

pub use align::{merge_read_pairs, query_fastq, query_read, query_smems};
pub use em::CsrLikelihood;
pub use group::Grouping;
pub use haystackfm::BidirFmIndex as RefIndex;
pub use haystackfm::SeqId;
pub use io::{build_index_from_bytes, build_index_from_bytes_with, load_index, QueryWriters};
pub use tree::{Partition, TaxTree};
pub use utils::{EMProb, QueryProgress, ReadPair};

use anyhow::{Context, Result};
use bio::stats::LogProb;
use chrono::Local;
use clap::{arg, Arg, ArgAction, Command};
use em::{assign_reads, estimate_proportions, EmParams, Posteriors};
use haystackfm::alphabet::decode_char;
use io::{
    all_read_ids, load_fastq_forward, load_fastq_reverse, pair_reads, read_index_file,
    write_alignments, write_aligns, write_matches, write_nodes, write_posteriors, write_props,
};
use phylo::prelude::RootedTree;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::Arc;
use std::thread;
use std::time::Instant;
use tree::{contract_to_observed, reference_tree, TaxNodes};

/// `build --lookup_depth` default, as the `'static` string clap wants.
static DEFAULT_LOOKUP_DEPTH_ARG: std::sync::LazyLock<String> =
    std::sync::LazyLock::new(|| io::DEFAULT_LOOKUP_DEPTH.to_string());

/// Resolve the worker-thread count
fn resolve_thread_count(threads: usize) -> usize {
    if threads != 0 {
        return threads;
    }
    match thread::available_parallelism() {
        Ok(n) => n.get(),
        Err(e) => {
            eprintln!(
                "Warning: could not determine available parallelism ({}); falling back to 1 thread",
                e
            );
            1
        }
    }
}

/// Run pairwise alignment of reads against the reference index and write a TSV.
pub fn run_alignment<W: Write>(
    ref_file: &str,
    r1_file: &str,
    r2_file: &str,
    mem_seed_length: usize,
    eps_2: EMProb,
    threads: usize,
    out: &mut W,
) -> Result<String> {
    let num_threads = resolve_thread_count(threads);
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global();

    let forward = load_fastq_forward(r1_file)?;
    let reverse = load_fastq_reverse(r2_file)?;
    let (all_reads, unpaired) = pair_reads(forward, reverse, r1_file, r2_file)?;
    let read_ids = all_read_ids(&all_reads, &unpaired);
    let num_reads = read_ids.len();

    let fmidx = read_index_file(ref_file)?;

    let log_str = format!(
        "Timestamp: {}\nNum Threads: {}\nMin matching run: {}\nEps_2: {:e}",
        Local::now().format("%Y-%m-%d %H:%M:%S"),
        num_threads,
        mem_seed_length,
        eps_2,
    );
    println!("{}", log_str);

    let (out_alignments, _) = query_fastq(
        &fmidx,
        &all_reads,
        mem_seed_length,
        LogProb(f64::NEG_INFINITY),
        LogProb(eps_2.ln()),
        None,
    )?;

    write_alignments(out, &fmidx, &read_ids, &out_alignments)?;
    println!(
        "{} of {} reads could not be classified.",
        num_reads - out_alignments.len(),
        num_reads
    );
    Ok(log_str)
}

/// Full PREMISE query pipeline. `seq2taxid` is an optional `seqid<whitespace>taxid` map
/// grouping sequences into taxa; `nodes` an optional taxonomy (`nodes.dmp`) over those
/// taxids, which makes the EM hierarchical (see [`tree`]); `gamma` the confidence a
/// read's subtree posterior must reach for the read to descend to a vertex.
#[allow(clippy::too_many_arguments)]
pub fn run_query<W: Write>(
    ref_file: &str,
    r1_file: &str,
    r2_file: &str,
    mem_seed_length: usize,
    eps_1: EMProb,
    eps_2: EMProb,
    num_iter: usize,
    rho: EMProb,
    omega: EMProb,
    em_threshold: EMProb,
    use_penalty: bool,
    threads: usize,
    seq2taxid: Option<&str>,
    nodes: Option<&str>,
    gamma: EMProb,
    progress: Option<Arc<QueryProgress>>,
    out: &mut QueryWriters<W>,
) -> Result<(Vec<EMProb>, String)> {
    if !(gamma > 0.5 && gamma <= 1.0) {
        anyhow::bail!("gamma must be in (0.5, 1], got {}", gamma);
    }
    let num_threads = resolve_thread_count(threads);
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global();

    let forward = load_fastq_forward(r1_file)?;
    let reverse = load_fastq_reverse(r2_file)?;
    let (all_reads, unpaired) = pair_reads(forward, reverse, r1_file, r2_file)?;
    let read_ids = all_read_ids(&all_reads, &unpaired);

    let fmidx = read_index_file(ref_file)?;
    let grouping = match seq2taxid {
        Some(path) => Grouping::from_map(&fmidx, path)?,
        None => Grouping::identity(&fmidx),
    };
    // The taxonomy is built before the expensive alignment so bad input fails early.
    let tax_nodes = nodes.map(TaxNodes::from_file).transpose()?;
    let ref_tree = reference_tree(&grouping, tax_nodes.as_ref())?;
    let tree_log = if tax_nodes.is_some() {
        format!(
            "\nTaxonomy: {} vertices, {} sequences at tips ({} taxids not in the nodes file, {} sequences without taxid)",
            ref_tree.n_nodes(),
            ref_tree.n_tips(),
            ref_tree.n_unknown,
            ref_tree.n_unmapped
        )
    } else if seq2taxid.is_some() {
        format!(
            "\nTaxonomy: {} vertices, {} sequences at tips under their taxids ({} sequences without taxid)",
            ref_tree.n_nodes(),
            ref_tree.n_tips(),
            ref_tree.n_unmapped
        )
    } else {
        String::new()
    };
    let rho_eff = if use_penalty { rho } else { 0.0 };

    let log_str = format!("Timestamp: {}\nNum Threads: {}\nMin matching run: {}\nEps_1: {:e} ({:.2})\nEps_2: {:e} ({:.2})\nEM Iterations: {}\nEM Threshold: {:e}\nRho: {}\nOmega: {:e}\nGamma: {}{}{}",
        Local::now().format("%Y-%m-%d %H:%M:%S"),
        num_threads,
        mem_seed_length,
        eps_1,
        eps_1.ln(),
        eps_2,
        eps_2.ln(),
        num_iter,
        em_threshold,
        rho_eff,
        omega,
        gamma,
        grouping.log_suffix(),
        tree_log,
    );
    println!("{}", log_str);

    let (out_aligns, _all_refs) = query_fastq(
        &fmidx,
        &all_reads,
        mem_seed_length,
        LogProb(eps_1.ln()),
        LogProb(eps_2.ln()),
        progress.clone(),
    )?;

    let csr = CsrLikelihood::build(&out_aligns);
    let (mut tree, mut part) = contract_to_observed(&ref_tree, csr.refs())?;
    if !tree_log.is_empty() {
        println!(
            "Contracted taxonomy: {} vertices, {} tips with reads",
            tree.num_nodes(),
            tree.get_leaf_ids().count()
        );
    }
    let fit = estimate_proportions(
        &csr,
        &mut tree,
        &mut part,
        &EmParams {
            num_iter,
            rho: rho_eff,
            omega,
            em_threshold,
        },
        progress.clone(),
    );
    if fit.pruned + fit.collapsed > 0 {
        println!(
            "Pruned {} vertices and collapsed {} into tips; {} tips remain",
            fit.pruned,
            fit.collapsed,
            tree.get_leaf_ids().count()
        );
    }
    if fit.root_collapsed {
        eprintln!(
            "Warning: every child of the root fell below rho = {}; the sample collapsed into one tip (rho is too large for these data)",
            rho_eff
        );
    }
    let posteriors = Posteriors {
        matrix: &csr,
        pi: fit.pi_tilde,
        d: fit.d,
    };
    let assignments = assign_reads(&csr, &tree, &part, &posteriors, gamma);

    write_posteriors(
        &mut out.posteriors,
        &grouping,
        &tree,
        &part,
        &read_ids,
        &out_aligns,
        &assignments,
        &posteriors,
    )?;
    let _classified_reads = write_matches(
        &mut out.matches,
        &grouping,
        &tree,
        &part,
        &read_ids,
        &out_aligns,
        &assignments,
        &posteriors,
    )?;

    write_props(&mut out.props, &tree, &fit.node_pi)?;
    write_nodes(&mut out.nodes, &tree, &part, &fit.node_pi)?;
    write_aligns(&mut out.aligns, &fmidx, &read_ids, &out_aligns)?;

    Ok((fit.ll_trace, log_str))
}

pub fn run() -> Result<()> {
    let matches = Command::new("Maximum Likelihood Metagenomic Classification")
        .version(env!("CARGO_PKG_VERSION"))
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
                .arg(arg!(--sa_sample_rate <RATE> "Keep one suffix-array entry in every RATE (1 = full suffix array)")
                    .default_value("1")
                    .value_parser(clap::value_parser!(u32).range(1..))
                )
                .arg(arg!(--lookup_depth <K> "k of the k-mer tables that seed the SMEM search (0 = none; each index half grows by ~12 x 4^K bytes)")
                    .default_value(DEFAULT_LOOKUP_DEPTH_ARG.as_str())
                    .value_parser(clap::value_parser!(u32).range(0..=i64::from(io::MAX_LOOKUP_DEPTH)))
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
                .about("Align paired-end reads to an index and write per-reference alignment likelihoods")
                .arg(arg!(-s --source <SRC_FILE> "Source index file with reference sequences(.fmidx)")
                    .required(true)
                    .value_parser(clap::value_parser!(String))
                    )
                .arg(arg!(-m --mem <MEM_SEED_LENGTH> "Minimum run of matching bases required to report an alignment")
                    .default_value("22")
                    .value_parser(clap::value_parser!(usize))
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
                .about("Classify paired-end reads: align, then estimate reference abundances by EM")
                .arg(arg!(-s --source <SRC_FILE> "Source index file with reference sequences(.fmidx)")
                    .required(true)
                    .value_parser(clap::value_parser!(String))
                    )
                .arg(arg!(-m --mem <MEM_SEED_LENGTH> "Minimum run of matching bases required to report an alignment")
                    .default_value("22")
                    .value_parser(clap::value_parser!(usize))
                    )
                .arg(arg!(--eps_1 <EPS_1>"Cutoff likelihood for dropping alignments (0 disables the cutoff)")
                    .default_value("0")
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
                    .default_value("1e-10")
                    .value_parser(clap::value_parser!(EMProb))
                    )
                .arg(arg!(--rho <RHO>"penalty weight")
                    .default_value("150")
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
                .arg(arg!(--em_threshold <EM_THRESHOLD>"EM convergence threshold: relative change in the penalized log-likelihood")
                    .default_value("1e-6")
                    .value_parser(clap::value_parser!(EMProb))
                    )
                .arg(arg!(-t --threads <THREADS>"Number of threads (defaults to 2; 0 uses maximum number of threads)")
                    .default_value("2")
                    .value_parser(clap::value_parser!(usize))
                    )
                .arg(arg!(--seq2taxid <MAP_FILE> "seqid2taxid.map (seqid<whitespace>taxid) assigning index sequences to taxa: every sequence becomes a leaf under its taxid, an internal vertex the EM can report; sequences without a taxid hang under the root")
                    .required(false)
                    .value_parser(clap::value_parser!(String))
                    )
                .arg(arg!(--nodes <NODES_FILE> "Taxonomy over the --seq2taxid taxids: NCBI nodes.dmp, or lines of `taxid parent [rank]`. The EM then runs hierarchically, pruning and reporting vertices of the tree")
                    .required(false)
                    .requires("seq2taxid")
                    .value_parser(clap::value_parser!(String))
                    )
                .arg(arg!(--gamma <GAMMA> "Confidence level in (0.5, 1]: a read descends to the deepest vertex whose subtree posterior is at least GAMMA")
                    .default_value("0.6")
                    .value_parser(clap::value_parser!(EMProb))
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

    match matches.subcommand() {
        Some(("build", sub_m)) => {
            let src_file = sub_m
                .get_one::<String>("source")
                .expect("required")
                .as_str();
            let outfile = sub_m.get_one::<String>("out").unwrap().as_str();
            let sa_sample_rate = *sub_m.get_one::<u32>("sa_sample_rate").expect("defaulted");
            let lookup_depth = *sub_m.get_one::<u32>("lookup_depth").expect("defaulted");

            let fasta_data = std::fs::read(src_file)
                .with_context(|| format!("failed to read reference FASTA '{}'", src_file))?;
            let (idx_bytes, _) =
                build_index_from_bytes_with(&fasta_data, sa_sample_rate, lookup_depth)?;

            let out_path = match outfile {
                "" => format!("{}.fmidx", src_file),
                p => p.to_string(),
            };
            std::fs::write(&out_path, &idx_bytes)
                .with_context(|| format!("failed to write index to '{}'", out_path))?;
            println!("Index written to {}", out_path);
        }
        Some(("align", sub_m)) => {
            let ref_file = sub_m
                .get_one::<String>("source")
                .expect("required")
                .as_str();
            let r1_file = sub_m.get_one::<String>("r1").expect("required").as_str();
            let r2_file = sub_m.get_one::<String>("r2").expect("required").as_str();
            let mem_seed_length = *sub_m.get_one::<usize>("mem").expect("required");
            let eps_2 = *sub_m.get_one::<EMProb>("eps_2").expect("required");
            let outfile = sub_m.get_one::<String>("out").unwrap().as_str();
            let threads = *sub_m.get_one::<usize>("threads").expect("required");

            let now = Instant::now();
            let mut out = BufWriter::new(
                File::create(outfile)
                    .with_context(|| format!("failed to create alignment output '{}'", outfile))?,
            );
            run_alignment(
                ref_file,
                r1_file,
                r2_file,
                mem_seed_length,
                eps_2,
                threads,
                &mut out,
            )
            .with_context(|| format!("alignment run failed (output '{}')", outfile))?;
            println!("Alignment written to {} ({:.2?})", outfile, now.elapsed());
        }
        Some(("query", sub_m)) => {
            let ref_file = sub_m
                .get_one::<String>("source")
                .expect("required")
                .as_str();
            let r1_file = sub_m.get_one::<String>("r1").expect("required").as_str();
            let r2_file = sub_m.get_one::<String>("r2").expect("required").as_str();
            let num_iter = *sub_m.get_one::<usize>("iter").expect("required");
            let mem_seed_length = *sub_m.get_one::<usize>("mem").expect("required");
            let eps_1 = *sub_m.get_one::<EMProb>("eps_1").expect("required");
            let eps_2 = *sub_m.get_one::<EMProb>("eps_2").expect("required");
            let omega = *sub_m.get_one::<EMProb>("omega").expect("required");
            let rho = *sub_m.get_one::<EMProb>("rho").expect("required");
            let outfile = sub_m.get_one::<String>("out").unwrap().as_str();
            let use_penalty = *sub_m.get_one::<bool>("no-penalty").unwrap();
            let em_threshold = *sub_m.get_one::<EMProb>("em_threshold").expect("required");
            let threads = *sub_m.get_one::<usize>("threads").expect("required");
            let seq2taxid = sub_m.get_one::<String>("seq2taxid").map(String::as_str);
            let nodes = sub_m.get_one::<String>("nodes").map(String::as_str);
            let gamma = *sub_m.get_one::<EMProb>("gamma").expect("defaulted");
            if !(gamma > 0.5 && gamma <= 1.0) {
                anyhow::bail!("--gamma must be in (0.5, 1], got {}", gamma);
            }

            let now = Instant::now();
            let open = |ext: &str| -> Result<BufWriter<File>> {
                let path = format!("{}.{}", outfile, ext);
                File::create(&path)
                    .map(BufWriter::new)
                    .with_context(|| format!("failed to create '{}'", path))
            };
            let mut tables = QueryWriters {
                matches: open("matches")?,
                posteriors: open("posteriors")?,
                props: open("props")?,
                nodes: open("nodes")?,
                aligns: open("aligns")?,
            };
            run_query(
                ref_file,
                r1_file,
                r2_file,
                mem_seed_length,
                eps_1,
                eps_2,
                num_iter,
                rho,
                omega,
                em_threshold,
                use_penalty,
                threads,
                seq2taxid,
                nodes,
                gamma,
                None,
                &mut tables,
            )
            .with_context(|| {
                format!(
                    "query run failed (output '{}.{{matches,posteriors,props,nodes,aligns}}')",
                    outfile
                )
            })?;
            println!(
                "Query written to {}.{{matches,posteriors,props,nodes,aligns}} ({:.2?})",
                outfile,
                now.elapsed()
            );
        }
        Some(("fasta", sub_m)) => {
            let index_file = sub_m.get_one::<String>("index").expect("required").as_str();

            let file_bytes = std::fs::read(index_file).with_context(|| {
                format!(
                    "failed to read FM-index '{}' (build it first with the `build` subcommand)",
                    index_file
                )
            })?;
            let fmidx = load_index(&file_bytes, index_file)?;

            let mut out = String::new();
            for ref_idx in (0..fmidx.num_sequences()).map(SeqId::new) {
                let header = fmidx.seq_header(ref_idx).unwrap_or("");
                let seq: String = fmidx
                    .sequence(ref_idx)
                    .unwrap_or(&[])
                    .iter()
                    .map(|&c| decode_char(c).unwrap_or('N'))
                    .collect();
                out.push_str(&format!(">{}\n{}\n", header, seq));
            }
            print!("{}", out);
        }
        Some(("inspect", sub_m)) => {
            let index_file = sub_m.get_one::<String>("index").expect("required").as_str();

            let file_bytes = std::fs::read(index_file).with_context(|| {
                format!(
                    "failed to read FM-index '{}' (build it first with the `build` subcommand)",
                    index_file
                )
            })?;
            let fmidx = load_index(&file_bytes, index_file)?;

            let ids: Vec<SeqId> = (0..fmidx.num_sequences()).map(SeqId::new).collect();
            let total_len: usize = ids
                .iter()
                .map(|id| fmidx.sequence(*id).map(|s| s.len()).unwrap_or(0))
                .sum();
            println!("Index file: {}", index_file);
            println!("Number of references: {}", ids.len());
            println!("Total sequence length: {} bp", total_len);
            println!("References:");
            for ref_idx in ids {
                let header = fmidx.seq_header(ref_idx).unwrap_or("");
                let len = fmidx.sequence(ref_idx).map(|s| s.len()).unwrap_or(0);
                println!("  [{}] {} ({} bp)", ref_idx.index(), header, len);
            }
        }
        Some(("server", sub_m)) => {
            let port = *sub_m.get_one::<u16>("port").unwrap();
            let ip = sub_m.get_one::<String>("ip").unwrap();
            server::serve(&format!("{}:{}", ip, port))?;
        }
        _ => {
            println!("No subcommand selected. Run with --help to see available commands.");
        }
    }

    Ok(())
}
