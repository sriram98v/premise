extern crate clap;
pub mod align;
pub mod analyze;
pub mod em;
pub mod group;
pub mod io;
pub mod refit;
pub mod report;
pub mod server;
pub mod tree;
pub mod utils;

pub use align::{merge_read_pairs, query_fastq, query_read, query_smems};
pub use em::CsrLikelihood;
pub use group::Grouping;
pub use haystackfm::BidirRlFmIndex as RefIndex;
pub use haystackfm::SeqId;
pub use io::{build_index_from_bytes, build_index_from_bytes_with, load_index};
pub use tree::{Partition, TaxTree};
pub use utils::{EMProb, QueryProgress, ReadPair};

use analyze::{run_analyze, AnalyzeOutputs, AnalyzeParams};
use anyhow::{Context, Result};
use bio::stats::LogProb;
use chrono::Local;
use clap::{arg, Arg, ArgAction, Command};
use em::Init;
use haystackfm::alphabet::decode_char;
use io::{
    all_read_ids, load_fastq_forward, load_fastq_reverse, pair_reads, read_index_file,
    write_alignments,
};
use refit::{analyze_report, RefitParams, TableRequest};
use report::nexus::Sink;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::thread;
use std::time::Instant;

/// `build --lookup_depth` default, as the `'static` string clap wants.
static DEFAULT_LOOKUP_DEPTH_ARG: std::sync::LazyLock<String> =
    std::sync::LazyLock::new(|| io::DEFAULT_LOOKUP_DEPTH.to_string());

/// Resolve the worker-thread count
pub(crate) fn resolve_thread_count(threads: usize) -> usize {
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

/// The EM arguments, without defaults: `query` adds its own, `analyze` takes the
/// report's for anything left unset.
fn em_args() -> [Arg; 5] {
    [
        arg!(-i --iter <ITER>"Number of iterations for EM")
            .value_parser(clap::value_parser!(usize)),
        arg!(--omega <OMEGA>"Penalty weight omega (used when --rho > 0)")
            .value_parser(clap::value_parser!(EMProb)),
        arg!(--rho <RHO>"Penalty weight rho: > 0 turns the penalty on, 0 leaves it off")
            .value_parser(clap::value_parser!(EMProb)),
        arg!(--em_threshold <EM_THRESHOLD>"EM convergence threshold: relative change in the penalized log-likelihood")
            .value_parser(clap::value_parser!(EMProb)),
        arg!(--init <INIT> "Initial conditionals: `split` (each read's best-scoring sequences share one count, ties split evenly) or `uniform` (1/|children| at every vertex)")
            .value_parser(["split", "uniform"]),
    ]
}

fn threads_arg() -> Arg {
    arg!(-t --threads <THREADS>"Number of threads (defaults to 2; 0 uses maximum number of threads)")
        .default_value("2")
        .value_parser(clap::value_parser!(usize))
}

/// The `query` arguments (everything but `--gamma` and `--matches`).
fn query_args(cmd: Command) -> Command {
    let [iter, omega, rho, em_threshold, init] = em_args();
    cmd.arg(arg!(-s --source <SRC_FILE> "Source index file with reference sequences(.fmidx)")
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
        .arg(iter.default_value("100"))
        .arg(omega.default_value("1e-10"))
        .arg(rho.default_value("0"))
        .arg(arg!(-o --out <OUT_FILE>"Output prefix")
            .default_value("out")
            .value_parser(clap::value_parser!(String))
            )
        .arg(arg!(--eps_2 <EPS_2>"Minimum match log-probability threshold")
            .default_value("1e-18")
            .value_parser(clap::value_parser!(EMProb))
            )
        .arg(em_threshold.default_value("1e-6"))
        .arg(threads_arg())
        .arg(arg!(--seq2taxid <MAP_FILE> "seqid2taxid.map (seqid<whitespace>taxid) assigning index sequences to taxa: every sequence becomes a leaf under its taxid, an internal vertex the EM can report; sequences without a taxid hang under the root")
            .required(false)
            .value_parser(clap::value_parser!(String))
            )
        .arg(arg!(--nodes <NODES_FILE> "Taxonomy over the --seq2taxid taxids: NCBI nodes.dmp, or lines of `taxid parent [rank]`. The EM then runs hierarchically, pruning and reporting vertices of the tree")
            .required(false)
            .requires("seq2taxid")
            .value_parser(clap::value_parser!(String))
            )
        .arg(init.default_value("split"))
}

/// The `analyze` arguments: the report, the output prefix, which tables, and the EM
/// overrides that make it a refit.
fn analyze_args(cmd: Command) -> Command {
    let [iter, omega, rho, em_threshold, init] = em_args();
    cmd.arg(arg!(-r --report <REPORT> "Report written by `query` or an earlier `analyze` (.report or .report.zst)")
            .required(true)
            .value_parser(clap::value_parser!(String))
            )
        .arg(arg!(-o --out <OUT_FILE> "Output prefix: <OUT_FILE>.matches / .posteriors / .props / .nodes / .report")
            .required(true)
            .value_parser(clap::value_parser!(String))
            )
        .arg(gamma_arg())
        .arg(matches_arg("Write <OUT_FILE>.matches: the vertex each read is assigned to at --gamma"))
        .arg(Arg::new("posteriors")
            .long("posteriors")
            .help("Write <OUT_FILE>.posteriors: the posterior of every classified read on each EM leaf")
            .num_args(0)
            .action(ArgAction::SetTrue))
        .arg(Arg::new("props")
            .long("props")
            .help("Write <OUT_FILE>.props, the proportions as a legacy table")
            .num_args(0)
            .action(ArgAction::SetTrue))
        .arg(Arg::new("nodes")
            .long("nodes")
            .help("Write <OUT_FILE>.nodes, the fitted vertices as a legacy table")
            .num_args(0)
            .action(ArgAction::SetTrue))
        .arg(Arg::new("write-report")
            .long("write-report")
            .help("Write the refit as <OUT_FILE>.report (needs an EM parameter below)")
            .num_args(0)
            .action(ArgAction::SetTrue))
        .arg(Arg::new("compress")
            .long("compress")
            .help("zstd-compress the refit report (<OUT_FILE>.report.zst)")
            .num_args(0)
            .action(ArgAction::SetTrue)
            .requires("write-report"))
        .arg(threads_arg())
        .arg(iter.help("Refit: number of EM iterations (default: the report's)"))
        .arg(omega.help("Refit: penalty weight omega (default: the report's)"))
        .arg(rho.help("Refit: penalty weight rho; > 0 turns the penalty on, 0 turns it off (default: the report's)")
            .conflicts_with("no-penalty"))
        .arg(em_threshold.help("Refit: EM convergence threshold (default: the report's)"))
        .arg(init.help("Refit: initial conditionals, split or uniform (default: the report's)"))
        .arg(Arg::new("no-penalty")
            .long("no-penalty")
            .help("Refit without the penalty")
            .num_args(0)
            .action(ArgAction::SetTrue))
}

fn gamma_arg() -> Arg {
    arg!(--gamma <GAMMA>)
        .help("Confidence level in [0.5, 1), used with --seq2taxid: a read descends to the deepest vertex whose subtree posterior exceeds GAMMA. Without a taxonomy every read goes to its most probable sequence and GAMMA has no effect")
        .default_value("0.5")
        .value_parser(clap::value_parser!(EMProb))
}

fn matches_arg(help: &'static str) -> Arg {
    Arg::new("matches")
        .long("matches")
        .help(help)
        .num_args(0)
        .action(ArgAction::SetTrue)
}

/// `--gamma`, validated.
fn parse_gamma(sub_m: &clap::ArgMatches) -> Result<EMProb> {
    let gamma = *sub_m.get_one::<EMProb>("gamma").expect("defaulted");
    if !(0.5..1.0).contains(&gamma) {
        anyhow::bail!("--gamma must be in [0.5, 1), got {}", gamma);
    }
    Ok(gamma)
}

/// The [`AnalyzeParams`] of a parsed `query` command line.
fn analyze_params(sub_m: &clap::ArgMatches) -> AnalyzeParams<'_> {
    AnalyzeParams {
        ref_file: sub_m
            .get_one::<String>("source")
            .expect("required")
            .as_str(),
        r1_file: sub_m.get_one::<String>("r1").expect("required").as_str(),
        r2_file: sub_m.get_one::<String>("r2").expect("required").as_str(),
        mem_seed_length: *sub_m.get_one::<usize>("mem").expect("defaulted"),
        eps_1: *sub_m.get_one::<EMProb>("eps_1").expect("defaulted"),
        eps_2: *sub_m.get_one::<EMProb>("eps_2").expect("defaulted"),
        num_iter: *sub_m.get_one::<usize>("iter").expect("defaulted"),
        rho: *sub_m.get_one::<EMProb>("rho").expect("defaulted"),
        omega: *sub_m.get_one::<EMProb>("omega").expect("defaulted"),
        em_threshold: *sub_m.get_one::<EMProb>("em_threshold").expect("defaulted"),
        use_penalty: *sub_m.get_one::<EMProb>("rho").expect("defaulted") > 0.0,
        threads: *sub_m.get_one::<usize>("threads").expect("defaulted"),
        seq2taxid: sub_m.get_one::<String>("seq2taxid").map(String::as_str),
        nodes: sub_m.get_one::<String>("nodes").map(String::as_str),
        init: Init::parse(sub_m.get_one::<String>("init").expect("defaulted"))
            .expect("clap restricts the values"),
        index_hash: !sub_m.get_flag("no-index-hash"),
    }
}

/// Delete the outputs of a run that failed midway, so no truncated file is mistaken for
/// a complete one.
fn remove_partial(paths: &[&str]) {
    for p in paths {
        let _ = std::fs::remove_file(p);
    }
}

/// `elapsed, peak RSS, report size` for the completion line.
fn run_summary(summary: &analyze::AnalyzeSummary, report_path: Option<&str>) -> String {
    let human = |bytes: u64| -> String {
        let gib = bytes as f64 / (1u64 << 30) as f64;
        if gib >= 1.0 {
            format!("{:.2} GiB", gib)
        } else {
            format!("{:.1} MiB", bytes as f64 / (1u64 << 20) as f64)
        }
    };
    let size = report_path
        .and_then(|p| std::fs::metadata(p).ok())
        .map(|m| format!(", report {}", human(m.len())))
        .unwrap_or_default();
    let rss = summary
        .profile
        .peak_rss_kb
        .map(|kb| format!(", peak RSS {}", human(kb * 1024)))
        .unwrap_or_default();
    format!("{:.2}s{}{}", summary.profile.total_s, rss, size)
}

pub fn run() -> Result<()> {
    let matches = Command::new(env!("CARGO_PKG_NAME"))
        .display_name("PREMISE")
        .version(env!("CARGO_PKG_VERSION"))
        .propagate_version(true)
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
            query_args(Command::new("query")
                .about("Classify paired-end reads: align, fit the EM and write the report; with --matches also the .matches table"))
                .arg(gamma_arg())
                .arg(matches_arg("Write <OUT>.matches: the vertex each read is assigned to at --gamma"))
                .arg(Arg::new("no-report")
                    .long("no-report")
                    .help("Write only <OUT>.matches (needs --matches)")
                    .num_args(0)
                    .action(ArgAction::SetTrue)
                    .requires("matches"))
                .arg(Arg::new("compress")
                    .long("compress")
                    .help("zstd-compress the report (writes <OUT>.report.zst)")
                    .num_args(0)
                    .action(ArgAction::SetTrue)
                    .conflicts_with("no-report"))
                .arg(Arg::new("no-index-hash")
                    .long("no-index-hash")
                    .help("Skip the SHA-256 of the index in the report (a full read of the index file, run alongside the alignment)")
                    .num_args(0)
                    .action(ArgAction::SetTrue))
        )
        .subcommand(analyze_args(
            Command::new("analyze")
                .about("From a report, without the index or the reads: the legacy tables at a gamma, or a refit of the EM with other penalty parameters"),
        ))
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
        .about("PREMISE: A probabilistic framework for source assignment of Illumina reads")
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
            let outfile = sub_m.get_one::<String>("out").unwrap().as_str();
            let gamma = parse_gamma(sub_m)?;
            let write_matches = sub_m.get_flag("matches");
            let no_report = sub_m.get_flag("no-report");
            let compress = sub_m.get_flag("compress");
            let params = analyze_params(sub_m);
            let report_path = if compress {
                format!("{}.report.zst", outfile)
            } else {
                format!("{}.report", outfile)
            };
            let matches_path = format!("{}.matches", outfile);
            let mut written: Vec<&str> = Vec::new();
            if !no_report {
                written.push(&report_path);
            }
            if write_matches {
                written.push(&matches_path);
            }
            let run = || -> Result<analyze::AnalyzeSummary> {
                let mut report = if no_report {
                    None
                } else {
                    Some(Sink::create(&report_path, compress)?)
                };
                let mut matches = match write_matches {
                    true => Some(
                        File::create(&matches_path)
                            .map(BufWriter::new)
                            .with_context(|| format!("failed to create '{}'", matches_path))?,
                    ),
                    false => None,
                };
                let summary = run_analyze(
                    &params,
                    None,
                    AnalyzeOutputs {
                        report: report.as_mut().map(|r| r as &mut dyn Write),
                        matches: matches.as_mut().map(|m| (m as &mut dyn Write, gamma)),
                    },
                )?;
                if let Some(r) = report {
                    r.finish()?;
                }
                if let Some(mut m) = matches {
                    m.flush()?;
                }
                Ok(summary)
            };
            let summary = run()
                .inspect_err(|_| remove_partial(&written))
                .with_context(|| format!("query run failed (output {})", written.join(", ")))?;
            println!(
                "Query written to {} ({})",
                written.join(", "),
                run_summary(&summary, if no_report { None } else { Some(&report_path) })
            );
        }
        Some(("analyze", sub_m)) => {
            let report_path = sub_m
                .get_one::<String>("report")
                .expect("required")
                .as_str();
            let req = TableRequest {
                out: sub_m.get_one::<String>("out").expect("required").as_str(),
                gamma: parse_gamma(sub_m)?,
                matches: sub_m.get_flag("matches"),
                posteriors: sub_m.get_flag("posteriors"),
                props: sub_m.get_flag("props"),
                nodes: sub_m.get_flag("nodes"),
                report: sub_m.get_flag("write-report"),
                compress: sub_m.get_flag("compress"),
            };
            let refit = RefitParams {
                rho: sub_m.get_one::<EMProb>("rho").copied(),
                omega: sub_m.get_one::<EMProb>("omega").copied(),
                iter: sub_m.get_one::<usize>("iter").copied(),
                em_threshold: sub_m.get_one::<EMProb>("em_threshold").copied(),
                init: sub_m
                    .get_one::<String>("init")
                    .map(|i| Init::parse(i).expect("clap restricts the values")),
                no_penalty: sub_m.get_flag("no-penalty"),
            };
            let threads = *sub_m.get_one::<usize>("threads").expect("defaulted");
            let now = Instant::now();
            let written = analyze_report(report_path, &req, &refit, threads)?;
            println!("Written {} ({:.2?})", written.join(", "), now.elapsed());
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
                    .unwrap_or_default()
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
            let total_len: usize = ids.iter().map(|id| io::ref_len(&fmidx, *id)).sum();
            println!("Index file: {}", index_file);
            println!("Number of references: {}", ids.len());
            println!("Total sequence length: {} bp", total_len);
            println!("References:");
            for ref_idx in ids {
                let header = fmidx.seq_header(ref_idx).unwrap_or("");
                let len = io::ref_len(&fmidx, ref_idx);
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
