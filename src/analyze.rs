//! The classification pipeline: align, fit the EM, then stream the report (and, for
//! `query`, the `.matches` table) read by read. Everything that starts from a report
//! instead lives in [`crate::refit`].

use crate::em::EmFit;
use crate::em::{estimate_proportions, EmParams, Init};
use crate::group::Grouping;
use crate::io::{
    all_read_ids, load_fastq_forward, load_fastq_reverse, pair_reads, read_index_file,
};
use crate::report::posterior::{leaf_by_seq, pi_tilde_by_seq, row_posteriors};
use crate::report::tables::{Candidate, MatchesWriter};
use crate::report::write::ReportWriter;
use crate::report::{
    AlignRow, EmIter, EmStats, PhaseTimes, Profile, RefMap, ReportHead, ReportMeta, ReportTree,
    RunParams, TaxonomyStats, FORMAT_VERSION,
};
use crate::tree::{
    contract_to_observed, label, reference_tree, Partition, RefTree, TaxNodes, TaxTree, NO_NODE,
};
use crate::utils::{EMProb, QueryProgress, ReadAlignments, ReadID};
use crate::{query_fastq, resolve_thread_count, CsrLikelihood};
use anyhow::Result;
use bio::stats::LogProb;
use chrono::Local;
use phylo::node::NodeID;
use phylo::prelude::*;
use sha2::{Digest, Sha256};
use std::collections::BTreeSet;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::sync::Arc;
use std::time::Instant;

/// Everything `query` / `analyze` take on the command line except `gamma`.
#[derive(Clone, Debug)]
pub struct AnalyzeParams<'a> {
    pub ref_file: &'a str,
    pub r1_file: &'a str,
    pub r2_file: &'a str,
    pub mem_seed_length: usize,
    pub eps_1: EMProb,
    pub eps_2: EMProb,
    pub num_iter: usize,
    pub rho: EMProb,
    pub omega: EMProb,
    pub em_threshold: EMProb,
    pub use_penalty: bool,
    pub threads: usize,
    pub seq2taxid: Option<&'a str>,
    pub nodes: Option<&'a str>,
    pub init: Init,
    /// Record the SHA-256 of the index (a full read of the file on a background thread).
    pub index_hash: bool,
}

/// Where the run writes: the report, and / or `.matches` at a `gamma`.
pub struct AnalyzeOutputs<'a> {
    pub report: Option<&'a mut dyn Write>,
    pub matches: Option<(&'a mut dyn Write, EMProb)>,
}

/// The file name without its directory (the report never records a path); `-` when the
/// path has no file name (`..`, `/`).
pub(crate) fn basename(path: &str) -> String {
    Path::new(path)
        .file_name()
        .map_or_else(|| "-".to_string(), |n| n.to_string_lossy().into_owned())
}

/// SHA-256 of a file, streamed (never loaded whole), as lowercase hex.
fn sha256_file(path: &str) -> std::io::Result<String> {
    let mut file = File::open(path)?;
    let mut hasher = Sha256::new();
    std::io::copy(&mut file, &mut hasher)?;
    Ok(format!("{:x}", hasher.finalize()))
}

/// The hash the background thread produced, or `None` with a warning when it failed.
fn join_index_hash(handle: std::thread::JoinHandle<std::io::Result<String>>) -> Option<String> {
    match handle.join() {
        Ok(Ok(h)) => Some(h),
        Ok(Err(e)) => {
            eprintln!(
                "Warning: could not hash the index ({}); INDEX_SHA256 left empty",
                e
            );
            None
        }
        Err(_) => {
            eprintln!("Warning: the index hashing thread failed; INDEX_SHA256 left empty");
            None
        }
    }
}

/// True when two paths name the same file (so an output would clobber an input).
pub(crate) fn same_file(a: &str, b: &str) -> bool {
    match (std::fs::canonicalize(a), std::fs::canonicalize(b)) {
        (Ok(x), Ok(y)) => x == y,
        _ => a == b,
    }
}

/// What a run leaves behind in memory.
pub struct AnalyzeSummary {
    pub head: ReportHead,
    pub profile: Profile,
    /// The parameter log printed to stdout.
    pub log: String,
}

/// Peak resident set size of this process in KiB, where the platform reports it.
pub fn peak_rss_kb() -> Option<u64> {
    #[cfg(unix)]
    {
        let mut usage = std::mem::MaybeUninit::<libc::rusage>::uninit();
        // SAFETY: getrusage writes a complete rusage into the buffer on success.
        let rc = unsafe { libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()) };
        if rc != 0 {
            return None;
        }
        // SAFETY: rc == 0, so the struct is initialised.
        let usage = unsafe { usage.assume_init() };
        let max = u64::try_from(usage.ru_maxrss).ok()?;
        // Linux reports KiB, macOS bytes.
        Some(if cfg!(target_os = "macos") {
            max / 1024
        } else {
            max
        })
    }
    #[cfg(not(unix))]
    {
        None
    }
}

/// The taxonomy as `run_analyze` logs it.
fn taxonomy_log(ref_tree: &RefTree, nodes: bool, seq2taxid: bool) -> String {
    if nodes {
        format!(
            "\nTaxonomy: {} vertices, {} sequences at tips ({} taxids not in the nodes file, {} sequences without taxid)",
            ref_tree.n_nodes(),
            ref_tree.n_tips(),
            ref_tree.n_unknown,
            ref_tree.n_unmapped
        )
    } else if seq2taxid {
        format!(
            "\nTaxonomy: {} vertices, {} sequences at tips under their taxids ({} sequences without taxid)",
            ref_tree.n_nodes(),
            ref_tree.n_tips(),
            ref_tree.n_unmapped
        )
    } else {
        String::new()
    }
}

/// The PREMISE block of a `query` run; the fit's own statistics are filled in by
/// [`fit_head`].
#[allow(clippy::too_many_arguments)]
fn query_meta(
    p: &AnalyzeParams<'_>,
    grouping: &Grouping,
    ref_tree: &RefTree,
    stored: &TaxTree,
    num_threads: usize,
    timestamp: String,
    index_sha256: Option<String>,
    n_reads: usize,
    runtime: PhaseTimes,
) -> ReportMeta {
    ReportMeta {
        version: FORMAT_VERSION,
        tool: format!("premise {}", env!("CARGO_PKG_VERSION")),
        timestamp,
        label_header: grouping.label_header().to_string(),
        index: basename(p.ref_file),
        index_sha256,
        r1: basename(p.r1_file),
        r2: basename(p.r2_file),
        seq2taxid: p.seq2taxid.map(basename),
        nodes: p.nodes.map(basename),
        source: None,
        params: RunParams {
            mem: p.mem_seed_length,
            eps_1: p.eps_1,
            eps_2: p.eps_2,
            iter: p.num_iter,
            rho: p.rho,
            omega: p.omega,
            em_threshold: p.em_threshold,
            init: p.init.name().to_string(),
            penalty: p.use_penalty,
        },
        threads: num_threads,
        reads_total: n_reads,
        reads_aligned: 0,
        taxonomy: TaxonomyStats {
            ref_vertices: ref_tree.n_nodes(),
            ref_tips: ref_tree.n_tips(),
            unknown_taxids: ref_tree.n_unknown,
            unmapped_seqs: ref_tree.n_unmapped,
            contracted_vertices: stored.num_nodes(),
            contracted_tips: stored.get_leaf_ids().count(),
            final_tips: 0,
        },
        em: EmStats::default(),
        runtime,
    }
}

/// Assemble the head blocks of a fit; consumes the fit's trace. `meta` is complete but
/// for what the fit itself decides (`reads_aligned`, `final_tips`, the EM statistics).
#[allow(clippy::too_many_arguments)]
pub(crate) fn fit_head(
    mut meta: ReportMeta,
    stored: &TaxTree,
    tip_of_col: &[NodeID],
    tree: &TaxTree,
    part: &Partition,
    csr: &CsrLikelihood,
    fit: &mut EmFit,
) -> ReportHead {
    let mut props: Vec<(NodeID, &str, EMProb)> = tree
        .get_leaf_ids()
        .filter_map(|leaf| {
            let pi = fit.node_pi.get(leaf).copied().unwrap_or(0.0);
            (pi > 0.0).then(|| (leaf, label(tree, leaf), pi))
        })
        .collect();
    props.sort_by(|a, b| a.1.cmp(b.1));
    let props: Vec<(NodeID, EMProb)> = props.into_iter().map(|(n, _, p)| (n, p)).collect();
    let em_trace: Vec<EmIter> = fit
        .ll_trace
        .iter()
        .zip(std::mem::take(&mut fit.pi_trace))
        .map(|(ll, leaves)| EmIter { ll: *ll, leaves })
        .collect();
    let report_tree = ReportTree::from_fit(stored, tip_of_col, tree, part, fit);
    meta.reads_aligned = csr.n_reads();
    meta.taxonomy.final_tips = tree.get_leaf_ids().count();
    meta.em = EmStats {
        iterations: fit.ll_trace.len().saturating_sub(1),
        pruned: fit.pruned,
        collapsed: fit.collapsed,
        root_collapsed: fit.root_collapsed,
        final_ll: fit.ll_trace.last().copied().unwrap_or(EMProb::NAN),
    };
    ReportHead {
        meta,
        em_trace,
        tree: report_tree,
        refmap: RefMap::new(csr.refs(), tip_of_col),
        props,
    }
}

/// What the EM did, as `query` and a refit report it.
pub(crate) fn log_fit(fit: &EmFit, tree: &TaxTree, rho_eff: EMProb) {
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
}

/// One pass over the reads: every finite alignment to the report and, when wanted, the
/// `.matches` row of each read from posteriors recomputed exactly as `analyze` does.
fn stream_alignments<W: Write>(
    mut writer: Option<&mut ReportWriter<W>>,
    mut matches: Option<MatchesWriter<&mut dyn Write>>,
    read_ids: &BTreeSet<&ReadID>,
    aligns: &ReadAlignments<'_>,
    pi_tilde: &[EMProb],
    leaf_of: &[NodeID],
) -> Result<()> {
    let mut rows: Vec<AlignRow> = Vec::new();
    let mut posts: Vec<EMProb> = Vec::new();
    let mut cands: Vec<Candidate> = Vec::new();
    for read_id in read_ids {
        rows.clear();
        if let Some(read_idx) = aligns.read_idx(read_id) {
            let (refs, als) = aligns.row(read_idx);
            for (seq, al) in refs.iter().zip(als) {
                let loglik = *al.get_full_match_ll() as EMProb;
                if loglik.is_finite() {
                    rows.push(AlignRow {
                        seq: *seq,
                        loglik,
                        pos: al.get_pos(),
                    });
                }
            }
        }
        if let Some(w) = writer.as_mut() {
            w.push_read(read_id, &rows)?;
        }
        if let Some(m) = matches.as_mut() {
            read_candidates(&rows, pi_tilde, leaf_of, &mut posts, &mut cands);
            m.push(read_id, &cands)?;
        }
    }
    if let Some(w) = writer.as_mut() {
        w.end_alignments()?;
    }
    if let Some(m) = matches {
        m.finish()?;
    }
    Ok(())
}

/// The candidates of one read for [`MatchesWriter::push`]: its unpruned rows with the
/// posteriors from [`row_posteriors`].
pub(crate) fn read_candidates(
    rows: &[AlignRow],
    pi_tilde: &[EMProb],
    leaf_of: &[NodeID],
    posts: &mut Vec<EMProb>,
    cands: &mut Vec<Candidate>,
) {
    row_posteriors(
        |s| pi_tilde.get(s.index()).copied().unwrap_or(0.0),
        rows,
        posts,
    );
    cands.clear();
    for (r, &posterior) in rows.iter().zip(posts.iter()) {
        let leaf = leaf_of.get(r.seq.index()).copied().unwrap_or(NO_NODE);
        if leaf == NO_NODE {
            continue;
        }
        cands.push(Candidate {
            leaf,
            posterior,
            loglik: r.loglik,
            pos: r.pos,
        });
    }
}

/// Run the pipeline. Returns the head and profile so the caller can keep serving them
/// (the web UI) after the report is on disk.
pub fn run_analyze(
    p: &AnalyzeParams<'_>,
    progress: Option<Arc<QueryProgress>>,
    out: AnalyzeOutputs<'_>,
) -> Result<AnalyzeSummary> {
    let t_start = Instant::now();
    let num_threads = resolve_thread_count(p.threads);
    // Fails when a pool already exists (the server runs several queries): the first
    // pool's size is then kept, which is what the server wants.
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global();

    let forward = load_fastq_forward(p.r1_file)?;
    let reverse = load_fastq_reverse(p.r2_file)?;
    let (all_reads, unpaired) = pair_reads(forward, reverse, p.r1_file, p.r2_file)?;
    let read_ids = all_read_ids(&all_reads, &unpaired);

    let fmidx = read_index_file(p.ref_file)?;
    // The index is in the page cache now; hash it while the alignment runs.
    let index_hash = p.index_hash.then(|| {
        let path = p.ref_file.to_string();
        std::thread::spawn(move || sha256_file(&path))
    });
    let grouping = match p.seq2taxid {
        Some(path) => Grouping::from_map(&fmidx, path)?,
        None => Grouping::identity(&fmidx),
    };
    // The taxonomy is built before the expensive alignment so bad input fails early.
    let tax_nodes = p.nodes.map(TaxNodes::from_file).transpose()?;
    let ref_tree = reference_tree(&grouping, tax_nodes.as_ref())?;
    let tree_log = taxonomy_log(&ref_tree, tax_nodes.is_some(), p.seq2taxid.is_some());
    let rho_eff = if p.use_penalty { p.rho } else { 0.0 };
    let timestamp = Local::now().format("%Y-%m-%d %H:%M:%S").to_string();

    let log_str = format!("Timestamp: {}\nNum Threads: {}\nMin matching run: {}\nEps_1: {:e} ({:.2})\nEps_2: {:e} ({:.2})\nEM Iterations: {}\nEM Threshold: {:e}\nRho: {}\nOmega: {:e}\nInit: {}{}{}",
        timestamp,
        num_threads,
        p.mem_seed_length,
        p.eps_1,
        p.eps_1.ln(),
        p.eps_2,
        p.eps_2.ln(),
        p.num_iter,
        p.em_threshold,
        rho_eff,
        p.omega,
        p.init.name(),
        grouping.log_suffix(),
        tree_log,
    );
    println!("{}", log_str);
    let load_s = t_start.elapsed().as_secs_f64();

    let t_align = Instant::now();
    let (out_aligns, _all_refs) = query_fastq(
        &fmidx,
        &all_reads,
        p.mem_seed_length,
        LogProb(p.eps_1.ln()),
        LogProb(p.eps_2.ln()),
        progress.clone(),
    )?;
    let align_s = t_align.elapsed().as_secs_f64();

    let t_em = Instant::now();
    let csr = CsrLikelihood::build(&out_aligns);
    let (mut tree, mut part) = contract_to_observed(&ref_tree, csr.refs())?;
    // The contracted tree and its tips before the EM prunes it: the report's tree.
    let stored = tree.clone();
    let tip_of_col: Vec<NodeID> = (0..part.n_cols()).map(|c| part.leaf_of(c)).collect();
    if !tree_log.is_empty() {
        println!(
            "Contracted taxonomy: {} vertices, {} tips with reads",
            stored.num_nodes(),
            stored.get_leaf_ids().count()
        );
    }
    let mut fit = estimate_proportions(
        &csr,
        &mut tree,
        &mut part,
        &EmParams {
            init: p.init,
            num_iter: p.num_iter,
            rho: rho_eff,
            omega: p.omega,
            em_threshold: p.em_threshold,
        },
        progress.clone(),
    );
    log_fit(&fit, &tree, rho_eff);
    let em_s = t_em.elapsed().as_secs_f64();

    let index_sha256 = index_hash.and_then(join_index_hash);
    let meta = query_meta(
        p,
        &grouping,
        &ref_tree,
        &stored,
        num_threads,
        timestamp,
        index_sha256,
        read_ids.len(),
        PhaseTimes {
            load_s,
            align_s,
            em_s,
        },
    );
    let head = fit_head(meta, &stored, &tip_of_col, &tree, &part, &csr, &mut fit);

    let t_write = Instant::now();
    if let Some(ref pr) = progress {
        pr.phase.store(3, std::sync::atomic::Ordering::Relaxed);
    }
    let pi_tilde = pi_tilde_by_seq(&head);
    let leaf_of = leaf_by_seq(&head);
    let matches = match out.matches {
        Some((w, gamma)) => Some(MatchesWriter::new(
            w,
            grouping.label_header(),
            head.tree.to_tax_tree(),
            gamma,
        )?),
        None => None,
    };
    let mut writer = match out.report {
        Some(w) => Some(ReportWriter::begin(w, &head, read_ids.len())?),
        None => None,
    };
    stream_alignments(
        writer.as_mut(),
        matches,
        &read_ids,
        &out_aligns,
        &pi_tilde,
        &leaf_of,
    )?;
    let profile = Profile {
        write_s: t_write.elapsed().as_secs_f64(),
        total_s: t_start.elapsed().as_secs_f64(),
        peak_rss_kb: peak_rss_kb(),
    };
    if let Some(w) = writer {
        w.finish(&profile)?;
    }

    Ok(AnalyzeSummary {
        head,
        profile,
        log: log_str,
    })
}
