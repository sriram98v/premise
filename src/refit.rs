//! `analyze`: everything that starts from a report instead of the index and the reads.
//! `.matches` at a `gamma` and the NEXUS `.props` from the stored fit, or from a refit of
//! the EM over the report's alignment likelihoods with other penalty parameters.

use crate::analyze::{basename, fit_head, log_fit, peak_rss_kb, read_candidates, same_file};
use crate::em::{estimate_proportions, CsrBuilder, EmParams, Init};
use crate::report::export::props_nexus;
use crate::report::nexus::{open_report, Sink};
use crate::report::posterior::{leaf_by_seq, pi_tilde_by_seq};
use crate::report::read::{RawRead, ReportReader};
use crate::report::tables::{Candidate, Classifier, Scratch};
use crate::report::write::ReportWriter;
use crate::report::{PhaseTimes, Profile, ReadRows, ReportHead, ReportMeta, RunParams};
use crate::tree::{children, Partition};
use crate::utils::{EMProb, ReadIdx};
use crate::{resolve_thread_count, SeqId};
use anyhow::{bail, Context, Result};
use chrono::Local;
use phylo::node::NodeID;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufWriter, Write};
use std::sync::Mutex;
use std::time::Instant;

/// EM parameters given on the `analyze` command line. Any of them asks for a refit;
/// the rest come from the report's own PARAMS.
#[derive(Clone, Copy, Debug, Default)]
pub struct RefitParams {
    pub rho: Option<EMProb>,
    pub omega: Option<EMProb>,
    pub iter: Option<usize>,
    pub em_threshold: Option<EMProb>,
    pub init: Option<Init>,
    pub no_penalty: bool,
}

impl RefitParams {
    pub fn is_refit(&self) -> bool {
        self.rho.is_some()
            || self.omega.is_some()
            || self.iter.is_some()
            || self.em_threshold.is_some()
            || self.init.is_some()
            || self.no_penalty
    }

    /// The report's parameters with these overrides, validated. `--rho` > 0 turns the
    /// penalty on and `--rho 0` or `--no-penalty` off; the alignment parameters are the
    /// report's.
    fn apply(&self, base: &RunParams) -> Result<RunParams> {
        let mut p = base.clone();
        if let Some(rho) = self.rho {
            p.rho = rho;
            p.penalty = rho > 0.0;
        }
        if self.no_penalty {
            p.penalty = false;
        }
        if let Some(omega) = self.omega {
            p.omega = omega;
        }
        if let Some(iter) = self.iter {
            p.iter = iter;
        }
        if let Some(t) = self.em_threshold {
            p.em_threshold = t;
        }
        if let Some(init) = self.init {
            p.init = init.name().to_string();
        }
        let nonneg = |v: EMProb, what: &str| -> Result<()> {
            if !(v.is_finite() && v >= 0.0) {
                bail!("{} must be a finite non-negative number, got {}", what, v);
            }
            Ok(())
        };
        nonneg(p.rho, "--rho")?;
        nonneg(p.omega, "--omega")?;
        nonneg(p.em_threshold, "--em_threshold")?;
        if p.iter == 0 {
            bail!("--iter must be at least 1");
        }
        Ok(p)
    }
}

/// What `analyze` writes, all as `<out>.<ext>`.
#[derive(Clone, Copy, Debug)]
pub struct TableRequest<'a> {
    pub out: &'a str,
    /// The confidence level `.matches` classifies reads at.
    pub gamma: EMProb,
    pub matches: bool,
    pub props: bool,
    /// The refit report.
    pub report: bool,
    pub compress: bool,
}

impl TableRequest<'_> {
    fn path(&self, ext: &str) -> String {
        format!("{}.{}", self.out, ext)
    }

    fn report_path(&self) -> String {
        self.path(if self.compress {
            "report.zst"
        } else {
            "report"
        })
    }

    /// Every file the request writes, in the order the completion line lists them.
    pub fn paths(&self) -> Vec<String> {
        let mut v = Vec::new();
        if self.report {
            v.push(self.report_path());
        }
        if self.matches {
            v.push(self.path("matches"));
        }
        if self.props {
            v.push(self.path("props"));
        }
        v
    }
}

type Reader = ReportReader<Box<dyn BufRead>>;

/// Reads of a parallel chunk, and the text a chunk may hold before it is parsed.
const CHUNK_READS: usize = 16_384;
const CHUNK_BYTES: usize = 64 << 20;

/// Every remaining read of ALIGNMENT_LIKELIHOODS, in chunks: reads are cut from the
/// stream in order, then parsed and handed to `map` on the rayon pool (with a per-thread
/// state from `init`), and `consume` sees each read and its result in file order, so the
/// outputs, and the first error, are the sequential ones. Returns the number of reads.
fn par_reads<S: Send, T: Send>(
    reader: &mut Reader,
    init: impl Fn() -> S,
    map: impl Fn(&mut S, &ReadRows) -> Result<T> + Sync,
    mut consume: impl FnMut(&ReadRows, T) -> Result<()>,
) -> Result<usize> {
    let states: Vec<Mutex<(S, Vec<EMProb>)>> = (0..rayon::current_num_threads())
        .map(|_| Mutex::new((init(), Vec::new())))
        .collect();
    let mut raws: Vec<RawRead> = Vec::new();
    let mut n = 0usize;
    loop {
        let (mut filled, mut bytes, mut done) = (0, 0, false);
        while filled < CHUNK_READS && bytes < CHUNK_BYTES {
            if filled == raws.len() {
                raws.push(RawRead::default());
            }
            if !reader.next_raw_read(&mut raws[filled])? {
                done = true;
                break;
            }
            bytes += raws[filled].text_len();
            filled += 1;
        }
        let parsed: Vec<Result<(ReadRows, T)>> = raws[..filled]
            .par_iter()
            .map(|raw| {
                let slot = rayon::current_thread_index().unwrap_or(0) % states.len();
                let mut state = states[slot].lock().unwrap_or_else(|e| e.into_inner());
                let (s, palette) = &mut *state;
                let mut rows = ReadRows::default();
                raw.parse(&mut rows, palette)?;
                let t = map(s, &rows)?;
                Ok((rows, t))
            })
            .collect();
        for r in parsed {
            let (rows, t) = r?;
            consume(&rows, t)?;
            n += 1;
        }
        if done {
            return Ok(n);
        }
    }
}

/// Parse a report's head; the reader is left on its first read.
fn parse(report_path: &str) -> Result<(ReportHead, usize, Reader)> {
    let reader = open_report(report_path)?;
    ReportReader::parse_head(reader)
        .with_context(|| format!("failed to parse report '{}'", report_path))
}

/// Size and modification time: enough to notice the report changing between the two
/// passes of a refit.
fn fingerprint(report_path: &str) -> Result<(u64, Option<std::time::SystemTime>)> {
    let m = std::fs::metadata(report_path)
        .with_context(|| format!("failed to stat report '{}'", report_path))?;
    Ok((m.len(), m.modified().ok()))
}

/// Run `analyze`. Returns the files written; on an error none is left behind.
pub fn analyze_report(
    report_path: &str,
    req: &TableRequest<'_>,
    refit: &RefitParams,
    threads: usize,
) -> Result<Vec<String>> {
    let t_start = Instant::now();
    // Every pass parses on this pool. Fails when a pool already exists (a library
    // caller); its size is then kept.
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(resolve_thread_count(threads))
        .build_global();
    let paths = req.paths();
    if paths.is_empty() {
        bail!("analyze: nothing to do; give --matches, --props or --write-report");
    }
    if req.report && !refit.is_refit() {
        bail!(
            "analyze: --write-report without EM parameters would only copy the input; give --rho, --omega, --iter, --em_threshold, --init or --no-penalty to refit"
        );
    }
    if let Some(clobber) = paths.iter().find(|o| same_file(o, report_path)) {
        bail!("analyze: output '{}' is the input report", clobber);
    }
    let stamp = fingerprint(report_path)?;
    let (head0, n_reads, reader) = parse(report_path)?;
    let source = Source {
        path: report_path,
        n_reads,
        stamp,
    };
    let (head, reader) = if refit.is_refit() {
        let head = refit_from_report(report_path, &head0, n_reads, reader, refit, threads)?;
        // Pass 2 re-reads the file (it must be a regular file, not a pipe).
        (head, source.reopen()?)
    } else {
        (head0, reader)
    };
    let mut created: Vec<String> = Vec::new();
    let result = write_outputs(&head, &source, reader, req, t_start, &mut created);
    if result.is_err() {
        // Never leave a half-written file behind (only what this run opened).
        for p in &created {
            let _ = std::fs::remove_file(p);
        }
    }
    result.map(|()| paths)
}

/// Every sequence of a read must be a REFMAP column.
fn check_refmap(head: &ReportHead, buf: &ReadRows) -> Result<()> {
    if let Some(r) = buf
        .rows
        .iter()
        .find(|r| head.refmap.tip_of(r.seq).is_none())
    {
        bail!(
            "report: sequence {} of read '{}' is not in REFMAP",
            r.seq.index(),
            buf.read_id
        );
    }
    Ok(())
}

fn check_count(n: usize, n_reads: usize) -> Result<()> {
    if n != n_reads {
        bail!("report announces {} reads but holds {}", n_reads, n);
    }
    Ok(())
}

/// Pass 1 of a refit: the likelihood matrix from the report, the EM over the stored
/// tree, and the new head.
fn refit_from_report(
    report_path: &str,
    head0: &ReportHead,
    n_reads: usize,
    mut reader: Reader,
    refit: &RefitParams,
    threads: usize,
) -> Result<ReportHead> {
    let t_load = Instant::now();
    let num_threads = resolve_thread_count(threads);
    // Fails when a pool already exists (a library caller); its size is then kept.
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global();
    let params = refit.apply(&head0.meta.params)?;
    let init = Init::parse(&params.init)
        .ok_or_else(|| anyhow::anyhow!("report PARAMS: unknown init '{}'", params.init))?;
    let rho_eff = if params.penalty { params.rho } else { 0.0 };
    let timestamp = Local::now().format("%Y-%m-%d %H:%M:%S").to_string();
    println!(
        "Refit of {}\nTimestamp: {}\nNum Threads: {}\nEM Iterations: {}\nEM Threshold: {:e}\nRho: {}\nOmega: {:e}\nInit: {}",
        basename(report_path),
        timestamp,
        num_threads,
        params.iter,
        params.em_threshold,
        rho_eff,
        params.omega,
        init.name()
    );

    let mut builder = CsrBuilder::new(n_reads);
    let mut n = 0usize;
    par_reads(
        &mut reader,
        || (),
        |_, _| Ok(()),
        |rows, ()| {
            check_refmap(head0, rows)?;
            builder.push(ReadIdx(n), rows.rows.iter().map(|r| (r.seq, r.loglik)));
            n += 1;
            Ok(())
        },
    )?;
    check_count(n, n_reads)?;
    drop(reader);
    let csr = builder.finish();
    let refmap_seqs: Vec<SeqId> = head0.refmap.entries.iter().map(|e| e.0).collect();
    if csr.refs() != refmap_seqs.as_slice() {
        let odd = refmap_seqs
            .iter()
            .find(|s| !csr.refs().contains(s))
            .or_else(|| csr.refs().iter().find(|s| !refmap_seqs.contains(s)));
        bail!(
            "report: REFMAP lists {} sequences but the alignment rows use {} (sequence {} is in one only)",
            refmap_seqs.len(),
            csr.refs().len(),
            odd.map_or(0, |s| s.index())
        );
    }
    let load_s = t_load.elapsed().as_secs_f64();

    let t_em = Instant::now();
    let mut tree = head0.tree.to_stored_tax_tree()?;
    let stored = tree.clone();
    let tip_of_col: Vec<NodeID> = head0.refmap.entries.iter().map(|e| e.1).collect();
    if let Some((seq, tip)) = head0
        .refmap
        .entries
        .iter()
        .find(|(_, tip)| !children(&stored, *tip).is_empty())
    {
        bail!(
            "report: REFMAP puts sequence {} at vertex {}, which is not a tip of the stored tree",
            seq.index(),
            tip
        );
    }
    let mut part = Partition::new(tip_of_col.clone());
    let mut fit = estimate_proportions(
        &csr,
        &mut tree,
        &mut part,
        &EmParams {
            init,
            num_iter: params.iter,
            rho: rho_eff,
            omega: params.omega,
            em_threshold: params.em_threshold,
        },
        None,
    );
    log_fit(&fit, &tree, rho_eff);
    let em_s = t_em.elapsed().as_secs_f64();

    let meta = ReportMeta {
        timestamp,
        source: Some(basename(report_path)),
        params,
        threads: num_threads,
        runtime: PhaseTimes {
            load_s,
            align_s: 0.0,
            em_s,
        },
        ..head0.meta.clone()
    };
    Ok(fit_head(
        meta,
        &stored,
        &tip_of_col,
        &tree,
        &part,
        &csr,
        &mut fit,
    ))
}

/// Open an output for writing, remembering it for cleanup.
fn create(path: String, created: &mut Vec<String>) -> Result<BufWriter<File>> {
    created.push(path.clone());
    File::create(&path)
        .map(BufWriter::new)
        .with_context(|| format!("failed to create '{}'", path))
}

/// The report being analysed: where it is, its read count and its fingerprint, so a
/// later pass can re-read it and notice a change.
struct Source<'a> {
    path: &'a str,
    n_reads: usize,
    stamp: (u64, Option<std::time::SystemTime>),
}

impl Source<'_> {
    /// A fresh reader on the first read of ALIGNMENT_LIKELIHOODS.
    fn reopen(&self) -> Result<Reader> {
        let (_, n_again, reader) = parse(self.path)?;
        if n_again != self.n_reads || fingerprint(self.path)? != self.stamp {
            bail!("report '{}' changed while it was being analysed", self.path);
        }
        Ok(reader)
    }
}

/// The tables from `head`: `.props` first, then one stream over the reads for the
/// report copy and `.matches`. `created` lists every file opened.
fn write_outputs(
    head: &ReportHead,
    source: &Source<'_>,
    mut reader: Reader,
    req: &TableRequest<'_>,
    t_start: Instant,
    created: &mut Vec<String>,
) -> Result<()> {
    let n_reads = source.n_reads;
    if req.props {
        let path = req.path("props");
        created.push(path.clone());
        std::fs::write(&path, props_nexus(head))
            .with_context(|| format!("failed to write '{}'", path))?;
    }
    if !req.matches && !req.report {
        return Ok(());
    }
    let t_write = Instant::now();
    let mut report = match req.report {
        true => {
            let path = req.report_path();
            created.push(path.clone());
            Some(ReportWriter::begin(
                Sink::create(&path, req.compress)?,
                head,
                n_reads,
            )?)
        }
        false => None,
    };
    let classifier = Classifier::new(head.tree.to_tax_tree(), req.gamma);
    let mut matches = match req.matches {
        true => {
            let mut out = create(req.path("matches"), created)?;
            classifier.write_header(&mut out, &head.meta.em_summary(), &head.meta.label_header)?;
            Some(out)
        }
        false => None,
    };
    let pi_tilde = pi_tilde_by_seq(head);
    let leaf_of = leaf_by_seq(head);
    let want_matches = matches.is_some();
    // Each read's `.matches` row is classified and formatted in parallel, then written
    // in order.
    let n = par_reads(
        &mut reader,
        || {
            (
                Scratch::new(&classifier),
                Vec::<EMProb>::new(),
                Vec::<Candidate>::new(),
            )
        },
        |(scratch, posts, cands), rows| {
            let mut line = Vec::new();
            if want_matches {
                read_candidates(&rows.rows, &pi_tilde, &leaf_of, posts, cands);
                classifier.write_row(scratch, &mut line, &rows.read_id, cands)?;
            }
            Ok(line)
        },
        |rows, line| {
            check_refmap(head, rows)?;
            if let Some(w) = report.as_mut() {
                w.push_read(&rows.read_id, &rows.rows)?;
            }
            if let Some(m) = matches.as_mut() {
                m.write_all(&line)?;
            }
            Ok(())
        },
    )?;
    check_count(n, n_reads)?;
    if let Some(mut m) = matches {
        m.flush()?;
    }
    if let Some(mut w) = report {
        w.end_alignments()?;
        let sink = w.finish(&Profile {
            write_s: t_write.elapsed().as_secs_f64(),
            total_s: t_start.elapsed().as_secs_f64(),
            peak_rss_kb: peak_rss_kb(),
        })?;
        sink.finish()?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn base() -> RunParams {
        RunParams {
            mem: 22,
            eps_1: 0.0,
            eps_2: 1e-18,
            iter: 100,
            rho: 150.0,
            omega: 1e-10,
            em_threshold: 1e-6,
            init: "split".into(),
            penalty: false,
        }
    }

    #[test]
    fn overrides_fall_back_to_the_report_and_rho_enables_the_penalty() {
        let none = RefitParams::default();
        assert!(!none.is_refit());
        assert_eq!(none.apply(&base()).unwrap(), base());
        let p = RefitParams {
            rho: Some(8.0),
            init: Some(Init::Uniform),
            ..Default::default()
        };
        assert!(p.is_refit());
        let got = p.apply(&base()).unwrap();
        assert_eq!(
            (got.rho, got.penalty, got.init.as_str()),
            (8.0, true, "uniform")
        );
        assert_eq!((got.omega, got.iter, got.em_threshold), (1e-10, 100, 1e-6));
        let off = RefitParams {
            no_penalty: true,
            ..Default::default()
        };
        assert!(off.is_refit());
        let mut on = base();
        on.penalty = true;
        assert!(!off.apply(&on).unwrap().penalty);
        let zero = RefitParams {
            rho: Some(0.0),
            ..Default::default()
        };
        assert!(zero.is_refit());
        assert!(!zero.apply(&on).unwrap().penalty);
    }

    #[test]
    fn overrides_are_validated() {
        let bad = |p: RefitParams| p.apply(&base()).unwrap_err().to_string();
        assert!(bad(RefitParams {
            rho: Some(EMProb::NAN),
            ..Default::default()
        })
        .contains("--rho"));
        assert!(bad(RefitParams {
            omega: Some(-1.0),
            ..Default::default()
        })
        .contains("--omega"));
        assert!(bad(RefitParams {
            em_threshold: Some(EMProb::INFINITY),
            ..Default::default()
        })
        .contains("--em_threshold"));
        assert!(bad(RefitParams {
            iter: Some(0),
            ..Default::default()
        })
        .contains("--iter"));
        let mut b = base();
        b.rho = -5.0;
        assert!(RefitParams::default().apply(&b).is_err());
    }

    #[test]
    fn request_paths_follow_the_flags() {
        let req = TableRequest {
            out: "x",
            gamma: 0.6,
            matches: true,
            props: true,
            report: true,
            compress: true,
        };
        assert_eq!(req.paths(), vec!["x.report.zst", "x.matches", "x.props"]);
        let props = TableRequest {
            out: "x",
            gamma: 0.6,
            matches: false,
            props: true,
            report: false,
            compress: false,
        };
        assert_eq!(props.paths(), vec!["x.props"]);
        let none = TableRequest {
            props: false,
            ..props
        };
        assert!(none.paths().is_empty());
    }
}
