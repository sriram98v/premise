//! Parse a report: the head blocks into memory, then ALIGNMENT_LIKELIHOODS one read at a
//! time.

use super::nexus::tokens;
use super::FORMAT_VERSION;
use super::{
    AlignRow, EmIter, EmStats, PhaseTimes, ReadRows, RefMap, ReportHead, ReportMeta, ReportTree,
    ReportVertex, RunParams, TaxonomyStats, VertexStatus,
};
use crate::utils::EMProb;
use crate::SeqId;
use anyhow::{anyhow, bail, Context, Result};
use phylo::node::NodeID;
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, Read};

/// Upper bound on a capacity hint taken from the file, so a corrupt count cannot
/// trigger a huge allocation before the rows prove it right.
const MAX_PREALLOC: usize = 1 << 16;

/// Longest line the reader accepts (a read with a million alignments is ~25 MB).
const MAX_LINE: u64 = 256 << 20;

/// A finite number.
fn finite(v: EMProb, what: &str) -> Result<EMProb> {
    if v.is_finite() {
        Ok(v)
    } else {
        bail!("{}: expected a finite number, got {}", what, v)
    }
}

/// A finite number in `[0, 1]`.
fn unit(v: EMProb, what: &str) -> Result<EMProb> {
    if v.is_finite() && (0.0..=1.0).contains(&v) {
        Ok(v)
    } else {
        bail!("{}: expected a number in [0, 1], got {}", what, v)
    }
}

/// A finite, non-negative probability or likelihood.
fn finite_nonneg(v: EMProb, what: &str) -> Result<EMProb> {
    if v.is_finite() && v >= 0.0 {
        Ok(v)
    } else {
        bail!("{}: expected a finite non-negative number, got {}", what, v)
    }
}

/// `key=value` pairs of a tokenised command (`k`, `=`, `v` triples).
fn kv(toks: &[String]) -> Vec<(&str, &str)> {
    let mut out = Vec::new();
    let mut i = 0;
    while i + 2 < toks.len() {
        if toks[i + 1] == "=" {
            out.push((toks[i].as_str(), toks[i + 2].as_str()));
            i += 3;
        } else {
            i += 1;
        }
    }
    out
}

fn get<'a>(pairs: &[(&'a str, &'a str)], key: &str, what: &str) -> Result<&'a str> {
    pairs
        .iter()
        .find(|(k, _)| *k == key)
        .map(|(_, v)| *v)
        .ok_or_else(|| anyhow!("{}: missing {}", what, key))
}

fn num<T: std::str::FromStr>(s: &str, what: &str) -> Result<T>
where
    T::Err: std::fmt::Display,
{
    s.parse::<T>()
        .map_err(|e| anyhow!("{}: bad number '{}' ({})", what, s, e))
}

fn opt_num<T: std::str::FromStr>(s: &str, what: &str) -> Result<Option<T>>
where
    T::Err: std::fmt::Display,
{
    if s == "-" {
        Ok(None)
    } else {
        num(s, what).map(Some)
    }
}

fn opt_str(s: &str) -> Option<String> {
    (s != "-").then(|| s.to_string())
}

fn flag(s: &str, what: &str) -> Result<bool> {
    match s {
        "0" => Ok(false),
        "1" => Ok(true),
        _ => bail!("{}: expected 0 or 1, got '{}'", what, s),
    }
}

/// Builds a [`ReportMeta`] command by command.
#[derive(Default)]
struct MetaBuilder {
    version: Option<u32>,
    tool: Option<String>,
    timestamp: Option<String>,
    label_header: Option<String>,
    index: Option<String>,
    index_sha256: Option<String>,
    r1: Option<String>,
    r2: Option<String>,
    seq2taxid: Option<String>,
    nodes: Option<String>,
    source: Option<String>,
    params: Option<RunParams>,
    threads: Option<usize>,
    reads: Option<(usize, usize)>,
    taxonomy: Option<TaxonomyStats>,
    em: Option<EmStats>,
    runtime: Option<PhaseTimes>,
}

impl MetaBuilder {
    fn command(&mut self, toks: &[String]) -> Result<()> {
        let key = toks[0].to_ascii_uppercase();
        let arg = toks.get(1).map(String::as_str).unwrap_or("");
        let pairs = kv(&toks[1..]);
        match key.as_str() {
            "VERSION" => self.version = Some(num(arg, "VERSION")?),
            "TOOL" => self.tool = Some(arg.to_string()),
            "TIMESTAMP" => self.timestamp = Some(arg.to_string()),
            "LABEL_HEADER" => self.label_header = Some(arg.to_string()),
            "INDEX" => self.index = Some(arg.to_string()),
            "INDEX_SHA256" => {
                self.index_sha256 = match arg {
                    "-" => None,
                    h if h.len() == 64 && h.bytes().all(|b| b.is_ascii_hexdigit()) => {
                        Some(h.to_ascii_lowercase())
                    }
                    other => bail!(
                        "INDEX_SHA256: expected 64 hex digits or '-', got '{}'",
                        other
                    ),
                }
            }
            "R1" => self.r1 = Some(arg.to_string()),
            "R2" => self.r2 = Some(arg.to_string()),
            "SEQ2TAXID" => self.seq2taxid = opt_str(arg),
            "NODES" => self.nodes = opt_str(arg),
            "SOURCE" => self.source = opt_str(arg),
            "PARAMS" => {
                self.params = Some(RunParams {
                    mem: num(get(&pairs, "mem", "PARAMS")?, "mem")?,
                    eps_1: finite(num(get(&pairs, "eps_1", "PARAMS")?, "eps_1")?, "eps_1")?,
                    eps_2: finite(num(get(&pairs, "eps_2", "PARAMS")?, "eps_2")?, "eps_2")?,
                    iter: num(get(&pairs, "iter", "PARAMS")?, "iter")?,
                    rho: finite(num(get(&pairs, "rho", "PARAMS")?, "rho")?, "rho")?,
                    omega: finite(num(get(&pairs, "omega", "PARAMS")?, "omega")?, "omega")?,
                    em_threshold: finite(
                        num(get(&pairs, "em_threshold", "PARAMS")?, "em_threshold")?,
                        "em_threshold",
                    )?,
                    init: get(&pairs, "init", "PARAMS")?.to_string(),
                    penalty: flag(get(&pairs, "penalty", "PARAMS")?, "penalty")?,
                })
            }
            "THREADS" => self.threads = Some(num(arg, "THREADS")?),
            "READS" => {
                self.reads = Some((
                    num(get(&pairs, "total", "READS")?, "total")?,
                    num(get(&pairs, "aligned", "READS")?, "aligned")?,
                ))
            }
            "TAXONOMY" => {
                let g = |k: &str| -> Result<usize> { num(get(&pairs, k, "TAXONOMY")?, k) };
                self.taxonomy = Some(TaxonomyStats {
                    ref_vertices: g("ref_vertices")?,
                    ref_tips: g("ref_tips")?,
                    unknown_taxids: g("unknown_taxids")?,
                    unmapped_seqs: g("unmapped_seqs")?,
                    contracted_vertices: g("contracted_vertices")?,
                    contracted_tips: g("contracted_tips")?,
                    final_tips: g("final_tips")?,
                })
            }
            "EM" => {
                self.em = Some(EmStats {
                    iterations: num(get(&pairs, "iterations", "EM")?, "iterations")?,
                    pruned: num(get(&pairs, "pruned", "EM")?, "pruned")?,
                    collapsed: num(get(&pairs, "collapsed", "EM")?, "collapsed")?,
                    root_collapsed: flag(get(&pairs, "root_collapsed", "EM")?, "root_collapsed")?,
                    final_ll: match get(&pairs, "final_ll", "EM")? {
                        "-" => EMProb::NAN,
                        v => finite(num(v, "final_ll")?, "final_ll")?,
                    },
                })
            }
            "RUNTIME" => {
                self.runtime = Some(PhaseTimes {
                    load_s: finite(num(get(&pairs, "load_s", "RUNTIME")?, "load_s")?, "load_s")?,
                    align_s: finite(
                        num(get(&pairs, "align_s", "RUNTIME")?, "align_s")?,
                        "align_s",
                    )?,
                    em_s: finite(num(get(&pairs, "em_s", "RUNTIME")?, "em_s")?, "em_s")?,
                })
            }
            // Later format versions may add commands; they carry nothing we need.
            _ => {}
        }
        Ok(())
    }

    fn build(self) -> Result<ReportMeta> {
        let need = |what: &str| anyhow!("PREMISE block: missing {}", what);
        let (reads_total, reads_aligned) = self.reads.ok_or_else(|| need("READS"))?;
        let version = self.version.ok_or_else(|| need("VERSION"))?;
        if version != FORMAT_VERSION {
            bail!(
                "report format version {} is not the version this PREMISE reads ({})",
                version,
                FORMAT_VERSION
            );
        }
        Ok(ReportMeta {
            version,
            tool: self.tool.ok_or_else(|| need("TOOL"))?,
            timestamp: self.timestamp.ok_or_else(|| need("TIMESTAMP"))?,
            label_header: self.label_header.ok_or_else(|| need("LABEL_HEADER"))?,
            index: self.index.ok_or_else(|| need("INDEX"))?,
            index_sha256: self.index_sha256,
            r1: self.r1.ok_or_else(|| need("R1"))?,
            r2: self.r2.ok_or_else(|| need("R2"))?,
            seq2taxid: self.seq2taxid,
            nodes: self.nodes,
            source: self.source,
            params: self.params.ok_or_else(|| need("PARAMS"))?,
            threads: self.threads.ok_or_else(|| need("THREADS"))?,
            reads_total,
            reads_aligned,
            taxonomy: self.taxonomy.ok_or_else(|| need("TAXONOMY"))?,
            em: self.em.ok_or_else(|| need("EM"))?,
            runtime: self.runtime.ok_or_else(|| need("RUNTIME"))?,
        })
    }
}

/// Line reader that tokenises and skips blank lines.
pub struct ReportReader<R: BufRead> {
    r: R,
    line: String,
    line_no: usize,
    /// The current read's distinct log-likelihoods.
    palette: Vec<EMProb>,
}

impl<R: BufRead> ReportReader<R> {
    /// Read the next line into `self.line`; `false` at end of input. Lines longer than
    /// [`MAX_LINE`] are an error rather than an allocation.
    fn read_line(&mut self) -> Result<bool> {
        self.line.clear();
        self.line_no += 1;
        let n = self.r.by_ref().take(MAX_LINE).read_line(&mut self.line)?;
        if n == 0 {
            return Ok(false);
        }
        if n as u64 >= MAX_LINE && !self.line.ends_with('\n') {
            bail!("line {} is longer than {} bytes", self.line_no, MAX_LINE);
        }
        Ok(true)
    }

    /// Next non-blank line as tokens; `None` at end of input.
    fn next_tokens(&mut self) -> Result<Option<Vec<String>>> {
        loop {
            if !self.read_line()? {
                return Ok(None);
            }
            let toks = tokens(&self.line).with_context(|| format!("line {}", self.line_no))?;
            if !toks.is_empty() {
                return Ok(Some(toks));
            }
        }
    }

    fn expect_tokens(&mut self, what: &str) -> Result<Vec<String>> {
        self.next_tokens()?
            .ok_or_else(|| anyhow!("unexpected end of report while reading {}", what))
    }

    fn is_end(toks: &[String]) -> bool {
        toks.len() == 2 && toks[0].eq_ignore_ascii_case("END") && toks[1] == ";"
    }

    fn is_matrix_end(toks: &[String]) -> bool {
        toks.len() == 1 && toks[0] == ";"
    }

    /// Skip to the end of the current block.
    fn skip_block(&mut self, name: &str) -> Result<()> {
        loop {
            let toks = self.expect_tokens(name)?;
            if Self::is_end(&toks) {
                return Ok(());
            }
        }
    }

    /// `DIMENSIONS key=n;` then `MATRIX`: returns `n`.
    fn matrix_header(&mut self, block: &str, key: &str) -> Result<usize> {
        let toks = self.expect_tokens(block)?;
        if !toks[0].eq_ignore_ascii_case("DIMENSIONS") {
            bail!("{} block: expected DIMENSIONS, got '{}'", block, toks[0]);
        }
        let n: usize = num(get(&kv(&toks[1..]), key, block)?, key)?;
        let toks = self.expect_tokens(block)?;
        if !toks[0].eq_ignore_ascii_case("MATRIX") {
            bail!("{} block: expected MATRIX, got '{}'", block, toks[0]);
        }
        Ok(n)
    }

    /// Rows of a MATRIX up to its `;`, then the block's `END;`. Never holds more than
    /// `max_rows` rows: a count larger than announced is an error, not an allocation.
    fn matrix_rows(&mut self, block: &str, max_rows: usize) -> Result<Vec<Vec<String>>> {
        let mut rows = Vec::new();
        loop {
            let toks = self.expect_tokens(block)?;
            if Self::is_matrix_end(&toks) {
                break;
            }
            if rows.len() >= max_rows {
                bail!("{}: more than {} rows", block, max_rows);
            }
            rows.push(toks);
        }
        let toks = self.expect_tokens(block)?;
        if !Self::is_end(&toks) {
            bail!("{} block: expected END; after the matrix", block);
        }
        Ok(rows)
    }

    fn premise_block(&mut self) -> Result<ReportMeta> {
        let mut b = MetaBuilder::default();
        loop {
            let toks = self.expect_tokens("PREMISE")?;
            if Self::is_end(&toks) {
                break;
            }
            let body: Vec<String> = toks.into_iter().filter(|t| t != ";").collect();
            if !body.is_empty() {
                b.command(&body)?;
            }
        }
        b.build()
    }

    fn em_trace_block(&mut self) -> Result<Vec<EmIter>> {
        let n = self.matrix_header("EM_TRACE", "NITER")?;
        // One header row per iteration plus its leaf rows; bound by a generous multiple.
        let rows = self.matrix_rows("EM_TRACE", n.saturating_mul(1 << 20).max(n))?;
        let mut trace = Vec::with_capacity(n.min(MAX_PREALLOC));
        let mut i = 0;
        while i < rows.len() {
            let r = &rows[i];
            if r.len() != 3 {
                bail!("EM_TRACE: expected 'iter ll n_leaves', got {:?}", r);
            }
            let ll: EMProb = finite(num(&r[1], "EM_TRACE ll")?, "EM_TRACE ll")?;
            let k: usize = num(&r[2], "EM_TRACE n_leaves")?;
            let mut leaves = Vec::with_capacity(k.min(rows.len()));
            for j in 1..=k {
                let row = rows
                    .get(i + j)
                    .ok_or_else(|| anyhow!("EM_TRACE: truncated iteration {}", r[0]))?;
                if row.len() != 2 {
                    bail!("EM_TRACE: expected 'node proportion', got {:?}", row);
                }
                leaves.push((
                    num(&row[0], "EM_TRACE node")?,
                    finite_nonneg(num(&row[1], "EM_TRACE pi")?, "EM_TRACE pi")?,
                ));
            }
            trace.push(EmIter { ll, leaves });
            i += 1 + k;
        }
        if trace.len() != n {
            bail!("EM_TRACE: {} iterations, {} announced", trace.len(), n);
        }
        Ok(trace)
    }

    fn nodes_block(&mut self) -> Result<ReportTree> {
        let n = self.matrix_header("NODES", "N")?;
        let rows = self.matrix_rows("NODES", n)?;
        if rows.len() != n {
            bail!("NODES: {} rows, {} announced", rows.len(), n);
        }
        let mut vertices = Vec::with_capacity(n);
        let mut status_of: HashMap<NodeID, VertexStatus> = HashMap::with_capacity(n);
        for r in rows {
            if r.len() != 8 {
                bail!("NODES: expected 8 columns, got {:?}", r);
            }
            let id: NodeID = num(&r[0], "NODES id")?;
            let parent: Option<NodeID> = opt_num(&r[1], "NODES parent")?;
            let status = VertexStatus::parse(&r[3])
                .ok_or_else(|| anyhow!("NODES: unknown status '{}'", r[3]))?;
            // The rows are a pre-order: the root first, every parent before its children,
            // and a vertex of the final tree only under another one.
            match parent {
                None if !vertices.is_empty() => bail!("NODES: vertex {} has no parent", id),
                Some(p) if vertices.is_empty() => {
                    bail!(
                        "NODES: the first vertex ({}) must be the root, not a child of {}",
                        id,
                        p
                    )
                }
                Some(p) if p == id => bail!("NODES: vertex {} is its own parent", id),
                Some(p) => match status_of.get(&p) {
                    None => bail!("NODES: vertex {} comes before its parent {}", id, p),
                    Some(ps) if status.is_final() && !ps.is_final() => bail!(
                        "NODES: {} vertex {} under {} vertex {}",
                        status.as_str(),
                        id,
                        ps.as_str(),
                        p
                    ),
                    Some(_) => {}
                },
                None => {}
            }
            if status_of.insert(id, status).is_some() {
                bail!("NODES: duplicate vertex id {}", id);
            }
            vertices.push(ReportVertex {
                id,
                parent,
                // Filled from TAXLABELS once both blocks are in.
                label: String::new(),
                rank: r[2].clone(),
                status,
                contigs: num(&r[4], "NODES contigs")?,
                count: opt_num(&r[5], "NODES count")?
                    .map(|c| finite_nonneg(c, "NODES count"))
                    .transpose()?,
                proportion: finite_nonneg(num(&r[6], "NODES proportion")?, "NODES proportion")?,
                conditional: opt_num(&r[7], "NODES conditional")?
                    .map(|q| unit(q, "NODES conditional"))
                    .transpose()?,
            });
        }
        Ok(ReportTree { vertices })
    }

    /// `DIMENSIONS NTAX=n;` then `TAXLABELS l1 l2 ... ;` (labels may span lines).
    fn taxa_block(&mut self) -> Result<Vec<String>> {
        let toks = self.expect_tokens("TAXA")?;
        if !toks[0].eq_ignore_ascii_case("DIMENSIONS") {
            bail!("TAXA block: expected DIMENSIONS, got '{}'", toks[0]);
        }
        let n: usize = num(get(&kv(&toks[1..]), "NTAX", "TAXA")?, "NTAX")?;
        let mut labels: Vec<String> = Vec::with_capacity(n.min(MAX_PREALLOC));
        let mut toks = self.expect_tokens("TAXA")?;
        if !toks[0].eq_ignore_ascii_case("TAXLABELS") {
            bail!("TAXA block: expected TAXLABELS, got '{}'", toks[0]);
        }
        toks.remove(0);
        loop {
            let mut done = false;
            for t in toks {
                if t == ";" {
                    done = true;
                    break;
                }
                if labels.len() >= n {
                    bail!("TAXA: more than {} labels", n);
                }
                labels.push(t);
            }
            if done {
                break;
            }
            toks = self.expect_tokens("TAXA")?;
        }
        let end = self.expect_tokens("TAXA")?;
        if !Self::is_end(&end) {
            bail!("TAXA block: expected END; after TAXLABELS");
        }
        if labels.len() != n {
            bail!("TAXA: {} labels, {} announced", labels.len(), n);
        }
        Ok(labels)
    }

    fn refmap_block(&mut self) -> Result<RefMap> {
        let n = self.matrix_header("REFMAP", "N")?;
        let rows = self.matrix_rows("REFMAP", n)?;
        if rows.len() != n {
            bail!("REFMAP: {} rows, {} announced", rows.len(), n);
        }
        let mut entries = Vec::with_capacity(n);
        for r in rows {
            if r.len() != 2 {
                bail!("REFMAP: expected 'seqid node', got {:?}", r);
            }
            let seq: u32 = num(&r[0], "REFMAP seqid")?;
            entries.push((SeqId::new(seq), num::<NodeID>(&r[1], "REFMAP node")?));
        }
        entries.sort_unstable_by_key(|e| e.0.index());
        if let Some(w) = entries.windows(2).find(|w| w[0].0 == w[1].0) {
            bail!("REFMAP: sequence {} listed twice", w[0].0.index());
        }
        Ok(RefMap { entries })
    }

    fn props_block(&mut self) -> Result<Vec<(NodeID, EMProb)>> {
        let n = self.matrix_header("PROPS", "N")?;
        let rows = self.matrix_rows("PROPS", n)?;
        if rows.len() != n {
            bail!("PROPS: {} rows, {} announced", rows.len(), n);
        }
        rows.into_iter()
            .map(|r| {
                if r.len() != 2 {
                    bail!("PROPS: expected 'node proportion', got {:?}", r);
                }
                Ok((
                    num(&r[0], "PROPS node")?,
                    unit(num(&r[1], "PROPS pi")?, "PROPS pi")?,
                ))
            })
            .collect()
    }

    /// Parse every block before ALIGNMENT_LIKELIHOODS and stop right after its `MATRIX`
    /// line. Returns the announced read count with the reader positioned on the first read.
    pub fn parse_head(r: R) -> Result<(ReportHead, usize, Self)> {
        let mut me = Self {
            r,
            line: String::new(),
            line_no: 0,
            palette: Vec::new(),
        };
        let first = me.expect_tokens("header")?;
        if first.len() != 1 || !first[0].eq_ignore_ascii_case("#NEXUS") {
            bail!("not a NEXUS report: first line is '{}'", me.line.trim_end());
        }
        let mut meta = None;
        let mut em_trace = None;
        let mut taxa: Option<Vec<String>> = None;
        let mut tree = None;
        let mut refmap = None;
        let mut props = None;
        loop {
            let Some(toks) = me.next_tokens()? else {
                bail!("report ends before an ALIGNMENT_LIKELIHOODS block");
            };
            if toks.len() != 3 || !toks[0].eq_ignore_ascii_case("BEGIN") || toks[2] != ";" {
                bail!(
                    "line {}: expected 'BEGIN <block>;', got {:?}",
                    me.line_no,
                    toks
                );
            }
            let name = toks[1].to_ascii_uppercase();
            match name.as_str() {
                "PREMISE" => meta = Some(me.premise_block()?),
                "EM_TRACE" => em_trace = Some(me.em_trace_block()?),
                "TAXA" => taxa = Some(me.taxa_block()?),
                "NODES" => tree = Some(me.nodes_block()?),
                "REFMAP" => refmap = Some(me.refmap_block()?),
                "PROPS" => props = Some(me.props_block()?),
                "ALIGNMENT_LIKELIHOODS" => {
                    let n = me.matrix_header("ALIGNMENT_LIKELIHOODS", "NREADS")?;
                    let need = |what: &str| anyhow!("report has no {} block", what);
                    let mut tree: ReportTree = tree.ok_or_else(|| need("NODES"))?;
                    let taxa = taxa.ok_or_else(|| need("TAXA"))?;
                    if taxa.len() != tree.vertices.len() {
                        bail!(
                            "TAXA lists {} labels but NODES has {} vertices",
                            taxa.len(),
                            tree.vertices.len()
                        );
                    }
                    for (v, label) in tree.vertices.iter_mut().zip(taxa) {
                        v.label = label;
                    }
                    let head = ReportHead {
                        meta: meta.ok_or_else(|| need("PREMISE"))?,
                        em_trace: em_trace.ok_or_else(|| need("EM_TRACE"))?,
                        tree,
                        refmap: refmap.ok_or_else(|| need("REFMAP"))?,
                        props: props.ok_or_else(|| need("PROPS"))?,
                    };
                    let known: HashSet<NodeID> = head.tree.vertices.iter().map(|v| v.id).collect();
                    if let Some((seq, node)) =
                        head.refmap.entries.iter().find(|(_, n)| !known.contains(n))
                    {
                        bail!(
                            "REFMAP: sequence {} maps to unknown vertex {}",
                            seq.index(),
                            node
                        );
                    }
                    let status_of: HashMap<NodeID, VertexStatus> = head
                        .tree
                        .vertices
                        .iter()
                        .map(|v| (v.id, v.status))
                        .collect();
                    let mut seen: HashSet<NodeID> = HashSet::new();
                    for (node, _) in &head.props {
                        match status_of.get(node) {
                            None => bail!("PROPS: unknown vertex {}", node),
                            Some(s) if !s.is_final() => {
                                bail!(
                                    "PROPS: vertex {} is {}, not a leaf of the final tree",
                                    node,
                                    s.as_str()
                                )
                            }
                            Some(_) => {}
                        }
                        if !seen.insert(*node) {
                            bail!("PROPS: vertex {} listed twice", node);
                        }
                    }
                    let total: EMProb = head.props.iter().map(|(_, p)| p).sum();
                    if !head.props.is_empty() && (total - 1.0).abs() > 1e-6 {
                        bail!("PROPS: proportions sum to {}, not 1", total);
                    }
                    return Ok((head, n, me));
                }
                _ => me.skip_block(&name)?,
            }
        }
    }

    /// One array line of a read: `n` whitespace tokens, each handed to `f` with its index.
    fn array_line(
        &mut self,
        n: usize,
        what: &str,
        read_id: &str,
        f: impl FnMut(usize, &str) -> Result<()>,
    ) -> Result<()> {
        self.array_line_of(n, false, what, read_id, f).map(|_| ())
    }

    /// [`Self::array_line`] that, with `allow_one`, also accepts a single token (a value
    /// shared by every row); returns the number of tokens.
    fn array_line_of(
        &mut self,
        n: usize,
        allow_one: bool,
        what: &str,
        read_id: &str,
        mut f: impl FnMut(usize, &str) -> Result<()>,
    ) -> Result<usize> {
        self.line.clear();
        self.line_no += 1;
        if self.r.read_line(&mut self.line)? == 0 {
            bail!(
                "ALIGNMENT_LIKELIHOODS: truncated {} of read '{}'",
                what,
                read_id
            );
        }
        let mut i = 0;
        for tok in self.line.split_ascii_whitespace() {
            if i == n {
                bail!(
                    "line {}: more than {} {} for read '{}'",
                    self.line_no,
                    n,
                    what,
                    read_id
                );
            }
            f(i, tok).with_context(|| {
                format!("line {}: {} of read '{}'", self.line_no, what, read_id)
            })?;
            i += 1;
        }
        if i != n && !(allow_one && i == 1) {
            bail!(
                "line {}: {} {} for read '{}', {} announced",
                self.line_no,
                i,
                what,
                read_id,
                n
            );
        }
        Ok(i)
    }

    /// The next read of ALIGNMENT_LIKELIHOODS into `buf`; `false` once the matrix is
    /// exhausted.
    pub fn next_read(&mut self, buf: &mut ReadRows) -> Result<bool> {
        // Read-id lines are `'id' n` or `id n`; the fast path skips the tokenizer when
        // nothing needs it.
        let toks: Vec<String> = loop {
            if !self.read_line()? {
                bail!("unexpected end of report while reading ALIGNMENT_LIKELIHOODS");
            }
            let t = self.line.trim();
            if t.is_empty() {
                continue;
            }
            if !t.contains('\'') && !t.contains('[') {
                break t.split_ascii_whitespace().map(str::to_string).collect();
            }
            let toks = tokens(t).with_context(|| format!("line {}", self.line_no))?;
            if !toks.is_empty() {
                break toks;
            }
        };
        if Self::is_matrix_end(&toks) {
            let end = self.expect_tokens("ALIGNMENT_LIKELIHOODS")?;
            if !Self::is_end(&end) {
                bail!("ALIGNMENT_LIKELIHOODS: expected END; after the matrix");
            }
            return Ok(false);
        }
        let n: usize = match toks.len() {
            2 | 3 => num(&toks[1], "ALIGNMENT_LIKELIHOODS n_rows")?,
            _ => bail!(
                "line {}: expected 'read_id 0' or 'read_id n_rows n_lls', got {:?}",
                self.line_no,
                toks
            ),
        };
        if (n == 0) != (toks.len() == 2) {
            bail!(
                "line {}: expected 'read_id 0' or 'read_id n_rows n_lls', got {:?}",
                self.line_no,
                toks
            );
        }
        buf.read_id.clear();
        buf.read_id.push_str(&toks[0]);
        buf.rows.clear();
        if n == 0 {
            return Ok(true);
        }
        let k: usize = num(&toks[2], "ALIGNMENT_LIKELIHOODS n_lls")?;
        if k == 0 || k > n {
            bail!(
                "line {}: read '{}' has {} distinct log-likelihoods for {} rows",
                self.line_no,
                toks[0],
                k,
                n
            );
        }
        buf.rows.reserve(n.min(MAX_PREALLOC));
        let ReadRows { read_id, rows } = buf;
        let read_id: &str = read_id;
        let mut prev: u32 = 0;
        self.array_line(n, "seqids", read_id, |i, tok| {
            let v: u32 = num(tok, "seqid")?;
            let seq = if i == 0 {
                v
            } else {
                if v == 0 {
                    bail!("seqid difference 0 (duplicate sequence)");
                }
                prev.checked_add(v)
                    .ok_or_else(|| anyhow!("seqid overflows u32"))?
            };
            prev = seq;
            rows.push(AlignRow {
                seq: SeqId::new(seq),
                loglik: 0.0,
                pos: (0, 0),
            });
            Ok(())
        })?;
        let mut palette = std::mem::take(&mut self.palette);
        palette.clear();
        let parsed = self
            .array_line(k, "log-likelihoods", read_id, |_, tok| {
                let ll: EMProb = num(tok, "log-likelihood")?;
                if !ll.is_finite() {
                    bail!("log-likelihood is not finite: {}", tok);
                }
                palette.push(ll);
                Ok(())
            })
            .and_then(|()| {
                if k == 1 {
                    rows.iter_mut().for_each(|r| r.loglik = palette[0]);
                    return Ok(());
                }
                self.array_line(n, "log-likelihood indices", read_id, |i, tok| {
                    let j: usize = num(tok, "log-likelihood index")?;
                    rows[i].loglik = *palette
                        .get(j)
                        .ok_or_else(|| anyhow!("log-likelihood index {} is not below {}", j, k))?;
                    Ok(())
                })
            });
        self.palette = palette;
        parsed?;
        self.array_line(n, "forward positions", read_id, |i, tok| {
            rows[i].pos.0 = num(tok, "forward position")?;
            Ok(())
        })?;
        let reverse = |fwd: usize, off: i64| {
            i64::try_from(fwd)
                .ok()
                .and_then(|f| f.checked_add(off))
                .and_then(|r| usize::try_from(r).ok())
                .ok_or_else(|| anyhow!("reverse position {} {:+} is out of range", fwd, off))
        };
        let mut shared: i64 = 0;
        let got = self.array_line_of(n, true, "insert offsets", read_id, |i, tok| {
            let off: i64 = num(tok, "insert offset")?;
            shared = off;
            rows[i].pos.1 = reverse(rows[i].pos.0, off)?;
            Ok(())
        })?;
        if got == 1 {
            for r in rows.iter_mut().skip(1) {
                r.pos.1 = reverse(r.pos.0, shared)
                    .with_context(|| format!("insert offset of read '{}'", read_id))?;
            }
        }
        Ok(true)
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::report::fixtures::fitted;
    use crate::report::nexus::Sink;
    use crate::report::write::ReportWriter;
    use crate::report::{Profile, ReportTree};
    use std::io::Read;

    pub(crate) fn sample_head() -> ReportHead {
        let (stored, tips, fit_tree, part, fit) = fitted();
        let tree = ReportTree::from_fit(&stored, &tips, &fit_tree, &part, &fit);
        ReportHead {
            meta: ReportMeta {
                version: FORMAT_VERSION,
                tool: concat!("PREMISE ", env!("CARGO_PKG_VERSION")).into(),
                timestamp: "2026-10-10 10:00:00".into(),
                label_header: "TaxID".into(),
                index: "idx with space.fmidx".into(),
                index_sha256: Some("a".repeat(64)),
                r1: "r1.fq.gz".into(),
                r2: "r2.fq.gz".into(),
                seq2taxid: Some("map 'quoted'".into()),
                nodes: None,
                source: None,
                params: RunParams {
                    mem: 22,
                    eps_1: 0.0,
                    eps_2: 1e-18,
                    iter: 100,
                    rho: 150.0,
                    omega: 1e-10,
                    em_threshold: 1e-6,
                    init: "split".into(),
                    penalty: true,
                },
                threads: 8,
                reads_total: 3,
                reads_aligned: 2,
                taxonomy: TaxonomyStats {
                    ref_vertices: 12,
                    ref_tips: 6,
                    unknown_taxids: 1,
                    unmapped_seqs: 0,
                    contracted_vertices: 10,
                    contracted_tips: 4,
                    final_tips: 2,
                },
                em: EmStats {
                    iterations: 1,
                    pruned: 3,
                    collapsed: 1,
                    root_collapsed: false,
                    final_ll: -8.0,
                },
                runtime: PhaseTimes {
                    load_s: 0.5,
                    align_s: 1.25,
                    em_s: 0.125,
                },
            },
            em_trace: fit
                .ll_trace
                .iter()
                .zip(&fit.pi_trace)
                .map(|(ll, leaves)| EmIter {
                    ll: *ll,
                    leaves: leaves.clone(),
                })
                .collect(),
            tree,
            // columns 0..4 (tips 6, 7, 8, 9) are sequences 4, 9, 2, 7
            refmap: RefMap::new(
                &[SeqId::new(4), SeqId::new(9), SeqId::new(2), SeqId::new(7)],
                &tips,
            ),
            props: vec![(3, 0.4), (8, 0.6)],
        }
    }

    fn row(seq: u32, loglik: EMProb, pos: (usize, usize)) -> AlignRow {
        AlignRow {
            seq: SeqId::new(seq),
            loglik,
            pos,
        }
    }

    pub(crate) fn sample_reads() -> Vec<ReadRows> {
        vec![
            ReadRows {
                read_id: "read_1/1".into(),
                rows: vec![row(2, -690.0, (10, 20)), row(7, -10.25, (0, 0))],
            },
            ReadRows {
                read_id: "no rows".into(),
                rows: vec![],
            },
            ReadRows {
                read_id: "it's".into(),
                rows: vec![row(9, -0.5, (7, 8))],
            },
        ]
    }

    fn write_sample<W: std::io::Write>(out: W) -> W {
        let head = sample_head();
        let reads = sample_reads();
        let mut w = ReportWriter::begin(out, &head, reads.len()).unwrap();
        for r in &reads {
            w.push_read(&r.read_id, &r.rows).unwrap();
        }
        w.end_alignments().unwrap();
        w.finish(&Profile {
            write_s: 0.25,
            total_s: 2.5,
            peak_rss_kb: Some(123456),
        })
        .unwrap()
    }

    fn read_all<R: BufRead>(r: R) -> (ReportHead, usize, Vec<ReadRows>) {
        let (head, n, mut rr) = ReportReader::parse_head(r).unwrap();
        let mut got = Vec::new();
        let mut buf = ReadRows::default();
        while rr.next_read(&mut buf).unwrap() {
            got.push(buf.clone());
        }
        (head, n, got)
    }

    #[test]
    fn head_and_reads_round_trip_exactly() {
        let bytes = write_sample(Vec::new());
        let text = String::from_utf8(bytes.clone()).unwrap();
        assert!(text.starts_with("#NEXUS\n[PREMISE report v3]\n"));
        assert!(text.contains("\tINDEX 'idx with space.fmidx';\n\tINDEX_SHA256 aaaa"));
        assert!(text.contains("\t'read_1/1' 2 2\n\t2 5\n\t-6.9e2 -1.025e1\n\t0 1\n\t10 0\n\t10 0\n\t'no rows' 0\n\t'it''s' 1 1\n\t9\n\t-5e-1\n\t7\n\t1\n\t;\n"));
        assert!(text.contains("BEGIN PROFILE;\n\tRUNTIME write_s=2.5e-1 total_s=2.5e0;\n\tMEMORY peak_rss_kb=123456;\nEND;\n"));
        let (head, n, reads) = read_all(bytes.as_slice());
        assert_eq!(head, sample_head());
        assert_eq!(n, 3);
        assert_eq!(reads, sample_reads());
    }

    #[test]
    fn source_round_trips_when_set() {
        let mut head = sample_head();
        head.meta.source = Some("base run.report".into());
        let mut w = ReportWriter::begin(Vec::new(), &head, 0).unwrap();
        w.end_alignments().unwrap();
        let bytes = w.finish(&Profile::default()).unwrap();
        let text = String::from_utf8(bytes.clone()).unwrap();
        assert!(text.contains("\tNODES -;\n\tSOURCE 'base run.report';\n\tPARAMS "));
        let (got, n, _) = ReportReader::parse_head(bytes.as_slice()).unwrap();
        assert_eq!(got, head);
        assert_eq!(n, 0);
        assert_eq!(sample_head().meta.source, None);
    }

    #[test]
    fn zstd_report_round_trips() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir
            .path()
            .join("r.report.zst")
            .to_str()
            .unwrap()
            .to_string();
        let sink = write_sample(Sink::create(&path, true).unwrap());
        sink.finish().unwrap();
        let reader = crate::report::nexus::open_report(&path).unwrap();
        let (head, _, reads) = read_all(reader);
        assert_eq!(head, sample_head());
        assert_eq!(reads, sample_reads());
        let mut plain = String::new();
        crate::report::nexus::open_report(&path)
            .unwrap()
            .read_to_string(&mut plain)
            .unwrap();
        assert_eq!(plain.as_bytes(), write_sample(Vec::new()).as_slice());
    }

    #[test]
    fn unknown_blocks_and_commands_are_skipped() {
        let bytes = write_sample(Vec::new());
        let text = String::from_utf8(bytes).unwrap();
        let with_extra = text
            .replacen(
                "BEGIN NODES;",
                "BEGIN SETS;\n\tTAXSET all = 1-3;\nEND;\n\nBEGIN NODES;",
                1,
            )
            .replacen("\tTHREADS 8;", "\tTHREADS 8;\n\tFUTURE_KEY x=1;", 1);
        let (head, _, _) = read_all(with_extra.as_bytes());
        assert_eq!(head, sample_head());
        let err = ReportReader::parse_head("hello\n".as_bytes())
            .err()
            .unwrap()
            .to_string();
        assert!(err.contains("not a NEXUS report"), "{err}");
    }

    fn head_err(edit: impl Fn(&str) -> String) -> String {
        let text = String::from_utf8(write_sample(Vec::new())).unwrap();
        let edited = edit(&text);
        assert_ne!(edited, text, "the edit must change the report");
        match ReportReader::parse_head(edited.as_bytes()) {
            Err(e) => format!("{:#}", e),
            Ok(_) => panic!("parse_head accepted a broken report"),
        }
    }

    #[test]
    fn malformed_heads_are_rejected_with_a_reason() {
        // Other format versions, old or new.
        let e = head_err(|t| t.replacen("\tVERSION 3;", "\tVERSION 2;", 1));
        assert!(e.contains("not the version this PREMISE reads"), "{e}");
        let e = head_err(|t| t.replacen("\tVERSION 3;", "\tVERSION 4;", 1));
        assert!(e.contains("not the version"), "{e}");
        let e = head_err(|t| t.replacen("\tINDEX_SHA256 aaaa", "\tINDEX_SHA256 zzzz", 1));
        assert!(e.contains("INDEX_SHA256"), "{e}");
        // A missing block.
        let e = head_err(|t| {
            let a = t.find("BEGIN PROPS;").unwrap();
            let b = t[a..].find("END;\n").unwrap() + a + 5;
            format!("{}{}", &t[..a], &t[b..])
        });
        assert!(e.contains("no PROPS block"), "{e}");
        // TAXA / NODES disagreement.
        let e = head_err(|t| t.replacen("DIMENSIONS NTAX=10;", "DIMENSIONS NTAX=9;", 1));
        assert!(e.contains("TAXA: more than 9 labels"), "{e}");
        let e = head_err(|t| {
            t.replacen("DIMENSIONS NTAX=10;", "DIMENSIONS NTAX=9;", 1)
                .replacen("\t\ts3\n", "", 1)
        });
        assert!(e.contains("TAXA lists 9 labels but NODES has 10"), "{e}");
        let e = head_err(|t| {
            t.replacen(
                "BEGIN NODES;\n\tDIMENSIONS N=10;",
                "BEGIN NODES;\n\tDIMENSIONS N=9;",
                1,
            )
        });
        assert!(e.contains("NODES: more than 9 rows"), "{e}");
        let e = head_err(|t| t.replacen(" collapsed ", " melted ", 1));
        assert!(e.contains("unknown status 'melted'"), "{e}");
        let e = head_err(|t| {
            t.replacen(
                "\t3 1 species collapsed 2 2e1 4e-1 ",
                "\t3 1 species collapsed 2 2e1 NaN ",
                1,
            )
        });
        assert!(e.contains("NODES proportion"), "{e}");
        // Tree structure: parent before child, no self-parent, root first, unique ids,
        // final vertices only under final vertices.
        let e = head_err(|t| t.replacen("\t3 1 species collapsed", "\t3 3 species collapsed", 1));
        assert!(e.contains("its own parent"), "{e}");
        let e = head_err(|t| t.replacen("\t3 1 species collapsed", "\t3 99 species collapsed", 1));
        assert!(e.contains("before its parent 99"), "{e}");
        let e = head_err(|t| t.replacen("\t1 0 genus active", "\t1 - genus active", 1));
        assert!(e.contains("has no parent"), "{e}");
        let e = head_err(|t| t.replacen("\t0 - root active", "\t0 5 root active", 1));
        assert!(e.contains("must be the root"), "{e}");
        let e = head_err(|t| t.replacen("\t4 1 species active", "\t3 1 species active", 1));
        assert!(e.contains("duplicate vertex id 3"), "{e}");
        let e = head_err(|t| t.replacen("\t9 5 sequence pruned", "\t9 5 sequence active", 1));
        assert!(e.contains("active vertex 9 under pruned vertex 5"), "{e}");
        // Cross-references.
        let e = head_err(|t| t.replacen("\t2 8\n", "\t2 80\n", 1));
        assert!(
            e.contains("REFMAP: sequence 2 maps to unknown vertex 80"),
            "{e}"
        );
        // (EM_TRACE holds the same `3 4e-1 / 8 6e-1` rows; aim at the PROPS matrix.)
        let e = head_err(|t| {
            t.replacen(
                "[node proportion]\n\tMATRIX\n\t3 4e-1\n\t8 6e-1\n",
                "[node proportion]\n\tMATRIX\n\t3 4e-1\n\t88 6e-1\n",
                1,
            )
        });
        assert!(e.contains("PROPS: unknown vertex 88"), "{e}");
        // Truncated EM_TRACE.
        let e = head_err(|t| t.replacen("\t0 -1e1 4\n", "\t0 -1e1 40\n", 1));
        assert!(e.contains("EM_TRACE"), "{e}");
    }

    #[test]
    fn malformed_reads_are_rejected() {
        let text = String::from_utf8(write_sample(Vec::new())).unwrap();
        let read_err = |edited: String| -> String {
            assert_ne!(edited, text, "the edit must change the report");
            let (_, _, mut rr) = ReportReader::parse_head(edited.as_bytes()).unwrap();
            let mut buf = ReadRows::default();
            loop {
                match rr.next_read(&mut buf) {
                    Ok(true) => continue,
                    Ok(false) => panic!("accepted a broken ALIGNMENT_LIKELIHOODS block"),
                    Err(e) => return format!("{:#}", e),
                }
            }
        };
        let e = read_err(text.replacen("\t-6.9e2 -1.025e1\n", "\tNaN -1.025e1\n", 1));
        assert!(e.contains("not finite"), "{e}");
        let e = read_err(text.replacen("\t-6.9e2 -1.025e1\n", "\t-6.9e2 inf\n", 1));
        assert!(e.contains("not finite"), "{e}");
        let e = read_err(text.replacen("\t2 5\n", "\t2 0\n", 1));
        assert!(e.contains("difference 0"), "{e}");
        let e = read_err(text.replacen("\t2 5\n", "\t2 4294967295\n", 1));
        assert!(e.contains("overflows"), "{e}");
        let e = read_err(text.replacen("\t2 5\n", "\t2 -3\n", 1));
        assert!(e.contains("bad number"), "{e}");
        let e = read_err(text.replacen("\t0 1\n\t10 0\n", "\t0 1\n\t10\n", 1));
        assert!(
            e.contains("1 forward positions for read 'read_1/1', 2 announced"),
            "{e}"
        );
        let e = read_err(text.replacen("\t0 1\n\t10 0\n", "\t0 1\n\t10 0 7\n", 1));
        assert!(e.contains("more than 2 forward positions"), "{e}");
        // A hostile row count errors on the first array line instead of allocating.
        let e = read_err(text.replacen(
            "\t'read_1/1' 2 2\n",
            "\t'read_1/1' 1000000000000000000 2\n",
            1,
        ));
        assert!(
            e.contains("2 seqids for read 'read_1/1', 1000000000000000000 announced"),
            "{e}"
        );
        // Truncated file.
        let cut = text.find("\t'it''s' 1 1").unwrap();
        let e = read_err(text[..cut].to_string());
        assert!(e.contains("unexpected end of report"), "{e}");
        let cut = text.find("\t-5e-1\n").unwrap();
        let e = read_err(text[..cut].to_string());
        assert!(
            e.contains("truncated log-likelihoods of read 'it's'"),
            "{e}"
        );
        // The likelihood palette and its indices.
        let rd = "\t'read_1/1' 2 2\n";
        for (hdr, why) in [
            (
                "\t'read_1/1' 2 0\n",
                "0 distinct log-likelihoods for 2 rows",
            ),
            (
                "\t'read_1/1' 2 3\n",
                "3 distinct log-likelihoods for 2 rows",
            ),
            (
                "\t'read_1/1' 2\n",
                "expected 'read_id 0' or 'read_id n_rows n_lls'",
            ),
            (
                "\t'read_1/1' 0 1\n",
                "expected 'read_id 0' or 'read_id n_rows n_lls'",
            ),
        ] {
            let e = read_err(text.replacen(rd, hdr, 1));
            assert!(e.contains(why), "{hdr:?}: {e}");
        }
        let e = read_err(text.replacen("\t0 1\n\t10 0\n", "\t0 2\n\t10 0\n", 1));
        assert!(e.contains("log-likelihood index 2 is not below 2"), "{e}");
        let e = read_err(text.replacen("\t0 1\n\t10 0\n", "\t0\n\t10 0\n", 1));
        assert!(e.contains("1 log-likelihood indices"), "{e}");
        // Insert offsets: one or n values, and the reverse position must exist.
        let e = read_err(text.replacen("\t10 0\n\t10 0\n", "\t10 0\n\t10 0 1\n", 1));
        assert!(e.contains("more than 2 insert offsets"), "{e}");
        let e = read_err(text.replacen("\t10 0\n\t10 0\n", "\t10 0\n\t-11\n", 1));
        assert!(e.contains("out of range"), "{e}");
        let e = read_err(text.replacen("\t7\n\t1\n", "\t7\n\t-8\n", 1));
        assert!(e.contains("out of range"), "{e}");
    }

    #[test]
    fn adversarial_rows_round_trip_bit_for_bit() {
        let tiny = f64::from_bits(1); // smallest subnormal
        let x: f64 = -31.36650193982768;
        let next = f64::from_bits(x.to_bits() + 1); // one ulp away
        let big = i64::MAX as usize;
        let reads = vec![
            ReadRows {
                read_id: "zeros and ulps".into(),
                rows: vec![
                    row(0, -0.0, (0, 0)),
                    row(1, 0.0, (5, 2)),
                    row(2, x, (big, 3)),
                    row(3, next, (7, big)),
                    row(4, x, (9, 9)),
                    row(5, tiny, (1, 0)),
                    row(9, f64::MIN, (0, big)),
                ],
            },
            ReadRows {
                read_id: "one palette, shared offset".into(),
                rows: vec![
                    row(3, x, (100, 221)),
                    row(4, x, (5, 126)),
                    row(70, x, (0, 121)),
                ],
            },
            ReadRows {
                read_id: "negative shared offset".into(),
                rows: vec![row(1, -2.0, (50, 10)), row(2, -1.0, (40, 0))],
            },
            ReadRows {
                read_id: "single".into(),
                rows: vec![row(u32::MAX, f64::MAX.copysign(-1.0), (3, 1))],
            },
            ReadRows {
                read_id: "none".into(),
                rows: vec![],
            },
        ];
        let head = sample_head();
        let mut w = ReportWriter::begin(Vec::new(), &head, reads.len()).unwrap();
        for r in &reads {
            w.push_read(&r.read_id, &r.rows).unwrap();
        }
        w.end_alignments().unwrap();
        let bytes = w
            .finish(&Profile {
                write_s: 0.0,
                total_s: 0.0,
                peak_rss_kb: None,
            })
            .unwrap();
        let (_, n, got) = read_all(bytes.as_slice());
        assert_eq!(n, reads.len());
        assert_eq!(got.len(), reads.len());
        for (g, r) in got.iter().zip(&reads) {
            assert_eq!(g.read_id, r.read_id);
            assert_eq!(g.rows.len(), r.rows.len());
            for (a, b) in g.rows.iter().zip(&r.rows) {
                assert_eq!(a.seq, b.seq, "{}", r.read_id);
                assert_eq!(a.pos, b.pos, "{}", r.read_id);
                assert_eq!(a.loglik.to_bits(), b.loglik.to_bits(), "{}", r.read_id);
            }
        }
        let text = String::from_utf8(bytes).unwrap();
        // Shared values are written once.
        assert!(text.contains("\t'one palette, shared offset' 3 1\n\t3 1 66\n\t-3.136650193982768e1\n\t100 5 0\n\t121\n"), "{text}");
        assert!(
            text.contains(
                "\t'negative shared offset' 2 2\n\t1 1\n\t-2e0 -1e0\n\t0 1\n\t50 40\n\t-40\n"
            ),
            "{text}"
        );
    }
}
