//! Serialise a report. The head blocks are written at once; ALIGNMENT_LIKELIHOODS is
//! streamed one read at a time so nothing per read is held in memory.

use super::nexus::quote;
use super::{AlignRow, Profile, ReportHead, ReportTree};
use crate::utils::EMProb;
use crate::SeqId;
use anyhow::{anyhow, bail, Result};
use phylo::node::NodeID;
use std::collections::HashMap;
use std::io::Write;

fn opt_f(v: Option<EMProb>) -> String {
    v.map_or_else(|| "-".to_string(), |x| format!("{:e}", x))
}

fn opt_s(v: Option<&str>) -> String {
    v.map_or_else(|| "-".to_string(), |s| quote(s).into_owned())
}

/// Write `items` separated by single spaces.
fn write_joined<W: Write, T>(
    out: &mut W,
    items: impl Iterator<Item = T>,
    mut f: impl FnMut(&mut W, T) -> std::io::Result<()>,
) -> std::io::Result<()> {
    for (i, v) in items.enumerate() {
        if i > 0 {
            out.write_all(b" ")?;
        }
        f(out, v)?;
    }
    Ok(())
}

/// Write `s` NEXUS-quoted without building a temporary string when it is bare.
fn write_quoted<W: Write>(out: &mut W, s: &str) -> std::io::Result<()> {
    out.write_all(quote(s).as_bytes())
}

/// The tree as Newick over 1-based TAXA indices (vertex `i` of TAXA / NODES is `i + 1`),
/// the standard NEXUS form without a TRANSLATE table. Rank, status and the other
/// per-vertex fields live in NODES.
pub fn newick(tree: &ReportTree) -> String {
    let row: HashMap<NodeID, usize> = tree
        .vertices
        .iter()
        .enumerate()
        .map(|(i, v)| (v.id, i))
        .collect();
    let mut kids: Vec<Vec<usize>> = vec![Vec::new(); tree.vertices.len()];
    for (i, v) in tree.vertices.iter().enumerate() {
        if let Some(p) = v.parent.and_then(|p| row.get(&p)) {
            kids[*p].push(i);
        }
    }
    let mut out = String::new();
    fn rec(kids: &[Vec<usize>], i: usize, out: &mut String) {
        if !kids[i].is_empty() {
            out.push('(');
            for (k, &c) in kids[i].iter().enumerate() {
                if k > 0 {
                    out.push(',');
                }
                rec(kids, c, out);
            }
            out.push(')');
        }
        use std::fmt::Write as _;
        let _ = write!(out, "{}", i + 1);
    }
    if !tree.vertices.is_empty() {
        rec(&kids, 0, &mut out);
    }
    out
}

fn write_premise<W: Write>(out: &mut W, head: &ReportHead) -> Result<()> {
    let m = &head.meta;
    writeln!(out, "BEGIN PREMISE;")?;
    writeln!(out, "\tVERSION {};", m.version)?;
    writeln!(out, "\tTOOL {};", quote(&m.tool))?;
    writeln!(out, "\tTIMESTAMP {};", quote(&m.timestamp))?;
    writeln!(out, "\tLABEL_HEADER {};", quote(&m.label_header))?;
    writeln!(out, "\tINDEX {};", quote(&m.index))?;
    writeln!(out, "\tINDEX_SHA256 {};", opt_s(m.index_sha256.as_deref()))?;
    writeln!(out, "\tR1 {};", quote(&m.r1))?;
    writeln!(out, "\tR2 {};", quote(&m.r2))?;
    writeln!(out, "\tSEQ2TAXID {};", opt_s(m.seq2taxid.as_deref()))?;
    writeln!(out, "\tNODES {};", opt_s(m.nodes.as_deref()))?;
    if let Some(src) = &m.source {
        writeln!(out, "\tSOURCE {};", quote(src))?;
    }
    let p = &m.params;
    writeln!(
        out,
        "\tPARAMS mem={} eps_1={:e} eps_2={:e} iter={} rho={:e} omega={:e} em_threshold={:e} em_stop={} init={} penalty={};",
        p.mem,
        p.eps_1,
        p.eps_2,
        p.iter,
        p.rho,
        p.omega,
        p.em_threshold,
        quote(&p.em_stop),
        quote(&p.init),
        u8::from(p.penalty)
    )?;
    writeln!(out, "\tTHREADS {};", m.threads)?;
    writeln!(
        out,
        "\tREADS total={} aligned={};",
        m.reads_total, m.reads_aligned
    )?;
    let t = &m.taxonomy;
    writeln!(
        out,
        "\tTAXONOMY ref_vertices={} ref_tips={} unknown_taxids={} unmapped_seqs={} contracted_vertices={} contracted_tips={} final_tips={};",
        t.ref_vertices,
        t.ref_tips,
        t.unknown_taxids,
        t.unmapped_seqs,
        t.contracted_vertices,
        t.contracted_tips,
        t.final_tips
    )?;
    let e = &m.em;
    writeln!(
        out,
        "\tEM iterations={} pruned={} collapsed={} root_collapsed={} final_ll={};",
        e.iterations,
        e.pruned,
        e.collapsed,
        u8::from(e.root_collapsed),
        if e.final_ll.is_finite() {
            format!("{:e}", e.final_ll)
        } else {
            "-".to_string()
        }
    )?;
    let r = &m.runtime;
    writeln!(
        out,
        "\tRUNTIME load_s={:e} align_s={:e} em_s={:e};",
        r.load_s, r.align_s, r.em_s
    )?;
    writeln!(out, "END;")?;
    writeln!(out)?;
    Ok(())
}

fn write_em_trace<W: Write>(out: &mut W, head: &ReportHead) -> Result<()> {
    writeln!(out, "BEGIN EM_TRACE;")?;
    writeln!(out, "\tDIMENSIONS NITER={};", head.em_trace.len())?;
    writeln!(out, "\t[iter penalized_ll n_leaves / node proportion]")?;
    writeln!(out, "\tMATRIX")?;
    for (i, it) in head.em_trace.iter().enumerate() {
        writeln!(out, "\t{} {:e} {}", i, it.ll, it.leaves.len())?;
        for (node, p) in &it.leaves {
            writeln!(out, "\t{} {:e}", node, p)?;
        }
    }
    writeln!(out, "\t;")?;
    writeln!(out, "END;")?;
    writeln!(out)?;
    Ok(())
}

/// TAXA (labels), TREES (topology over TAXA indices) and NODES (every vertex field).
fn write_tree_blocks<W: Write>(out: &mut W, tree: &ReportTree) -> Result<()> {
    writeln!(out, "BEGIN TAXA;")?;
    writeln!(out, "\tDIMENSIONS NTAX={};", tree.vertices.len())?;
    writeln!(out, "\tTAXLABELS")?;
    for v in &tree.vertices {
        writeln!(out, "\t\t{}", quote(&v.label))?;
    }
    writeln!(out, "\t;")?;
    writeln!(out, "END;")?;
    writeln!(out)?;

    writeln!(out, "BEGIN TREES;")?;
    writeln!(out, "\tTREE taxonomy = [&R] {};", newick(tree))?;
    writeln!(out, "END;")?;
    writeln!(out)?;

    writeln!(out, "BEGIN NODES;")?;
    writeln!(out, "\tDIMENSIONS N={};", tree.vertices.len())?;
    writeln!(
        out,
        "\t[id parent rank status contigs count proportion conditional; the label of row i is TAXLABELS i]"
    )?;
    writeln!(out, "\tMATRIX")?;
    for v in &tree.vertices {
        writeln!(
            out,
            "\t{} {} {} {} {} {} {:e} {}",
            v.id,
            v.parent.map_or_else(|| "-".to_string(), |p| p.to_string()),
            quote(&v.rank),
            v.status.as_str(),
            v.contigs,
            opt_f(v.count),
            v.proportion,
            opt_f(v.conditional),
        )?;
    }
    writeln!(out, "\t;")?;
    writeln!(out, "END;")?;
    writeln!(out)?;
    Ok(())
}

/// REFMAP. The proportions of the final leaves are not written: they are the NODES
/// proportions of the final tree's leaves, and the reader takes them from there.
fn write_refmap<W: Write>(out: &mut W, head: &ReportHead) -> Result<()> {
    writeln!(out, "BEGIN REFMAP;")?;
    writeln!(out, "\tDIMENSIONS N={};", head.refmap.entries.len())?;
    writeln!(out, "\t[seqid node]")?;
    writeln!(out, "\tMATRIX")?;
    for (seq, node) in &head.refmap.entries {
        writeln!(out, "\t{} {}", seq.index(), node)?;
    }
    writeln!(out, "\t;")?;
    writeln!(out, "END;")?;
    writeln!(out)?;
    Ok(())
}

/// Write every block before ALIGNMENT_LIKELIHOODS.
pub fn write_head<W: Write>(out: &mut W, head: &ReportHead) -> Result<()> {
    writeln!(out, "#NEXUS")?;
    writeln!(out, "[PREMISE report v{}]", head.meta.version)?;
    writeln!(out)?;
    write_premise(out, head)?;
    write_em_trace(out, head)?;
    write_tree_blocks(out, &head.tree)?;
    write_refmap(out, head)
}

/// Write `seqs` as the first id then the difference to the previous one; the ids must
/// be strictly ascending.
fn write_seqids<W: Write>(
    out: &mut W,
    read_id: &str,
    seqs: impl Iterator<Item = SeqId>,
) -> Result<()> {
    let mut prev: u32 = 0;
    for (i, s) in seqs.enumerate() {
        let seq = u32::try_from(s.index())
            .map_err(|_| anyhow!("sequence id {} does not fit u32", s.index()))?;
        if i == 0 {
            write!(out, "{}", seq)?;
        } else {
            let delta = seq.checked_sub(prev).filter(|&d| d > 0).ok_or_else(|| {
                anyhow!(
                    "read '{}': rows are not in strictly ascending sequence order ({} after {})",
                    read_id,
                    seq,
                    prev
                )
            })?;
            write!(out, " {}", delta)?;
        }
        prev = seq;
    }
    Ok(())
}

/// Where a [`ReportWriter`] is: streaming reads, or past ALIGNMENT_LIKELIHOODS.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Phase {
    Alignments,
    Done,
}

/// Streams ALIGNMENT_LIKELIHOODS after [`write_head`].
pub struct ReportWriter<W: Write> {
    out: W,
    expected: usize,
    pushed: usize,
    phase: Phase,
    /// Per-read scratch: the distinct log-likelihoods (as bits, first-appearance order),
    /// each row's index into them, and each row's insert offset.
    palette: Vec<u64>,
    slot: HashMap<u64, usize>,
    index: Vec<usize>,
    offsets: Vec<i64>,
}

impl<W: Write> ReportWriter<W> {
    /// Write the head and open ALIGNMENT_LIKELIHOODS for `n_reads` reads.
    pub fn begin(mut out: W, head: &ReportHead, n_reads: usize) -> Result<Self> {
        write_head(&mut out, head)?;
        writeln!(out, "BEGIN ALIGNMENT_LIKELIHOODS;")?;
        writeln!(out, "\tDIMENSIONS NREADS={};", n_reads)?;
        writeln!(
            out,
            "\t[per read: read_id 0, or read_id n_rows n_lls then: seqids (first absolute, then differences); the n_lls distinct log-likelihoods in order of first appearance; each row's index into them (omitted when n_lls is 1); forward positions; insert offsets reverse - forward (a single value when every row has the same)]"
        )?;
        writeln!(out, "\tMATRIX")?;
        Ok(Self {
            out,
            expected: n_reads,
            pushed: 0,
            phase: Phase::Alignments,
            palette: Vec::new(),
            slot: HashMap::new(),
            index: Vec::new(),
            offsets: Vec::new(),
        })
    }

    fn expect_phase(&self, phase: Phase, what: &str) -> Result<()> {
        if self.phase != phase {
            bail!("report writer: {} out of order", what);
        }
        Ok(())
    }

    /// One read's alignments, in [`SeqId`] order.
    pub fn push_read(&mut self, read_id: &str, rows: &[AlignRow]) -> Result<()> {
        self.expect_phase(Phase::Alignments, "ALIGNMENT_LIKELIHOODS read")?;
        self.encode(read_id, rows)?;
        self.out.write_all(b"\t")?;
        write_quoted(&mut self.out, read_id)?;
        if rows.is_empty() {
            self.out.write_all(b" 0\n")?;
            self.pushed += 1;
            return Ok(());
        }
        writeln!(self.out, " {} {}", rows.len(), self.palette.len())?;
        self.out.write_all(b"\t")?;
        write_seqids(&mut self.out, read_id, rows.iter().map(|r| r.seq))?;
        self.out.write_all(b"\n\t")?;
        write_joined(&mut self.out, self.palette.iter(), |o, &b| {
            write!(o, "{:e}", EMProb::from_bits(b))
        })?;
        if self.palette.len() > 1 {
            self.out.write_all(b"\n\t")?;
            write_joined(&mut self.out, self.index.iter(), |o, v| write!(o, "{}", v))?;
        }
        self.out.write_all(b"\n\t")?;
        write_joined(&mut self.out, rows.iter().map(|r| r.pos.0), |o, v| {
            write!(o, "{}", v)
        })?;
        self.out.write_all(b"\n\t")?;
        let first = self.offsets[0];
        if self.offsets.iter().all(|&o| o == first) {
            write!(self.out, "{}", first)?;
        } else {
            write_joined(&mut self.out, self.offsets.iter(), |o, v| {
                write!(o, "{}", v)
            })?;
        }
        self.out.write_all(b"\n")?;
        self.pushed += 1;
        Ok(())
    }

    /// Fill the palette, index and offset scratch for `rows`.
    fn encode(&mut self, read_id: &str, rows: &[AlignRow]) -> Result<()> {
        self.palette.clear();
        self.slot.clear();
        self.index.clear();
        self.offsets.clear();
        for r in rows {
            if !r.loglik.is_finite() {
                bail!(
                    "read '{}': log-likelihood {} is not finite",
                    read_id,
                    r.loglik
                );
            }
            let bits = r.loglik.to_bits();
            let k = *self.slot.entry(bits).or_insert_with(|| {
                self.palette.push(bits);
                self.palette.len() - 1
            });
            self.index.push(k);
            let pos = |p: usize| {
                i64::try_from(p)
                    .map_err(|_| anyhow!("read '{}': position {} does not fit i64", read_id, p))
            };
            self.offsets.push(pos(r.pos.1)? - pos(r.pos.0)?);
        }
        Ok(())
    }

    /// Close ALIGNMENT_LIKELIHOODS, checking every announced read was pushed.
    pub fn end_alignments(&mut self) -> Result<()> {
        self.expect_phase(Phase::Alignments, "end of ALIGNMENT_LIKELIHOODS")?;
        if self.pushed != self.expected {
            bail!(
                "ALIGNMENT_LIKELIHOODS: {} reads pushed, {} announced",
                self.pushed,
                self.expected
            );
        }
        writeln!(self.out, "\t;")?;
        writeln!(self.out, "END;")?;
        writeln!(self.out)?;
        self.phase = Phase::Done;
        Ok(())
    }

    /// The PROFILE block, after [`Self::end_alignments`]; returns the sink so the caller
    /// can finish it.
    pub fn finish(mut self, profile: &Profile) -> Result<W> {
        self.expect_phase(Phase::Done, "PROFILE")?;
        writeln!(self.out, "BEGIN PROFILE;")?;
        writeln!(
            self.out,
            "\tRUNTIME write_s={:e} total_s={:e};",
            profile.write_s, profile.total_s
        )?;
        writeln!(
            self.out,
            "\tMEMORY peak_rss_kb={};",
            profile
                .peak_rss_kb
                .map_or_else(|| "-".to_string(), |k| k.to_string())
        )?;
        writeln!(self.out, "END;")?;
        self.out.flush()?;
        Ok(self.out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::report::fixtures::fitted;

    #[test]
    fn newick_is_over_one_based_taxa_indices() {
        let (stored, tips, fit_tree, part, fit) = fitted();
        let rt = ReportTree::from_fit(&stored, &tips, &fit_tree, &part, &fit);
        // Rewrite the expected tree over node ids into TAXA indices.
        let ix = |id: usize| {
            1 + rt
                .vertices
                .iter()
                .position(|v| v.id == id)
                .expect("vertex in the tree")
        };
        let nw = newick(&rt);
        assert!(!nw.contains('['), "{nw}");
        assert!(nw.ends_with(&format!("){}", ix(0))), "{nw}");
        assert!(
            nw.contains(&format!("({},{}){}", ix(6), ix(7), ix(3))),
            "{nw}"
        );
        assert!(
            nw.contains(&format!("(({}){}){}", ix(9), ix(5), ix(2))),
            "{nw}"
        );
        let mut head = Vec::new();
        write_tree_blocks(&mut head, &rt).unwrap();
        let text = String::from_utf8(head).unwrap();
        assert!(!text.contains("TRANSLATE"), "{text}");
        assert!(
            text.contains(&format!("\tTREE taxonomy = [&R] {};\n", nw)),
            "{text}"
        );
    }

    #[test]
    fn per_read_counts_are_checked() {
        use crate::report::read::tests::sample_head;
        let head = sample_head();
        let mut w = ReportWriter::begin(Vec::new(), &head, 2).unwrap();
        w.push_read("only one", &[]).unwrap();
        let err = w.end_alignments().unwrap_err().to_string();
        assert!(err.contains("1 reads pushed, 2 announced"), "{err}");
        let mut w = ReportWriter::begin(Vec::new(), &head, 0).unwrap();
        w.end_alignments().unwrap();
        let out = w
            .finish(&Profile {
                write_s: 0.0,
                total_s: 0.0,
                peak_rss_kb: None,
            })
            .unwrap();
        let text = String::from_utf8(out).unwrap();
        assert!(text.contains("peak_rss_kb=-;"));
        // Only what the EM cannot regenerate: no PROPS, no READ_POSTERIORS.
        assert!(
            !text.contains("BEGIN PROPS;") && !text.contains("READ_POSTERIORS"),
            "{text}"
        );
    }

    #[test]
    fn blocks_are_written_in_order() {
        use crate::report::read::tests::sample_head;
        let head = sample_head();
        let profile = Profile::default();
        let w = ReportWriter::begin(Vec::new(), &head, 0).unwrap();
        assert!(w.finish(&profile).is_err(), "PROFILE before the reads");
        let mut w = ReportWriter::begin(Vec::new(), &head, 1).unwrap();
        w.push_read("r", &[]).unwrap();
        w.end_alignments().unwrap();
        assert!(w.push_read("late", &[]).is_err());
        assert!(w.end_alignments().is_err());
    }

    #[test]
    fn unsorted_rows_are_an_error() {
        use crate::report::read::tests::sample_head;
        use crate::SeqId;
        let head = sample_head();
        let r = |seq: u32| AlignRow {
            seq: SeqId::new(seq),
            loglik: -1.0,
            pos: (0, 0),
        };
        let mut w = ReportWriter::begin(Vec::new(), &head, 1).unwrap();
        let err = w.push_read("x", &[r(5), r(5)]).unwrap_err().to_string();
        assert!(err.contains("strictly ascending"), "{err}");
        let mut w = ReportWriter::begin(Vec::new(), &head, 1).unwrap();
        let err = w.push_read("x", &[r(5), r(2)]).unwrap_err().to_string();
        assert!(err.contains("strictly ascending"), "{err}");
    }

    #[test]
    fn rows_are_written_with_a_likelihood_palette_and_insert_offsets() {
        use crate::report::read::tests::sample_head;
        use crate::SeqId;
        let head = sample_head();
        let mut w = ReportWriter::begin(Vec::new(), &head, 4).unwrap();
        let r = |seq: u32, loglik: EMProb, pos: (usize, usize)| AlignRow {
            seq: SeqId::new(seq),
            loglik,
            pos,
        };
        // Two distinct likelihoods, one insert offset.
        w.push_read(
            "a b",
            &[r(2, -1.5, (1, 2)), r(7, -300.0, (3, 4)), r(8, -1.5, (5, 6))],
        )
        .unwrap();
        // One likelihood (no index line), offsets differ (one negative).
        w.push_read("mixed", &[r(0, -2.0, (10, 4)), r(3, -2.0, (20, 25))])
            .unwrap();
        // -0.0 and 0.0 are different values.
        w.push_read("zeros", &[r(1, -0.0, (0, 0)), r(4, 0.0, (0, 0))])
            .unwrap();
        w.push_read("empty", &[]).unwrap();
        w.end_alignments().unwrap();
        let text = String::from_utf8(w.out).unwrap();
        let start = text.find("BEGIN ALIGNMENT_LIKELIHOODS;").unwrap();
        let block = &text[start + text[start..].find("\tMATRIX\n").unwrap() + 8..];
        assert!(
            block.starts_with(concat!(
                "\t'a b' 3 2\n\t2 5 1\n\t-1.5e0 -3e2\n\t0 1 0\n\t1 3 5\n\t1\n",
                "\tmixed 2 1\n\t0 3\n\t-2e0\n\t10 20\n\t-6 5\n",
                "\tzeros 2 2\n\t1 3\n\t-0e0 0e0\n\t0 1\n\t0 0\n\t0\n",
                "\tempty 0\n\t;\nEND;\n"
            )),
            "{block}"
        );
    }

    /// `cargo test --release -- --ignored --nocapture report_writer_timing`: per-read cost
    /// of ALIGNMENT_LIKELIHOODS for reads with many alignments.
    #[test]
    #[ignore]
    fn report_writer_timing() {
        use crate::report::read::tests::sample_head;
        use crate::SeqId;
        let head = sample_head();
        let (n_reads, per_read) = (20_000, 1000);
        let rows: Vec<AlignRow> = (0..per_read)
            .map(|i| AlignRow {
                seq: SeqId::new(3 * i as u32 + 1),
                loglik: -((i % 7) as EMProb) - 0.5,
                pos: (i, i + 300),
            })
            .collect();
        /// Counts bytes, so every value is really formatted (`io::sink` may skip it).
        struct Count(usize);
        impl Write for Count {
            fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
                self.0 += buf.len();
                Ok(buf.len())
            }
            fn flush(&mut self) -> std::io::Result<()> {
                Ok(())
            }
        }
        let mut w = ReportWriter::begin(Count(0), &head, n_reads).unwrap();
        let t = std::time::Instant::now();
        for _ in 0..n_reads {
            w.push_read("read", &rows).unwrap();
        }
        let aligns = t.elapsed().as_secs_f64();
        w.end_alignments().unwrap();
        eprintln!(
            "{n_reads} reads x {per_read} rows: ALIGNMENT_LIKELIHOODS {:.1} us/read",
            aligns * 1e6 / n_reads as f64
        );
    }

    #[test]
    fn positions_beyond_i64_are_an_error() {
        use crate::report::read::tests::sample_head;
        use crate::SeqId;
        let head = sample_head();
        let mut w = ReportWriter::begin(Vec::new(), &head, 1).unwrap();
        let row = AlignRow {
            seq: SeqId::new(0),
            loglik: -1.0,
            pos: (0, usize::MAX),
        };
        let err = w.push_read("x", &[row]).unwrap_err().to_string();
        assert!(err.contains("position"), "{err}");
    }

    #[test]
    fn non_finite_likelihoods_are_an_error() {
        use crate::report::read::tests::sample_head;
        use crate::SeqId;
        let head = sample_head();
        for ll in [EMProb::NAN, EMProb::INFINITY, EMProb::NEG_INFINITY] {
            let mut w = ReportWriter::begin(Vec::new(), &head, 1).unwrap();
            let row = AlignRow {
                seq: SeqId::new(0),
                loglik: ll,
                pos: (0, 0),
            };
            let err = w.push_read("x", &[row]).unwrap_err().to_string();
            assert!(err.contains("not finite"), "{err}");
        }
    }
}
