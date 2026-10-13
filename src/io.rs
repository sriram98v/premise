//! File input and output

use crate::em::Posteriors;
use crate::group::Grouping;
use crate::tree::{
    conditional, count, is_synthetic_root, label, preorder, rank, Partition, SubtreeIndex, TaxTree,
    NO_NODE,
};
use crate::utils::{
    get_extension_from_filename, EMProb, Read, ReadAlignments, ReadID, ReadIdx, ReadPair,
};
use crate::RefIndex;
use anyhow::{Context, Result};
use bio::io::{fasta, fastq};
use chrono::Local;
use flate2::read::GzDecoder;
use haystackfm::alphabet::{self, encode_byte};
use haystackfm::occ::OccEncoding;
use haystackfm::{DnaSequence, FmIndexConfig as RefIndexConfig};
use phylo::node::NodeID;
use phylo::prelude::*;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::fs::File;
use std::io::{BufReader, Cursor, Write};

/// Default depth of the k-mer lookup tables that seed the SMEM search: fastest in a
/// depth sweep on the benchmark reference (deeper tables cost more cache misses than
/// the extensions they save).
pub const DEFAULT_LOOKUP_DEPTH: u32 = 8;

/// Largest accepted lookup depth: each index half then holds about 12 × 4^13 bytes (~0.8 GB).
pub const MAX_LOOKUP_DEPTH: u32 = 13;

/// Build a serialized FM-index from raw FASTA bytes, keeping every suffix-array entry and
/// the default k-mer lookup tables.
pub fn build_index_from_bytes(fasta_data: &[u8]) -> Result<(Vec<u8>, String)> {
    build_index_from_bytes_with(fasta_data, 1, DEFAULT_LOOKUP_DEPTH)
}

/// Build a serialized FM-index from raw FASTA bytes, keeping one suffix-array entry in
/// every `sa_sample_rate` (1 = the full suffix array). Larger rates shrink the index and
/// slow down locating occurrences. `lookup_depth` is the k of the k-mer tables that let
/// the SMEM search skip its first k extensions (0 = no tables); each index half grows by
/// about 12 × 4^k bytes.
pub fn build_index_from_bytes_with(
    fasta_data: &[u8],
    sa_sample_rate: u32,
    lookup_depth: u32,
) -> Result<(Vec<u8>, String)> {
    if sa_sample_rate == 0 {
        anyhow::bail!("sa_sample_rate must be at least 1");
    }
    if lookup_depth > MAX_LOOKUP_DEPTH {
        anyhow::bail!("lookup_depth must be at most {MAX_LOOKUP_DEPTH}");
    }
    let cursor = Cursor::new(fasta_data);
    let reader = BufReader::new(cursor);
    let records = fasta::Reader::new(reader).records();

    let mut headers: Vec<String> = vec![];
    let mut refs_texts: Vec<Vec<u8>> = vec![];

    for result in records {
        let record = result.map_err(|e| anyhow::anyhow!("FASTA parse error: {}", e))?;
        headers.push(record.id().to_string());
        refs_texts.push(record.seq().to_vec());
    }

    if refs_texts.is_empty() {
        return Err(anyhow::anyhow!("No sequences found in FASTA file"));
    }

    let mut seen: HashMap<&str, usize> = HashMap::new();
    for (idx, header) in headers.iter().enumerate() {
        if let Some(first) = seen.insert(header.as_str(), idx) {
            return Err(anyhow::anyhow!(
                "duplicate FASTA header '{}' (records {} and {}); \
                 reference headers must be unique",
                header,
                first + 1,
                idx + 1
            ));
        }
    }

    let log_str = format!(
        "Timestamp: {}\nNum references: {}",
        Local::now().format("%Y-%m-%d %H:%M:%S"),
        headers.len(),
    );
    println!("{}", log_str);

    let sequences: Vec<DnaSequence> = (0..refs_texts.len())
        .map(|i| -> anyhow::Result<DnaSequence> {
            let header = headers[i].as_str();
            let seq_str: String = str::from_utf8(&refs_texts[i])
                .map_err(|e| {
                    anyhow::anyhow!("FASTA sequence for record {} is not valid UTF-8: {}", i, e)
                })?
                .chars()
                .map(|c| {
                    let upper = c.to_ascii_uppercase();
                    match upper.is_ascii().then(|| encode_byte(upper as u8)).flatten() {
                        Some(code) if code != alphabet::SENTINEL => upper,
                        _ => 'N',
                    }
                })
                .collect();
            DnaSequence::from_str_with_header(&seq_str, header).map_err(|e| {
                anyhow::anyhow!(
                    "could not build DNA sequence for record '{}': {}",
                    header,
                    e
                )
            })
        })
        .collect::<anyhow::Result<Vec<DnaSequence>>>()?;

    let config = RefIndexConfig {
        sa_sample_rate,
        use_gpu: false,
        occ_encoding: OccEncoding::OneHot,
        build_lcp: false,
        lookup_depth,
        ..Default::default()
    };
    let fmidx = RefIndex::build_cpu(&sequences, &config)
        .map_err(|e| anyhow::anyhow!("FM-index construction error: {}", e))?;

    let bytes = fmidx
        .to_bytes()
        .map_err(|e| anyhow::anyhow!("FmIndex serialize error: {}", e))?;
    Ok((bytes, log_str))
}

/// Load a serialized [`RefIndex`] from `.fmidx` bytes.
pub fn load_index(bytes: &[u8], path: &str) -> Result<RefIndex> {
    RefIndex::from_bytes(bytes).map_err(|e| {
        anyhow::anyhow!(
            "failed to parse FM-index '{}': {}. Indexes built by a much older or a newer \
             version of premise are not readable; rebuild it with this version's `build` \
             subcommand.",
            path,
            e
        )
    })
}

/// Read and load the `.fmidx` file at `path`.
pub fn read_index_file(path: &str) -> Result<RefIndex> {
    let file_bytes = std::fs::read(path).with_context(|| {
        format!(
            "failed to read FM-index '{}' (build it first with the `build` subcommand)",
            path
        )
    })?;
    load_index(&file_bytes, path)
}

/// Collect FASTQ records from an iterator into compact [`Read`]s
fn collect_fastq_reads<I, E>(
    records: I,
    mate: &str,
    path: &str,
    key: impl Fn(&str) -> ReadID,
) -> HashMap<ReadID, Read>
where
    I: Iterator<Item = std::result::Result<fastq::Record, E>>,
{
    let mut dropped = 0usize;
    let mut out = HashMap::new();
    for rec in records {
        match rec {
            Ok(rec) => {
                out.insert(key(rec.id()), Read::from(rec));
            }
            Err(_) => dropped += 1,
        }
    }
    if dropped > 0 {
        eprintln!(
            "Warning: skipped {} malformed record(s) while reading {} file '{}'",
            dropped, mate, path
        );
    }
    out
}

/// Load one mate's reads from a FASTQ or gzipped FASTQ file into a read-ID map.
pub fn load_fastq(path: &str, mate: &str, suffix: &str) -> Result<HashMap<ReadID, Read>> {
    let f =
        File::open(path).with_context(|| format!("failed to open {mate} reads file '{path}'"))?;
    match get_extension_from_filename(path) {
        Some("gz") => {
            let records = fastq::Reader::from_bufread(BufReader::new(GzDecoder::new(f))).records();
            Ok(collect_fastq_reads(records, mate, path, |id| {
                ReadID(id.to_string())
            }))
        }
        Some("fastq") | Some("fq") => {
            let records = fastq::Reader::from_bufread(BufReader::new(f)).records();
            Ok(collect_fastq_reads(records, mate, path, |id| {
                ReadID(id.strip_suffix(suffix).unwrap_or(id).to_string())
            }))
        }
        _ => Err(anyhow::anyhow!("Unsupported {mate} file type: {path}")),
    }
}

/// Load forward (R1) reads; see [`load_fastq`].
pub fn load_fastq_forward(path: &str) -> Result<HashMap<ReadID, Read>> {
    load_fastq(path, "R1", "/1")
}

/// Load reverse (R2) reads; see [`load_fastq`].
pub fn load_fastq_reverse(path: &str) -> Result<HashMap<ReadID, Read>> {
    load_fastq(path, "R2", "/2")
}

/// Join R1 and R2 into pairs by read id
pub fn pair_reads(
    mut forward: HashMap<ReadID, Read>,
    reverse: HashMap<ReadID, Read>,
    r1_file: &str,
    r2_file: &str,
) -> Result<(Vec<ReadPair>, Vec<ReadID>)> {
    let mut pairs = Vec::with_capacity(reverse.len().min(forward.len()));
    for (read_id, r2) in reverse {
        if let Some(r1) = forward.remove(&read_id) {
            pairs.push(ReadPair { read_id, r1, r2 });
        }
    }
    if pairs.is_empty() {
        return Err(anyhow::anyhow!(
            "no valid read pairs after matching R1 '{}' with R2 '{}'; \
             R1 and R2 read IDs do not correspond (check for mismatched files or /1,/2 suffix handling)",
            r1_file,
            r2_file
        ));
    }
    let unpaired: Vec<ReadID> = forward.into_keys().collect();
    Ok((pairs, unpaired))
}

/// All R1 read ids in output order
pub fn all_read_ids<'a>(pairs: &'a [ReadPair], unpaired: &'a [ReadID]) -> BTreeSet<&'a ReadID> {
    pairs
        .iter()
        .map(|p| &p.read_id)
        .chain(unpaired.iter())
        .collect()
}

/// The output tables of [`crate::run_query`], each streamed to its own writer:
/// - `matches`: per-read assignment table (TSV)
/// - `posteriors`: per-read posterior over the tips (TSV)
/// - `props`: estimated abundance proportions, one row per reported tip (TSV)
/// - `nodes`: every vertex of the fitted taxonomy (TSV)
/// - `aligns`: raw per-read alignment likelihoods (TSV, same format as `run_alignment`)
pub struct QueryWriters<W: Write> {
    pub matches: W,
    pub posteriors: W,
    pub props: W,
    pub nodes: W,
    pub aligns: W,
}

/// `align` table
pub fn write_alignments<W: Write>(
    out: &mut W,
    fmidx: &RefIndex,
    read_ids: &BTreeSet<&ReadID>,
    aligns: &ReadAlignments<'_>,
) -> Result<()> {
    writeln!(
        out,
        "ReadID\tRefID\tProbability\tForward Positions\tReverse Position"
    )?;
    for read_id in read_ids {
        let Some(read_idx) = aligns.read_idx(read_id) else {
            writeln!(out, "{}\tunclassified\t-\t-\t-", **read_id)?;
            continue;
        };
        let (refs, row) = aligns.row(read_idx);
        for (ref_idx, likelihood) in refs.iter().zip(row) {
            let ref_id = fmidx.seq_header(*ref_idx).unwrap_or("");
            writeln!(
                out,
                "{}\t{}\t{:.5e}\t{}\t{}",
                **read_id,
                ref_id,
                likelihood.get_full_match_ll().exp(),
                likelihood.get_pos().0,
                likelihood.get_pos().1,
            )?;
        }
    }
    out.flush()?;
    Ok(())
}

/// Where each read ended up: the vertex and its subtree posterior.
pub type Assignments = HashMap<ReadIdx, (NodeID, EMProb)>;

/// `.posteriors` table: one row per (read, leaf): a sequence, or a collapsed taxon over
/// its sequences.
#[allow(clippy::too_many_arguments)]
pub fn write_posteriors<W: Write>(
    out: &mut W,
    grouping: &Grouping,
    tree: &TaxTree,
    part: &Partition,
    read_ids: &BTreeSet<&ReadID>,
    aligns: &ReadAlignments<'_>,
    assignments: &Assignments,
    posteriors: &Posteriors<'_>,
) -> Result<()> {
    writeln!(out, "ReadID\t{}\tPosterior", grouping.label_header())?;
    for read_id in read_ids.iter() {
        let assigned = aligns
            .read_idx(read_id)
            .filter(|read_idx| assignments.contains_key(read_idx));
        let Some(read_idx) = assigned else {
            writeln!(out, "{}\tunclassified\t-", **read_id)?;
            continue;
        };
        let mut rows: Vec<(&str, EMProb)> = posteriors
            .per_leaf(read_idx, part)
            .into_iter()
            .map(|(leaf, p)| (label(tree, leaf), p))
            .collect();
        rows.sort_by(|a, b| a.0.cmp(b.0));
        for (name, p) in rows {
            writeln!(out, "{}\t{}\t{:.5e}", **read_id, name, p)?;
        }
    }
    out.flush()?;
    Ok(())
}

/// `.matches` table: the vertex each read was assigned to (a sequence, a taxid, or `root`),
/// its posterior there, and the position of the read's best-scoring sequence under that
/// vertex. Reads that stopped at the synthetic root are unclassified. Returns the
/// classified reads.
#[allow(clippy::too_many_arguments)]
pub fn write_matches<W: Write>(
    out: &mut W,
    grouping: &Grouping,
    tree: &TaxTree,
    part: &Partition,
    read_ids: &BTreeSet<&ReadID>,
    aligns: &ReadAlignments<'_>,
    assignments: &Assignments,
    posteriors: &Posteriors<'_>,
) -> Result<HashSet<ReadIdx>> {
    let mut classified_reads: HashSet<ReadIdx> = HashSet::new();
    let subtree = SubtreeIndex::new(tree);
    writeln!(
        out,
        "ReadID\t{}\tPosterior\tForward Position\tReverse Position",
        grouping.label_header()
    )?;
    for read_id in read_ids.iter() {
        let assignment = aligns
            .read_idx(read_id)
            .and_then(|read_idx| assignments.get(&read_idx).map(|a| (read_idx, *a)));
        let Some((read_idx, (vertex, h))) =
            assignment.filter(|(_, (v, _))| !is_synthetic_root(tree, *v))
        else {
            writeln!(out, "{}\tunclassified\t-\t-\t-", **read_id)?;
            continue;
        };
        let (refs, als) = aligns.row(read_idx);
        let best = refs
            .iter()
            .zip(als)
            .filter(|(r, _)| {
                let leaf = posteriors.leaf_of(**r, part);
                leaf != NO_NODE && subtree.contains(vertex, leaf)
            })
            .max_by(|(_, a), (_, b)| {
                EMProb::total_cmp(&a.get_full_match_ll().exp(), &b.get_full_match_ll().exp())
            });
        let Some((_, alignment)) = best else {
            writeln!(out, "{}\tunclassified\t-\t-\t-", **read_id)?;
            continue;
        };
        writeln!(
            out,
            "{}\t{}\t{:.5e}\t{}\t{}",
            **read_id,
            label(tree, vertex),
            h,
            alignment.get_pos().0,
            alignment.get_pos().1,
        )?;
        classified_reads.insert(read_idx);
    }
    out.flush()?;
    Ok(classified_reads)
}

/// `.props` table: one row per leaf with its proportion `π_ℓ`: a sequence header, or the
/// taxid of a collapsed vertex. The leaves partition the mass, so the rows sum to one.
pub fn write_props<W: Write>(out: &mut W, tree: &TaxTree, node_pi: &[EMProb]) -> Result<()> {
    let mut rows: Vec<(&str, EMProb)> = tree
        .get_leaf_ids()
        .filter_map(|leaf| {
            let p = node_pi.get(leaf).copied().unwrap_or(0.0);
            (p > 0.0).then(|| (label(tree, leaf), p))
        })
        .collect();
    rows.sort_by(|a, b| a.0.cmp(b.0));
    for (name, p) in rows {
        writeln!(out, "{}\t{:.10e}", name, p)?;
    }
    out.flush()?;
    Ok(())
}

/// `.nodes` table: every vertex of the fitted taxonomy in pre-order. `Proportion` is the
/// vertex's marginal `π_v` (nested: a parent's includes its children's), `Conditional` its
/// share of its parent's mass, `Count` its expected reads, `Contigs` the sequences a leaf
/// holds (1 for a sequence, the pooled set for a collapsed vertex, 0 for an internal
/// vertex), `Rank` the taxonomic rank (`sequence` for an index record).
pub fn write_nodes<W: Write>(
    out: &mut W,
    tree: &TaxTree,
    part: &Partition,
    node_pi: &[EMProb],
) -> Result<()> {
    writeln!(
        out,
        "TaxID\tParent\tRank\tContigs\tCount\tProportion\tConditional"
    )?;
    let root = tree.get_root_id();
    for v in preorder(tree) {
        let contigs = part.size_of(v);
        let parent = if v == root {
            "-"
        } else {
            tree.get_node_parent_id(v).map_or("-", |p| label(tree, p))
        };
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{:.3}\t{:.10e}\t{:.10e}",
            label(tree, v),
            parent,
            rank(tree, v),
            contigs,
            count(tree, v),
            node_pi.get(v).copied().unwrap_or(0.0),
            conditional(tree, v),
        )?;
    }
    out.flush()?;
    Ok(())
}
/// The `.aligns` table
pub fn write_aligns<W: Write>(
    out: &mut W,
    fmidx: &RefIndex,
    read_ids: &BTreeSet<&ReadID>,
    aligns: &ReadAlignments<'_>,
) -> Result<()> {
    writeln!(out, "ReadID\tRefID\tProbability")?;
    for read_id in read_ids.iter() {
        let Some(read_idx) = aligns.read_idx(read_id) else {
            writeln!(out, "{}\tunclassified\t-", **read_id)?;
            continue;
        };
        let (refs, row) = aligns.row(read_idx);
        for (ref_idx, likelihood) in refs.iter().zip(row) {
            let ref_id = fmidx.seq_header(*ref_idx).unwrap_or("");
            writeln!(
                out,
                "{}\t{}\t{:.5e}",
                **read_id,
                ref_id,
                likelihood.get_full_match_ll().exp(),
            )?;
        }
    }
    out.flush()?;
    Ok(())
}
