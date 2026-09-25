//! File input and output

use crate::em::Posteriors;
use crate::utils::{
    get_extension_from_filename, EMProb, Read, ReadAlignments, ReadID, ReadIdx, ReadPair,
};
use crate::{RefIndex, SeqId};
use anyhow::{Context, Result};
use bio::io::{fasta, fastq};
use chrono::Local;
use flate2::read::GzDecoder;
use haystackfm::alphabet::{self, encode_byte};
use haystackfm::occ::OccEncoding;
use haystackfm::{DnaSequence, FmIndexConfig as RefIndexConfig};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs::File;
use std::io::{BufReader, Cursor, Write};

/// Build a serialized FM-index from raw FASTA bytes, keeping every suffix-array entry.
pub fn build_index_from_bytes(fasta_data: &[u8]) -> Result<(Vec<u8>, String)> {
    build_index_from_bytes_with(fasta_data, 1)
}

/// Build a serialized FM-index from raw FASTA bytes, keeping one suffix-array entry in
/// every `sa_sample_rate` (1 = the full suffix array). Larger rates shrink the index and
/// slow down locating occurrences.
pub fn build_index_from_bytes_with(
    fasta_data: &[u8],
    sa_sample_rate: u32,
) -> Result<(Vec<u8>, String)> {
    if sa_sample_rate == 0 {
        anyhow::bail!("sa_sample_rate must be at least 1");
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

/// The four output tables of [`crate::run_query`], each streamed to its own writer:
/// - `matches`: per-read MAP assignment table (TSV)
/// - `posteriors`: full posterior probability matrix (TSV)
/// - `props`: estimated reference abundance proportions (TSV)
/// - `aligns`: raw per-read alignment likelihoods (TSV, same format as `run_alignment`)
pub struct QueryWriters<W: Write> {
    pub matches: W,
    pub posteriors: W,
    pub props: W,
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

/// `.posteriors` table
pub fn write_posteriors<W: Write>(
    out: &mut W,
    fmidx: &RefIndex,
    read_ids: &BTreeSet<&ReadID>,
    aligns: &ReadAlignments<'_>,
    read_assignments: &HashMap<ReadIdx, SeqId>,
    posteriors: &Posteriors<'_>,
) -> Result<()> {
    writeln!(out, "ReadID\tRefID\tPosterior")?;
    for read_id in read_ids.iter() {
        let classified = aligns
            .read_idx(read_id)
            .filter(|read_idx| read_assignments.contains_key(read_idx));
        let Some(read_idx) = classified else {
            writeln!(out, "{}\tunclassified\t-", **read_id)?;
            continue;
        };
        let (refs, _) = aligns.row(read_idx);
        for ref_idx in refs {
            let ref_id = fmidx.seq_header(*ref_idx).unwrap_or("");
            writeln!(
                out,
                "{}\t{}\t{:.5e}",
                **read_id,
                ref_id,
                posteriors.get(read_idx, *ref_idx),
            )?;
        }
    }
    out.flush()?;
    Ok(())
}

/// `.matches` table
pub fn write_matches<W: Write>(
    out: &mut W,
    fmidx: &RefIndex,
    read_ids: &BTreeSet<&ReadID>,
    aligns: &ReadAlignments<'_>,
    read_assignments: &HashMap<ReadIdx, SeqId>,
    posteriors: &Posteriors<'_>,
    props_refs: &HashSet<SeqId>,
) -> Result<HashSet<ReadIdx>> {
    let mut classified_reads: HashSet<ReadIdx> = HashSet::new();
    writeln!(
        out,
        "ReadID\tRefID\tPosterior\tForward Position\tReverse Position"
    )?;
    for read_id in read_ids.iter() {
        let assignment = aligns
            .read_idx(read_id)
            .and_then(|read_idx| read_assignments.get(&read_idx).map(|r| (read_idx, *r)));
        let Some((read_idx, ref_idx)) = assignment.filter(|(_, r)| props_refs.contains(r)) else {
            writeln!(out, "{}\tunclassified\t-\t-\t-", **read_id)?;
            continue;
        };
        let Some(alignment) = aligns.get(read_idx, ref_idx) else {
            continue;
        };
        let ref_id = fmidx.seq_header(ref_idx).unwrap_or("");
        writeln!(
            out,
            "{}\t{}\t{:.5e}\t{}\t{}",
            **read_id,
            ref_id,
            posteriors.get(read_idx, ref_idx),
            alignment.get_pos().0,
            alignment.get_pos().1,
        )?;
        classified_reads.insert(read_idx);
    }
    out.flush()?;
    Ok(classified_reads)
}

/// The `.props` table: one line per reference with its proportion, in reference order.
pub fn write_props<W: Write>(
    out: &mut W,
    fmidx: &RefIndex,
    props: &HashMap<SeqId, EMProb>,
) -> Result<()> {
    // BTreeMap only to keep the output order deterministic.
    for (ref_idx, prop) in props.iter().collect::<BTreeMap<_, _>>() {
        let ref_id = fmidx.seq_header(*ref_idx).unwrap_or("");
        writeln!(out, "{}\t{:.10e}", ref_id, prop)?;
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
