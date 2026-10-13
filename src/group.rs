//! Grouping of index sequences into the taxa the EM estimates
//!
//! The FM-index holds one sequence per FASTA record, and every record is a leaf of the
//! EM's tree ([`crate::tree`]). Sequences absent from the map hang under the root.

use crate::{RefIndex, SeqId};
use anyhow::{bail, Context, Result};
use std::collections::HashMap;

pub struct Grouping {
    /// `group_of[seq.index()]` = ordinal of the column that sequence belongs to.
    group_of: Vec<u32>,
    /// Reported id of each column, by ordinal: a taxid, or the header of an unmapped sequence.
    names: Vec<String>,
    /// Header of every index sequence, by [`SeqId`].
    seq_names: Vec<String>,
    /// The taxid of each column, by ordinal; `None` for an unmapped sequence.
    taxids: Vec<Option<u32>>,
    /// Sequences the map gave a taxid (the rest are their own column).
    n_mapped: usize,
    /// Built from a `--seq2taxid` map (tables then label the column `TaxID`).
    from_map: bool,
}

impl Grouping {
    /// Every sequence is its own column: the legacy per-record EM.
    pub fn identity(fmidx: &RefIndex) -> Self {
        let mut g = Self::from_pairs(fmidx, &HashMap::new());
        g.from_map = false;
        g
    }

    /// Group by a seqid2taxid map (`seqid<whitespace>taxid[<whitespace>…]`, see
    /// [`parse_map`]). Warns on stderr when index sequences have no taxid.
    pub fn from_map(fmidx: &RefIndex, path: &str) -> Result<Self> {
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("failed to read seq2taxid map '{}'", path))?;
        let map = parse_map(&text).with_context(|| format!("seq2taxid map '{}'", path))?;
        let g = Self::from_pairs(fmidx, &map);
        let unmapped = g.n_unmapped();
        if unmapped > 0 {
            eprintln!(
                "Warning: {} of {} index sequences have no taxid in '{}'; each is kept as its own column named by its header",
                unmapped,
                g.n_sequences(),
                path
            );
        }
        Ok(g)
    }

    /// Group by an in-memory `seq_id -> taxid` map. A sequence matches on its full header,
    /// else on the header's first whitespace-delimited token.
    pub fn from_pairs(fmidx: &RefIndex, map: &HashMap<String, u32>) -> Self {
        let n = fmidx.num_sequences() as usize;
        let mut group_of = Vec::with_capacity(n);
        let mut names: Vec<String> = Vec::new();
        let mut taxids: Vec<Option<u32>> = Vec::new();
        let mut seq_names: Vec<String> = Vec::with_capacity(n);
        let mut ordinal_of: HashMap<u32, u32> = HashMap::new();
        let mut n_mapped = 0;
        for i in 0..n {
            let header = fmidx.seq_header(SeqId::new(i as u32)).unwrap_or("");
            seq_names.push(header.to_string());
            let taxid = map
                .get(header)
                .or_else(|| map.get(header.split_whitespace().next().unwrap_or("")));
            let ordinal = match taxid {
                Some(&taxid) => {
                    n_mapped += 1;
                    *ordinal_of.entry(taxid).or_insert_with(|| {
                        names.push(taxid.to_string());
                        taxids.push(Some(taxid));
                        (names.len() - 1) as u32
                    })
                }
                None => {
                    names.push(header.to_string());
                    taxids.push(None);
                    (names.len() - 1) as u32
                }
            };
            group_of.push(ordinal);
        }
        Self {
            group_of,
            names,
            seq_names,
            taxids,
            n_mapped,
            from_map: true,
        }
    }

    pub fn n_groups(&self) -> usize {
        self.names.len()
    }

    pub fn n_sequences(&self) -> usize {
        self.group_of.len()
    }

    /// Index sequences the map gave no taxid.
    pub fn n_unmapped(&self) -> usize {
        self.n_sequences() - self.n_mapped
    }

    /// True when no two sequences share a column (the EM column of a sequence is itself).
    pub fn is_identity(&self) -> bool {
        self.names.len() == self.group_of.len()
    }

    /// The EM column of a sequence.
    pub fn group(&self, seq: SeqId) -> SeqId {
        SeqId::new(self.group_of[seq.index()])
    }

    /// The reported id of an EM column.
    pub fn name(&self, col: SeqId) -> &str {
        &self.names[col.index()]
    }

    /// The header of an index sequence (the label of its leaf).
    pub fn seq_name(&self, seq: SeqId) -> &str {
        &self.seq_names[seq.index()]
    }

    /// The taxid of an EM column; `None` for a sequence the map did not name.
    pub fn taxid_of(&self, col: SeqId) -> Option<u32> {
        self.taxids[col.index()]
    }

    /// Header of the id column in `.matches` (and `LABEL_HEADER` of the report).
    pub fn label_header(&self) -> &'static str {
        if self.from_map {
            "TaxID"
        } else {
            "RefID"
        }
    }

    /// Extra line for the query log; empty without a map so the log is unchanged.
    pub fn log_suffix(&self) -> String {
        if !self.from_map {
            return String::new();
        }
        format!(
            "\nTaxa: {} from {} sequences ({} with a taxid, {} without, kept under the root)",
            self.n_groups(),
            self.n_sequences(),
            self.n_mapped,
            self.n_unmapped(),
        )
    }
}

/// `seqid<whitespace>taxid[<whitespace>…]` lines into a map (kraken / centrifuge
/// `seqid2taxid.map` style). Blank lines and lines starting with `#` are skipped; extra
/// columns are ignored. A taxid must be a positive integer, and a seqid must not name two
/// taxids (an exact repeat is fine).
pub fn parse_map(text: &str) -> Result<HashMap<String, u32>> {
    let mut map: HashMap<String, u32> = HashMap::new();
    for (lineno, line) in text.lines().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let mut cols = trimmed.split_whitespace();
        let (Some(seq), Some(tax)) = (cols.next(), cols.next()) else {
            bail!(
                "line {}: expected seqid<whitespace>taxid, got {:?}",
                lineno + 1,
                line
            );
        };
        let taxid = match tax.parse::<u32>() {
            Ok(t) if t > 0 => t,
            _ => bail!(
                "line {}: taxid '{}' for '{}' is not a positive integer",
                lineno + 1,
                tax,
                seq
            ),
        };
        if let Some(&prev) = map.get(seq) {
            if prev != taxid {
                bail!(
                    "line {}: '{}' mapped to both taxid {} and {}",
                    lineno + 1,
                    seq,
                    prev,
                    taxid
                );
            }
            continue;
        }
        map.insert(seq.to_string(), taxid);
    }
    Ok(map)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::{build_index_from_bytes, load_index};

    const FASTA: &[u8] =
        b">g1|c1 first contig\nACGTACGTACGTACGTACGTACGTACGT\n>g1|c2\nTTAGCCTTAGCCTTAGCCTTAGCCTTAG\n>g2|c1\nGGGACCCTTTGGGACCCTTTGGGACCCT\n";

    fn index() -> RefIndex {
        let (bytes, _) = build_index_from_bytes(FASTA).expect("build");
        load_index(&bytes, "test.fmidx").expect("load")
    }

    fn pairs(items: &[(&str, u32)]) -> HashMap<String, u32> {
        items.iter().map(|(a, b)| (a.to_string(), *b)).collect()
    }

    #[test]
    fn identity_keeps_one_column_per_sequence_named_by_header() {
        let idx = index();
        let g = Grouping::identity(&idx);
        assert!(g.is_identity());
        assert_eq!(g.n_groups(), 3);
        assert_eq!(g.n_unmapped(), 3);
        assert_eq!(g.label_header(), "RefID");
        assert!(g.log_suffix().is_empty());
        for i in 0..3 {
            assert_eq!(g.group(SeqId::new(i)), SeqId::new(i));
            assert_eq!(
                g.name(SeqId::new(i)),
                idx.seq_header(SeqId::new(i)).unwrap()
            );
        }
    }

    #[test]
    fn map_merges_sequences_and_keeps_unmapped_as_own_column() {
        let idx = index();
        let g = Grouping::from_pairs(&idx, &pairs(&[("g1|c1", 7), ("g1|c2", 7)]));
        assert!(!g.is_identity());
        assert_eq!(g.n_groups(), 2);
        assert_eq!(g.n_unmapped(), 1);
        assert_eq!(g.group(SeqId::new(0)), g.group(SeqId::new(1)));
        assert_eq!(g.name(g.group(SeqId::new(0))), "7");
        assert_eq!(g.taxid_of(g.group(SeqId::new(0))), Some(7));
        assert_eq!(g.taxid_of(g.group(SeqId::new(2))), None);
        assert_ne!(g.group(SeqId::new(2)), g.group(SeqId::new(0)));
        assert_eq!(g.name(g.group(SeqId::new(2))), "g2|c1");
        assert_eq!(g.seq_name(SeqId::new(0)), "g1|c1");
        assert_eq!(g.seq_name(SeqId::new(2)), "g2|c1");
        assert_eq!(g.label_header(), "TaxID");
        assert!(g
            .log_suffix()
            .contains("2 from 3 sequences (2 with a taxid, 1 without"));
    }

    #[test]
    fn map_matches_first_header_token_when_full_header_absent() {
        // ">g1|c1 first contig" is stored with its description; the map names the token.
        let idx = index();
        let g = Grouping::from_pairs(&idx, &pairs(&[("g1|c1", 42)]));
        assert_eq!(g.name(g.group(SeqId::new(0))), "42");
    }

    #[test]
    fn map_entries_without_a_sequence_are_ignored() {
        let idx = index();
        let g = Grouping::from_pairs(&idx, &pairs(&[("nope", 1)]));
        assert!(g.is_identity());
        assert_eq!(g.n_unmapped(), 3);
        // The flag was given, so tables still label the column by taxid.
        assert_eq!(g.label_header(), "TaxID");
    }

    #[test]
    fn parse_map_splits_on_whitespace_skips_comments_ignores_extra_columns() {
        let m = parse_map("# c\nA 1 extra\n\n  # indented comment\nB\t1\nA\t1\n").unwrap();
        assert_eq!(m.len(), 2);
        assert_eq!(m["A"], 1);
        assert_eq!(m["B"], 1);
    }

    #[test]
    fn parse_map_rejects_conflicts_nonnumeric_and_zero() {
        let err = parse_map("A\t1\nA\t2\n").unwrap_err().to_string();
        assert!(err.contains("mapped to both"), "{err}");
        let err = parse_map("A\tfoo\n").unwrap_err().to_string();
        assert!(err.contains("not a positive integer"), "{err}");
        assert!(parse_map("A\t0\n").is_err());
        assert!(parse_map("A\t-3\n").is_err());
        assert!(parse_map("A\n").is_err());
    }
}
