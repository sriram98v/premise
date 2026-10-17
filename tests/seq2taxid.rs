//! `query --seq2taxid`: every taxid becomes an internal vertex over its index sequences.

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
}

fn run(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_premise"))
        .args(args)
        .output()
        .expect("failed to run premise")
}

fn ok(output: &std::process::Output) {
    assert!(
        output.status.success(),
        "exit {:?}\nstdout: {}\nstderr: {}",
        output.status.code(),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

fn build_index(dir: &Path) -> PathBuf {
    let out = dir.join("ref.fmidx");
    ok(&run(&[
        "build",
        "-s",
        fixtures().join("ref.fasta").to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]));
    out
}

fn query_args(index: &Path, out: &Path, map: Option<&Path>) -> Vec<String> {
    let r1 = fixtures().join("r1.fastq");
    let r2 = fixtures().join("r2.fastq");
    let mut args: Vec<String> = [
        "query",
        "-s",
        index.to_str().unwrap(),
        "-1",
        r1.to_str().unwrap(),
        "-2",
        r2.to_str().unwrap(),
        "-m",
        "5",
        "-o",
        out.to_str().unwrap(),
        "-t",
        "1",
        "-i",
        "5",
        "--matches",
        "--gamma",
        "0.6",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    if let Some(m) = map {
        args.push("--seq2taxid".to_string());
        args.push(m.to_str().unwrap().to_string());
    }
    args
}

/// Run a query with fixed EM settings; `map` adds `--seq2taxid`. Returns the output
/// prefix and the process output (for stderr checks).
fn query(
    dir: &Path,
    index: &Path,
    prefix: &str,
    map: Option<&Path>,
) -> (PathBuf, std::process::Output) {
    let out = dir.join(prefix);
    let args = query_args(index, &out, map);
    let argv: Vec<&str> = args.iter().map(String::as_str).collect();
    let output = run(&argv);
    ok(&output);
    (out, output)
}

fn read(prefix: &Path, ext: &str) -> String {
    fs::read_to_string(format!("{}.{}", prefix.display(), ext)).expect("output file")
}

/// One block of the report, from `BEGIN <name>;` to its `END;`.
fn block(prefix: &Path, name: &str) -> String {
    let report = read(prefix, "report");
    let start = report
        .find(&format!("BEGIN {name};"))
        .unwrap_or_else(|| panic!("no {name} block"));
    let end = report[start..].find("END;").expect("block end") + start;
    report[start..end].to_string()
}

/// The ALIGNMENT_LIKELIHOODS block without the two position lines of every read (tied
/// placements are not deterministic between runs).
fn aligns_no_pos(prefix: &Path) -> String {
    let block = block(prefix, "ALIGNMENT_LIKELIHOODS");
    let mut out = Vec::new();
    let mut lines = block.lines().skip_while(|l| !l.contains("MATRIX")).skip(1);
    while let Some(hdr) = lines.next() {
        out.push(hdr.trim().to_string());
        if hdr.trim() == ";" {
            break;
        }
        let n: usize = hdr
            .split_whitespace()
            .last()
            .and_then(|v| v.parse().ok())
            .unwrap_or(0);
        if n > 0 {
            out.push(lines.next().unwrap().trim().to_string());
            out.push(lines.next().unwrap().trim().to_string());
            lines.next();
            lines.next();
        }
    }
    out.join("\n")
}

/// Split a NEXUS row into tokens, keeping `'quoted labels'` together (`''` inside quotes
/// is a literal quote).
fn nexus_tokens(line: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut quoted = false;
    let mut chars = line.chars().peekable();
    while let Some(ch) = chars.next() {
        match ch {
            '\'' if quoted && chars.peek() == Some(&'\'') => {
                chars.next();
                cur.push('\'');
            }
            '\'' => quoted = !quoted,
            c if c.is_whitespace() && !quoted => {
                if !cur.is_empty() {
                    out.push(std::mem::take(&mut cur));
                }
            }
            c => cur.push(c),
        }
    }
    if !cur.is_empty() {
        out.push(cur);
    }
    out
}

/// The token rows of a block's MATRIX, up to its closing `;`.
fn matrix_rows(block: &str) -> Vec<Vec<String>> {
    block
        .lines()
        .skip_while(|l| !l.contains("MATRIX"))
        .skip(1)
        .take_while(|l| l.trim() != ";")
        .map(nexus_tokens)
        .collect()
}

/// The TAXLABELS of a TAXA block, in order.
fn taxlabels(block: &str) -> Vec<String> {
    block
        .lines()
        .skip_while(|l| !l.contains("TAXLABELS"))
        .skip(1)
        .take_while(|l| l.trim() != ";")
        .flat_map(nexus_tokens)
        .collect()
}

/// The final tree of the report's NODES block (the `active` and `collapsed` vertices),
/// keyed by label: (parent label or `-`, rank, contigs, proportion).
fn nodes_table(prefix: &Path) -> HashMap<String, (String, String, usize, f64)> {
    let labels = taxlabels(&block(prefix, "TAXA"));
    let rows = matrix_rows(&block(prefix, "NODES"));
    assert_eq!(labels.len(), rows.len());
    let label_of: HashMap<&str, &str> = rows
        .iter()
        .zip(&labels)
        .map(|(r, label)| (r[0].as_str(), label.as_str()))
        .collect();
    rows.iter()
        .zip(&labels)
        .filter(|(r, _)| r[3] == "active" || r[3] == "collapsed")
        .map(|(r, label)| {
            assert_eq!(r.len(), 8, "{r:?}");
            let parent = if r[1] == "-" {
                "-"
            } else {
                label_of[r[1].as_str()]
            };
            (
                label.clone(),
                (
                    parent.to_string(),
                    r[2].clone(),
                    r[4].parse().unwrap(),
                    r[6].parse().unwrap(),
                ),
            )
        })
        .collect()
}

/// Every row of the NEXUS `.props` tree in pre-order: (label, parent label or `-`, rank,
/// proportion). Checks that the first row is the root with all the mass and that every
/// other row's parent is an earlier row.
fn props_tree(prefix: &Path) -> Vec<(String, String, String, f64)> {
    let text = read(prefix, "props");
    assert!(
        text.starts_with("#NEXUS\n[PREMISE proportions v1]\n"),
        "{text}"
    );
    let start = text.find("BEGIN PROPORTIONS;").expect("PROPORTIONS block");
    let rows: Vec<(String, String, String, f64)> = matrix_rows(&text[start..])
        .into_iter()
        .map(|r| {
            assert_eq!(r.len(), 4, "{r:?}");
            (
                r[0].clone(),
                r[1].clone(),
                r[2].clone(),
                r[3].parse().unwrap(),
            )
        })
        .collect();
    let (root, rest) = rows.split_first().expect("a root row");
    assert_eq!((root.1.as_str(), root.2.as_str()), ("-", "root"), "{text}");
    assert!((root.3 - 1.0).abs() < 1e-9, "{text}");
    for (i, (label, parent, ..)) in rest.iter().enumerate() {
        assert!(
            rows[..=i].iter().any(|r| r.0 == *parent),
            "{label}: parent {parent} not an earlier row\n{text}"
        );
    }
    rows
}

/// The leaves of the `.props` tree (rows that are no other row's parent): the final EM
/// leaves with mass, label → proportion.
fn props(prefix: &Path) -> HashMap<String, f64> {
    let rows = props_tree(prefix);
    rows.iter()
        .filter(|(label, ..)| !rows.iter().any(|r| r.1 == *label))
        .map(|(label, _, _, p)| (label.clone(), *p))
        .collect()
}

/// Two leaf maps agree on their labels and, up to float summation order, their values.
fn assert_same_props(a: &HashMap<String, f64>, b: &HashMap<String, f64>) {
    assert_eq!(a.len(), b.len(), "{a:?} vs {b:?}");
    for (k, v) in a {
        assert!((b[k] - v).abs() <= 1e-9, "{k}: {a:?} vs {b:?}");
    }
}

/// The rows of a `.matches` table after its leading `#` EM summary: the column header,
/// then one row per read.
fn table_lines(table: &str) -> impl Iterator<Item = &str> {
    table.lines().filter(|l| !l.starts_with('#'))
}

fn ids(table: &str) -> Vec<String> {
    table_lines(table)
        .skip(1)
        .filter_map(|l| l.split('\t').nth(1))
        .filter(|r| *r != "unclassified")
        .map(str::to_string)
        .collect()
}

/// `.matches` without its position columns (tie placements are nondeterministic) or its
/// EM summary (whose full-precision log-likelihood two runs need not share).
fn no_pos(table: &str) -> String {
    table_lines(table)
        .map(|l| l.split('\t').take(3).collect::<Vec<_>>().join("\t"))
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn no_flag_header_is_refid() {
    let tmp = tempfile::tempdir().unwrap();
    let index = build_index(tmp.path());
    let (plain, output) = query(tmp.path(), &index, "plain", None);
    assert!(block(&plain, "PREMISE").contains("LABEL_HEADER RefID;"));
    // One `#` EM summary line, then the header.
    let matches = read(&plain, "matches");
    let (summary, rest) = matches.split_once('\n').unwrap();
    assert!(summary.starts_with("# PREMISE "), "{summary}");
    assert!(
        summary.contains(" EM: init=split iter=5 ") && summary.ends_with("; gamma=0.6"),
        "{summary}"
    );
    assert!(rest.starts_with("ReadID\tRefID\tPosterior\t"), "{matches}");
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(!stdout.contains("Taxa:"), "{stdout}");
}

#[test]
fn distinct_taxids_match_no_map() {
    let tmp = tempfile::tempdir().unwrap();
    let index = build_index(tmp.path());
    let map = tmp.path().join("distinct.map");
    fs::write(&map, "ref_A\t1\nref_B\t2\n").unwrap();

    let (plain, _) = query(tmp.path(), &index, "plain", None);
    let (taxed, output) = query(tmp.path(), &index, "taxed", Some(&map));

    // the taxid vertices are unary and change nothing but the table header
    let relabel = |t: String| t.replace("\tRefID\t", "\tTaxID\t");
    assert_eq!(aligns_no_pos(&plain), aligns_no_pos(&taxed));
    assert_same_props(&props(&plain), &props(&taxed));
    assert_eq!(
        no_pos(&relabel(read(&plain, "matches"))),
        no_pos(&read(&taxed, "matches"))
    );
    let nodes = nodes_table(&taxed);
    let n1 = &nodes["1"];
    assert_eq!(
        (n1.0.as_str(), n1.1.as_str(), n1.2),
        ("root", "-", 0),
        "{nodes:?}"
    );
    let a = &nodes["ref_A"];
    assert_eq!(
        (a.0.as_str(), a.1.as_str(), a.2),
        ("1", "sequence", 1),
        "{nodes:?}"
    );
    // .props holds the taxid vertices too, each over its sequence with the same mass
    let tree = props_tree(&taxed);
    for (seq, taxid) in [("ref_A", "1"), ("ref_B", "2")] {
        let row = tree.iter().find(|r| r.0 == seq).expect("sequence row");
        let up = tree.iter().find(|r| r.0 == taxid).expect("taxid row");
        assert_eq!(row.1, taxid, "{tree:?}");
        assert_eq!(up.1, "root", "{tree:?}");
        assert!((row.3 - up.3).abs() < 1e-12, "unary edge: {tree:?}");
    }
    assert!(block(&taxed, "PREMISE").contains("LABEL_HEADER TaxID;"));
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        stdout.contains("Taxa: 2 from 2 sequences (2 with a taxid, 0 without"),
        "{stdout}"
    );
}

#[test]
fn shared_taxid_is_an_internal_vertex_over_its_sequences() {
    let tmp = tempfile::tempdir().unwrap();
    let index = build_index(tmp.path());
    let map = tmp.path().join("one.map");
    fs::write(&map, "# seqid taxid genome\nref_A 9 extra\nref_B\t9\n").unwrap();

    let (taxed, output) = query(tmp.path(), &index, "taxed", Some(&map));
    assert!(
        output.stderr.is_empty() || !String::from_utf8_lossy(&output.stderr).contains("Warning")
    );

    // the leaves are the sequences; the taxid is the vertex above them
    let props = props(&taxed);
    let mut total = 0.0;
    for (name, prop) in &props {
        assert!(name == "ref_A" || name == "ref_B", "{props:?}");
        total += prop;
    }
    assert_eq!(props.len(), 2, "{props:?}");
    assert!((total - 1.0).abs() < 1e-9);
    let tree = props_tree(&taxed);
    assert!(
        tree.iter()
            .any(|r| r.0 == "9" && r.1 == "root" && (r.3 - 1.0).abs() < 1e-9),
        "{tree:?}"
    );
    let nodes = nodes_table(&taxed);
    let nine = &nodes["9"];
    assert_eq!(
        (nine.0.as_str(), nine.1.as_str(), nine.2),
        ("root", "-", 0),
        "{nodes:?}"
    );
    let b = &nodes["ref_B"];
    assert_eq!(
        (b.0.as_str(), b.1.as_str(), b.2),
        ("9", "sequence", 1),
        "{nodes:?}"
    );
    assert!((nine.3 - 1.0).abs() < 1e-9, "{nine:?}");

    for id in ids(&read(&taxed, "matches")) {
        assert!(id == "9" || id == "ref_A" || id == "ref_B", "{id}");
    }
    // every alignment names one of the two index sequences, each at most once per read
    let aligns = block(&taxed, "ALIGNMENT_LIKELIHOODS");
    let mut lines = aligns.lines().skip_while(|l| !l.contains("MATRIX")).skip(1);
    let mut n_reads = 0;
    while let Some(hdr) = lines.next() {
        if hdr.trim() == ";" {
            break;
        }
        n_reads += 1;
        let n: usize = hdr.split_whitespace().last().unwrap().parse().unwrap();
        assert!(n <= 2, "{hdr}");
        if n > 0 {
            let seqs: Vec<u32> = lines
                .next()
                .unwrap()
                .split_whitespace()
                .map(|t| t.parse().unwrap())
                .collect();
            assert!(seqs[0] <= 1, "{seqs:?}");
            assert!(seqs.iter().skip(1).all(|&d| d == 1), "{seqs:?}");
            for _ in 0..3 {
                lines.next();
            }
        }
    }
    assert_eq!(n_reads, 5);
    assert!(block(&taxed, "REFMAP").contains("DIMENSIONS N=2;"));
}

#[test]
fn partial_map_keeps_unmapped_sequence_as_own_column_and_warns() {
    let tmp = tempfile::tempdir().unwrap();
    let index = build_index(tmp.path());
    let map = tmp.path().join("partial.map");
    fs::write(&map, "ref_A 9\n").unwrap();

    let (taxed, output) = query(tmp.path(), &index, "taxed", Some(&map));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("1 of 2 index sequences have no taxid"),
        "{stderr}"
    );
    let props = props(&taxed);
    for name in props.keys() {
        assert!(name == "ref_A" || name == "ref_B", "{props:?}");
    }
    let nodes = nodes_table(&taxed);
    let a = &nodes["ref_A"];
    assert_eq!((a.0.as_str(), a.1.as_str()), ("9", "sequence"), "{nodes:?}");
    let b = &nodes["ref_B"];
    assert_eq!(
        (b.0.as_str(), b.1.as_str()),
        ("root", "sequence"),
        "{nodes:?}"
    );
    // .props places the unmapped sequence directly under the root too
    let tree = props_tree(&taxed);
    assert!(
        tree.iter().any(|r| r.0 == "ref_B" && r.1 == "root"),
        "{tree:?}"
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        stdout.contains("Taxa: 2 from 2 sequences (1 with a taxid, 1 without"),
        "{stdout}"
    );
}

#[test]
fn conflicting_map_is_an_error() {
    let tmp = tempfile::tempdir().unwrap();
    let index = build_index(tmp.path());
    let map = tmp.path().join("conflict.map");
    fs::write(&map, "ref_A 9\nref_A 10\n").unwrap();

    let out = tmp.path().join("bad");
    let args = query_args(&index, &out, Some(&map));
    let argv: Vec<&str> = args.iter().map(String::as_str).collect();
    let output = run(&argv);
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("mapped to both"), "{stderr}");
}

#[test]
fn non_numeric_or_zero_taxid_is_an_error() {
    let tmp = tempfile::tempdir().unwrap();
    let index = build_index(tmp.path());
    for (name, text) in [("alpha.map", "ref_A x\n"), ("zero.map", "ref_A 0\n")] {
        let map = tmp.path().join(name);
        fs::write(&map, text).unwrap();
        let out = tmp.path().join("bad");
        let args = query_args(&index, &out, Some(&map));
        let argv: Vec<&str> = args.iter().map(String::as_str).collect();
        let output = run(&argv);
        assert!(!output.status.success(), "{name}");
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(
            stderr.contains("not a positive integer"),
            "{name}: {stderr}"
        );
    }
}
