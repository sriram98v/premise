//! `query --seq2taxid`: every taxid becomes an internal vertex over its index sequences.

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
        "--gamma",
        "0.6",
        "--no-penalty",
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
    export_tables(&out);
    (out, output)
}

fn read(prefix: &Path, ext: &str) -> String {
    fs::read_to_string(format!("{}.{}", prefix.display(), ext)).expect("output file")
}

/// Regenerate the legacy `.props` / `.nodes` tables from the report.
fn export_tables(prefix: &Path) {
    let p = prefix.display();
    ok(&run(&[
        "analyze",
        "-r",
        &format!("{p}.report"),
        "-o",
        &p.to_string(),
        "--props",
        "--nodes",
    ]));
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

fn ids(table: &str) -> Vec<String> {
    table
        .lines()
        .skip(1)
        .filter_map(|l| l.split('\t').nth(1))
        .filter(|r| *r != "unclassified")
        .map(str::to_string)
        .collect()
}

/// `.matches` without its position columns (tie placements are nondeterministic).
fn no_pos(table: &str) -> String {
    table
        .lines()
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
    assert!(read(&plain, "matches").starts_with("ReadID\tRefID\tPosterior\t"));
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
    assert_eq!(read(&plain, "props"), read(&taxed, "props"));
    assert_eq!(
        no_pos(&relabel(read(&plain, "matches"))),
        no_pos(&read(&taxed, "matches"))
    );
    let nodes = read(&taxed, "nodes");
    assert!(
        nodes.lines().any(|l| l.starts_with("1\troot\t-\t0\t")),
        "{nodes}"
    );
    assert!(
        nodes
            .lines()
            .any(|l| l.starts_with("ref_A\t1\tsequence\t1\t")),
        "{nodes}"
    );
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
    let props = read(&taxed, "props");
    let mut total = 0.0;
    for row in props.lines() {
        let (name, prop) = row.split_once('\t').unwrap();
        assert!(name == "ref_A" || name == "ref_B", "{props}");
        total += prop.parse::<f64>().unwrap();
    }
    assert_eq!(props.lines().count(), 2, "{props}");
    assert!((total - 1.0).abs() < 1e-9);
    let nodes = read(&taxed, "nodes");
    assert!(
        nodes.lines().any(|l| l.starts_with("9\troot\t-\t0\t")),
        "{nodes}"
    );
    assert!(
        nodes
            .lines()
            .any(|l| l.starts_with("ref_B\t9\tsequence\t1\t")),
        "{nodes}"
    );
    let nine = nodes.lines().find(|l| l.starts_with("9\t")).unwrap();
    let pi: f64 = nine.split('\t').nth(5).unwrap().parse().unwrap();
    assert!((pi - 1.0).abs() < 1e-9, "{nine}");

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
    let props = read(&taxed, "props");
    for row in props.lines() {
        let name = row.split('\t').next().unwrap();
        assert!(name == "ref_A" || name == "ref_B", "{props}");
    }
    let nodes = read(&taxed, "nodes");
    assert!(
        nodes.lines().any(|l| l.starts_with("ref_A\t9\tsequence\t")),
        "{nodes}"
    );
    assert!(
        nodes
            .lines()
            .any(|l| l.starts_with("ref_B\troot\tsequence\t")),
        "{nodes}"
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
