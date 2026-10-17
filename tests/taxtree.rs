//! `query --nodes`: the hierarchical EM over a taxonomy whose leaves are the index
//! sequences, contraction, pruning / collapse and the `--gamma` read assignment.

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

fn build_index_from(dir: &Path, fasta: &Path) -> PathBuf {
    let out = dir.join(format!(
        "{}.fmidx",
        fasta.file_stem().unwrap().to_str().unwrap()
    ));
    ok(&run(&[
        "build",
        "-s",
        fasta.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]));
    out
}

fn build_index(dir: &Path) -> PathBuf {
    build_index_from(dir, &fixtures().join("ref.fasta"))
}

#[derive(Default)]
struct Opts<'a> {
    map: Option<&'a Path>,
    nodes: Option<&'a Path>,
    gamma: Option<&'a str>,
    /// `Some(rho)` runs with the penalty at that strength; `None` without (the default).
    rho: Option<&'a str>,
    init: Option<&'a str>,
}

fn query_args(index: &Path, out: &Path, opts: &Opts<'_>) -> Vec<String> {
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
        "20",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    if let Some(rho) = opts.rho {
        args.extend(["--rho", rho, "--omega", "1e-10"].map(str::to_string));
    }
    if let Some(m) = opts.map {
        args.push("--seq2taxid".to_string());
        args.push(m.to_str().unwrap().to_string());
    }
    if let Some(n) = opts.nodes {
        args.push("--nodes".to_string());
        args.push(n.to_str().unwrap().to_string());
    }
    // `query` writes `.matches` only with --matches; the tests read it.
    args.push("--matches".to_string());
    args.push("--gamma".to_string());
    args.push(opts.gamma.unwrap_or("0.6").to_string());
    if let Some(i) = opts.init {
        args.push("--init".to_string());
        args.push(i.to_string());
    }
    args
}

fn query(
    dir: &Path,
    index: &Path,
    prefix: &str,
    opts: &Opts<'_>,
) -> (PathBuf, std::process::Output) {
    let out = dir.join(prefix);
    let args = query_args(index, &out, opts);
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

/// The report's NODES rows with their TAXLABELS label: `(label, row tokens)`.
fn nodes_rows(prefix: &Path) -> Vec<(String, Vec<String>)> {
    let labels = taxlabels(&block(prefix, "TAXA"));
    let rows = matrix_rows(&block(prefix, "NODES"));
    assert_eq!(labels.len(), rows.len());
    assert!(rows.iter().all(|r| r.len() == 8), "{rows:?}");
    labels.into_iter().zip(rows).collect()
}

/// Status of every vertex of the report's NODES block, by label (labels come from
/// TAXLABELS, in NODES order).
fn vertex_status(prefix: &Path) -> HashMap<String, String> {
    nodes_rows(prefix)
        .into_iter()
        .map(|(label, r)| (label, r[3].clone()))
        .collect()
}

/// The final tree of the report's NODES block (the `active` and `collapsed` vertices),
/// keyed by label: (parent label or `-`, rank, contigs, count, proportion, conditional).
fn nodes_table(prefix: &Path) -> HashMap<String, (String, String, usize, f64, f64, f64)> {
    let rows = nodes_rows(prefix);
    let label_of: HashMap<&str, &str> = rows
        .iter()
        .map(|(label, r)| (r[0].as_str(), label.as_str()))
        .collect();
    let num = |v: &str| -> f64 {
        if v == "-" {
            0.0
        } else {
            v.parse().unwrap()
        }
    };
    rows.iter()
        .filter(|(_, r)| r[3] == "active" || r[3] == "collapsed")
        .map(|(label, r)| {
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
                    num(&r[5]),
                    num(&r[6]),
                    num(&r[7]),
                ),
            )
        })
        .collect()
}

/// Every row of the NEXUS `.props` tree in pre-order: (label, parent label or `-`, rank,
/// proportion). Checks the shape on the way: TAXLABELS and DIMENSIONS agree with the
/// rows, the first row is the root with all the mass, every other row's parent is an
/// earlier row, and no proportion exceeds its parent's.
fn props_tree(prefix: &Path) -> Vec<(String, String, String, f64)> {
    let text = read(prefix, "props");
    assert!(
        text.starts_with("#NEXUS\n[PREMISE proportions v1]\n"),
        "{text}"
    );
    let section = |name: &str| -> String {
        let start = text
            .find(&format!("BEGIN {name};"))
            .unwrap_or_else(|| panic!("no {name} block"));
        let end = text[start..].find("END;").expect("block end") + start;
        text[start..end].to_string()
    };
    let rows: Vec<(String, String, String, f64)> = matrix_rows(&section("PROPORTIONS"))
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
    let labels = taxlabels(&section("TAXA"));
    assert_eq!(
        labels,
        rows.iter().map(|r| r.0.clone()).collect::<Vec<_>>(),
        "{text}"
    );
    assert!(
        section("PROPORTIONS").contains(&format!("DIMENSIONS N={};", rows.len())),
        "{text}"
    );
    let (root, rest) = rows.split_first().expect("a root row");
    assert_eq!((root.1.as_str(), root.2.as_str()), ("-", "root"), "{text}");
    assert!((root.3 - 1.0).abs() < 1e-9, "{text}");
    let mut seen: HashMap<&str, f64> = HashMap::from([(root.0.as_str(), root.3)]);
    for (label, parent, _, p) in rest {
        let up = seen
            .get(parent.as_str())
            .unwrap_or_else(|| panic!("{label}: parent {parent} not an earlier row\n{text}"));
        assert!(
            *p > 0.0 && *p <= up + 1e-9,
            "{label}: {p} under {up}\n{text}"
        );
        seen.insert(label.as_str(), *p);
    }
    rows
}

/// The leaves of the `.props` tree (the rows that are no other row's parent): label →
/// proportion, i.e. the final EM leaves with mass.
fn props(prefix: &Path) -> HashMap<String, f64> {
    let rows = props_tree(prefix);
    rows.iter()
        .filter(|(label, ..)| !rows.iter().any(|r| r.1 == *label))
        .map(|(label, _, _, p)| (label.clone(), *p))
        .collect()
}

/// The rows of a `.matches` table after its leading `#` EM summary: the column header,
/// then one row per read.
fn table_lines(table: &str) -> impl Iterator<Item = &str> {
    table.lines().filter(|l| !l.starts_with('#'))
}

/// Assigned ids of `.matches`, excluding unclassified reads.
fn matches(prefix: &Path) -> Vec<String> {
    table_lines(&read(prefix, "matches"))
        .skip(1)
        .filter_map(|l| l.split('\t').nth(1))
        .filter(|r| *r != "unclassified")
        .map(str::to_string)
        .collect()
}

fn nodes_fixture() -> PathBuf {
    fixtures().join("nodes.dmp")
}

#[test]
fn nodes_requires_seq2taxid() {
    let tmp = tempfile::tempdir().unwrap();
    let index = build_index(tmp.path());
    let out = tmp.path().join("bad");
    let nodes = nodes_fixture();
    let opts = Opts {
        nodes: Some(&nodes),
        ..Default::default()
    };
    let args = query_args(&index, &out, &opts);
    let argv: Vec<&str> = args.iter().map(String::as_str).collect();
    let output = run(&argv);
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("--seq2taxid"), "{stderr}");
}

#[test]
fn gamma_outside_half_open_unit_interval_is_an_error() {
    let tmp = tempfile::tempdir().unwrap();
    let index = build_index(tmp.path());
    for g in ["0", "0.3", "0.49", "1", "1.5"] {
        let out = tmp.path().join("bad");
        let opts = Opts {
            gamma: Some(g),
            ..Default::default()
        };
        let args = query_args(&index, &out, &opts);
        let argv: Vec<&str> = args.iter().map(String::as_str).collect();
        let output = run(&argv);
        assert!(!output.status.success(), "gamma {g}");
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(stderr.contains("gamma"), "{stderr}");
    }
    let fine = Opts {
        gamma: Some("0.5"),
        ..Default::default()
    };
    query(tmp.path(), &index, "fine", &fine);
}

#[test]
fn plain_query_fits_a_star_tree() {
    let tmp = tempfile::tempdir().unwrap();
    let index = build_index(tmp.path());
    let (out, output) = query(tmp.path(), &index, "plain", &Opts::default());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(!stdout.contains("Taxonomy:"), "{stdout}");
    let nodes = nodes_table(&out);
    assert_eq!(nodes.len(), 3, "{nodes:?}");
    let root = &nodes["root"];
    assert_eq!(root.0, "-");
    assert_eq!(root.1, "root");
    assert!((root.4 - 1.0).abs() < 1e-9);
    assert_eq!(root.2, 0);
    for r in ["ref_A", "ref_B"] {
        assert_eq!(nodes[r].0, "root");
        assert_eq!(nodes[r].1, "sequence");
        assert_eq!(nodes[r].2, 1);
        assert!(
            (nodes[r].4 - nodes[r].5).abs() < 1e-12,
            "root conditional = marginal"
        );
    }
    let p = props(&out);
    assert!((p.values().sum::<f64>() - 1.0).abs() < 1e-6, "{p:?}");
    assert!((p["ref_A"] - nodes["ref_A"].4).abs() < 1e-9);
    // .props is the same star: the root over the two sequences.
    let tree: Vec<(String, String, String)> = props_tree(&out)
        .into_iter()
        .map(|(l, parent, rank, _)| (l, parent, rank))
        .collect();
    let row = |l: &str, p: &str, r: &str| (l.to_string(), p.to_string(), r.to_string());
    assert_eq!(
        tree,
        vec![
            row("root", "-", "root"),
            row("ref_A", "root", "sequence"),
            row("ref_B", "root", "sequence"),
        ]
    );
}

#[test]
fn sibling_species_share_genus() {
    let tmp = tempfile::tempdir().unwrap();
    let index = build_index(tmp.path());
    let map = tmp.path().join("sib.map");
    fs::write(&map, "ref_A\t11\nref_B\t12\n").unwrap();
    let nodes = nodes_fixture();

    let flat = Opts {
        map: Some(&map),
        ..Default::default()
    };
    let (flat_out, flat_output) = query(tmp.path(), &index, "flat", &flat);
    let stdout = String::from_utf8_lossy(&flat_output.stdout);
    assert!(
        stdout.contains("Taxonomy: 5 vertices, 2 sequences at tips under their taxids (0 sequences without taxid)"),
        "{stdout}"
    );
    let tree = Opts {
        map: Some(&map),
        nodes: Some(&nodes),
        ..Default::default()
    };
    let (tree_out, output) = query(tmp.path(), &index, "tree", &tree);

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        stdout.contains(
            "Taxonomy: 7 vertices, 2 sequences at tips (0 taxids not in the nodes file, 0 sequences without taxid)"
        ),
        "{stdout}"
    );
    assert!(
        stdout.contains("Contracted taxonomy: 7 vertices, 2 tips with reads"),
        "{stdout}"
    );

    let t = nodes_table(&tree_out);
    assert_eq!(t.len(), 7, "root paths are kept whole: {t:?}");
    assert_eq!(t["root"].0, "-");
    assert_eq!(t["1"].0, "root");
    assert_eq!(t["10"].0, "1");
    assert_eq!(t["10"].1, "genus");
    assert!((t["10"].4 - 1.0).abs() < 1e-9);
    assert_eq!(t["11"].0, "10");
    assert_eq!(t["12"].1, "species");
    assert_eq!(t["ref_A"].0, "11");
    assert_eq!(t["ref_A"].1, "sequence");
    assert_eq!(t["ref_A"].2, 1);
    assert_eq!(t["11"].2, 0, "a taxid holds no sequence directly");
    assert!((t["11"].4 + t["12"].4 - t["10"].4).abs() < 1e-9);
    assert!((t["ref_A"].4 - t["11"].4).abs() < 1e-12, "unary edge");
    assert!(
        (t["11"].3 + t["12"].3 - t["10"].3).abs() < 1e-6,
        "counts add up"
    );

    // with rho = 0 the tree fit is the fit without the lineages: the leaves are the same
    let (pf, pt) = (props(&flat_out), props(&tree_out));
    assert_eq!(pf.len(), 2);
    assert!(
        pf.contains_key("ref_A") && pf.contains_key("ref_B"),
        "{pf:?}"
    );
    for (k, v) in &pf {
        assert!((pt[k] - v).abs() < 1e-6, "{k}: {pt:?} vs {pf:?}");
    }
    assert_eq!(aligns_no_pos(&flat_out), aligns_no_pos(&tree_out));
    assert_eq!(matches(&flat_out), matches(&tree_out));
    assert!(block(&tree_out, "PREMISE").contains("LABEL_HEADER TaxID;"));
}

#[test]
fn sequences_of_an_internal_taxid_are_its_own_tips() {
    let tmp = tempfile::tempdir().unwrap();
    let index = build_index(tmp.path());
    let map = tmp.path().join("genus.map");
    // ref_A is resolved to the genus itself, ref_B to a species below it
    fs::write(&map, "ref_A\t10\nref_B\t11\n").unwrap();
    let nodes = nodes_fixture();
    let opts = Opts {
        map: Some(&map),
        nodes: Some(&nodes),
        ..Default::default()
    };
    let (out, _) = query(tmp.path(), &index, "genus", &opts);

    let t = nodes_table(&out);
    assert_eq!(t.len(), 6, "root, 1, 10, 11, ref_A, ref_B: {t:?}");
    assert!(t.values().all(|r| r.1 != "faux"));
    assert_eq!(t["ref_A"].0, "10", "ref_A is a leaf child of the genus");
    assert_eq!(t["ref_B"].0, "11");
    assert_eq!(t["10"].2, 0);
    assert!((t["10"].4 - 1.0).abs() < 1e-9);
    assert!((t["ref_A"].4 + t["11"].4 - t["10"].4).abs() < 1e-9);

    let p = props(&out);
    assert_eq!(p.len(), 2, "{p:?}");
    assert!((p["ref_A"] + p["ref_B"] - 1.0).abs() < 1e-6);
    assert!((p["ref_B"] - t["ref_B"].4).abs() < 1e-9);
    for id in matches(&out) {
        assert!(
            ["10", "11", "ref_A", "ref_B"].contains(&id.as_str()),
            "{id}"
        );
    }
    let status = vertex_status(&out);
    assert_eq!(status["ref_A"], "active");
    assert_eq!(status["ref_B"], "active");
    assert_eq!(status.len(), 6, "{status:?}");
}

#[test]
fn contraction_drops_reference_without_reads() {
    let tmp = tempfile::tempdir().unwrap();
    let mut fasta = fs::read_to_string(fixtures().join("ref.fasta")).unwrap();
    fasta.push_str(
        ">ref_C\nGGGTTTCCCAAAGGGTTTCCCAAAGGGTTTCCCAAAGGGTTTCCCAAAGGGTTTCCCAAAGGGTTTCCCAAA\n",
    );
    let fasta_path = tmp.path().join("three.fasta");
    fs::write(&fasta_path, fasta).unwrap();
    let index = build_index_from(tmp.path(), &fasta_path);
    let map = tmp.path().join("three.map");
    fs::write(&map, "ref_A\t11\nref_B\t12\nref_C\t13\n").unwrap();
    let mut nodes_text = fs::read_to_string(nodes_fixture()).unwrap();
    nodes_text.push_str("13\t|\t10\t|\tspecies\t|\n");
    let nodes = tmp.path().join("nodes.dmp");
    fs::write(&nodes, nodes_text).unwrap();
    let opts = Opts {
        map: Some(&map),
        nodes: Some(&nodes),
        ..Default::default()
    };
    let (out, output) = query(tmp.path(), &index, "three", &opts);
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        stdout.contains("Taxonomy: 9 vertices, 3 sequences at tips"),
        "{stdout}"
    );
    assert!(
        stdout.contains("Contracted taxonomy: 7 vertices, 2 tips with reads"),
        "{stdout}"
    );
    let t = nodes_table(&out);
    assert!(!t.contains_key("13") && !t.contains_key("ref_C"), "{t:?}");
    assert!(!props(&out).contains_key("ref_C"));
}

#[test]
fn ties_never_descend_below_the_genus() {
    let tmp = tempfile::tempdir().unwrap();
    // two identical references: every read is a perfect tie between species 11 and 12
    let ref_a = fs::read_to_string(fixtures().join("ref.fasta")).unwrap();
    let seq_a: String = ref_a
        .lines()
        .skip(1)
        .take_while(|l| !l.starts_with('>'))
        .collect();
    let fasta_path = tmp.path().join("twins.fasta");
    fs::write(&fasta_path, format!(">twin_1\n{seq_a}\n>twin_2\n{seq_a}\n")).unwrap();
    let index = build_index_from(tmp.path(), &fasta_path);
    let map = tmp.path().join("twins.map");
    fs::write(&map, "twin_1\t11\ntwin_2\t12\n").unwrap();
    let nodes = nodes_fixture();

    for gamma in ["0.99", "0.5"] {
        let opts = Opts {
            map: Some(&map),
            nodes: Some(&nodes),
            gamma: Some(gamma),
            ..Default::default()
        };
        let (out, _) = query(tmp.path(), &index, &format!("g{gamma}"), &opts);
        let ids = matches(&out);
        assert!(!ids.is_empty());
        assert!(ids.iter().all(|i| i == "10"), "gamma {gamma}: {ids:?}");
        let t = nodes_table(&out);
        assert!((t["10"].4 - 1.0).abs() < 1e-9);
        assert!((t["twin_1"].4 - 0.5).abs() < 1e-6, "{t:?}");
    }
}

/// Without a taxonomy the tree is the root over the sequences, and `.matches` descends
/// from the root like any other tree: twin references split every read 0.5 / 0.5, which
/// no gamma in [0.5, 1) descends into, so every read is unclassified. With the twins as
/// sibling species the same reads stop at their genus instead.
#[test]
fn without_a_taxonomy_twins_stay_unclassified_at_every_gamma() {
    let tmp = tempfile::tempdir().unwrap();
    let ref_a = fs::read_to_string(fixtures().join("ref.fasta")).unwrap();
    let seq_a: String = ref_a
        .lines()
        .skip(1)
        .take_while(|l| !l.starts_with('>'))
        .collect();
    let fasta_path = tmp.path().join("twins.fasta");
    fs::write(&fasta_path, format!(">twin_1\n{seq_a}\n>twin_2\n{seq_a}\n")).unwrap();
    let index = build_index_from(tmp.path(), &fasta_path);

    let (out, _) = query(tmp.path(), &index, "plain", &Opts::default());
    assert!(matches(&out).is_empty(), "{}", read(&out, "matches"));
    let p = props(&out);
    assert_eq!(p.len(), 2, "{p:?}");
    assert!((p["twin_1"] - 0.5).abs() < 1e-6 && (p["twin_2"] - 0.5).abs() < 1e-6);
    // the same report at two gammas: byte-identical tables, every read unclassified
    let report = format!("{}.report", out.display());
    let at = |gamma: &str| {
        let o = tmp.path().join(format!("a{gamma}"));
        ok(&run(&[
            "analyze",
            "-r",
            &report,
            "-o",
            o.to_str().unwrap(),
            "--matches",
            "--gamma",
            gamma,
        ]));
        read(&o, "matches")
    };
    let (lo, hi) = (at("0.5"), at("0.99"));
    // The same EM summary but for the gamma it names, then identical tables.
    let summary = |t: &str| t.lines().next().unwrap().to_string();
    let (s_lo, s_hi) = (summary(&lo), summary(&hi));
    assert!(
        s_lo.starts_with("# PREMISE ") && s_lo.contains(" EM: init="),
        "{s_lo}"
    );
    assert_eq!(
        s_lo.strip_suffix("; gamma=0.5"),
        s_hi.strip_suffix("; gamma=0.99"),
        "{s_lo}\n{s_hi}"
    );
    assert!(s_lo.strip_suffix("; gamma=0.5").is_some(), "{s_lo}");
    let body = |t: &str| table_lines(t).collect::<Vec<_>>().join("\n");
    assert_eq!(body(&lo), body(&hi));
    assert_eq!(lo.lines().count(), 7, "{lo}");
    for line in table_lines(&lo).skip(1) {
        assert!(line.ends_with("\tunclassified\t-\t-\t-"), "{line}");
    }

    // the twins as sibling species of genus 10: the reads of the twins go to the genus
    let map = tmp.path().join("twins.map");
    fs::write(&map, "twin_1\t11\ntwin_2\t12\n").unwrap();
    let nodes = nodes_fixture();
    for gamma in ["0.5", "0.99"] {
        let opts = Opts {
            map: Some(&map),
            nodes: Some(&nodes),
            gamma: Some(gamma),
            ..Default::default()
        };
        let (out, _) = query(tmp.path(), &index, &format!("tree{gamma}"), &opts);
        let ids = matches(&out);
        // reads 4 and 5 come from ref_B, absent from the twins
        assert_eq!(ids, vec!["10"; 3], "gamma {gamma}");
    }
}

#[test]
fn penalty_prunes_a_weak_species_and_reports_it() {
    let tmp = tempfile::tempdir().unwrap();
    let index = build_index(tmp.path());
    let map = tmp.path().join("sib.map");
    fs::write(&map, "ref_A\t11\nref_B\t12\n").unwrap();
    let nodes = nodes_fixture();
    // 3 reads hit ref_A, 2 hit ref_B; rho = 2.5 prunes species 12 and its sequence
    let opts = Opts {
        map: Some(&map),
        nodes: Some(&nodes),
        rho: Some("2.5"),
        ..Default::default()
    };
    let (out, output) = query(tmp.path(), &index, "pruned", &opts);
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("Pruned 2 vertices"), "{stdout}");
    let t = nodes_table(&out);
    assert!(!t.contains_key("12") && !t.contains_key("ref_B"), "{t:?}");
    assert_eq!(t["ref_A"].0, "11");
    let p = props(&out);
    assert_eq!(p.keys().collect::<Vec<_>>(), vec!["ref_A"], "{p:?}");
    assert!((p["ref_A"] - 1.0).abs() < 1e-9);
    assert!(matches(&out).iter().all(|i| i == "ref_A"));
}

#[test]
fn collapse_reports_the_taxid_over_its_sequences() {
    let tmp = tempfile::tempdir().unwrap();
    let index = build_index(tmp.path());
    let map = tmp.path().join("same.map");
    fs::write(&map, "ref_A\t11\nref_B\t11\n").unwrap();
    let nodes = nodes_fixture();
    // 3 + 2 reads: each sequence is below rho = 3.5, the species above it
    let opts = Opts {
        map: Some(&map),
        nodes: Some(&nodes),
        rho: Some("3.5"),
        ..Default::default()
    };
    let (out, output) = query(tmp.path(), &index, "collapsed", &opts);
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        stdout.contains("collapsed 1 into tips; 1 tips remain"),
        "{stdout}"
    );
    let t = nodes_table(&out);
    assert_eq!(t.len(), 4, "root, 1, 10, 11: {t:?}");
    assert_eq!(t["11"].2, 2, "the collapsed species pools both sequences");
    assert_eq!(t["11"].1, "species");
    assert!((t["11"].4 - 1.0).abs() < 1e-9);
    let p = props(&out);
    assert_eq!(p.keys().collect::<Vec<_>>(), vec!["11"], "{p:?}");
    let ids = matches(&out);
    assert!(!ids.is_empty());
    assert!(ids.iter().all(|i| i == "11"), "{ids:?}");
    let status = vertex_status(&out);
    assert_eq!(status["11"], "collapsed");
    assert_eq!(status["ref_A"], "pooled");
    assert_eq!(status["ref_B"], "pooled");
}

#[test]
fn uniform_init_matches_split_init_without_penalty() {
    let tmp = tempfile::tempdir().unwrap();
    let index = build_index(tmp.path());
    let map = tmp.path().join("sib.map");
    fs::write(&map, "ref_A\t11\nref_B\t12\n").unwrap();
    let nodes = nodes_fixture();
    let split = Opts {
        map: Some(&map),
        nodes: Some(&nodes),
        ..Default::default()
    };
    let (s_out, s_output) = query(tmp.path(), &index, "split", &split);
    assert!(String::from_utf8_lossy(&s_output.stdout).contains("Init: split"));
    let unif = Opts {
        map: Some(&map),
        nodes: Some(&nodes),
        init: Some("uniform"),
        ..Default::default()
    };
    let (u_out, u_output) = query(tmp.path(), &index, "unif", &unif);
    assert!(String::from_utf8_lossy(&u_output.stdout).contains("Init: uniform"));
    let (ps, pu) = (props(&s_out), props(&u_out));
    assert_eq!(ps.len(), pu.len());
    for (k, v) in &ps {
        assert!((pu[k] - v).abs() < 1e-4, "{k}: {pu:?} vs {ps:?}");
    }
    assert_eq!(matches(&s_out), matches(&u_out));
    let bad = Opts {
        init: Some("random"),
        ..Default::default()
    };
    let args = query_args(&index, &tmp.path().join("bad"), &bad);
    let argv: Vec<&str> = args.iter().map(String::as_str).collect();
    assert!(!run(&argv).status.success());
}

#[test]
fn root_collapse_warns_and_classifies_nothing() {
    let tmp = tempfile::tempdir().unwrap();
    let index = build_index(tmp.path());
    let map = tmp.path().join("sib.map");
    fs::write(&map, "ref_A\t11\nref_B\t12\n").unwrap();
    let nodes = nodes_fixture();
    // rho above the number of reads: every child of the root fails
    let opts = Opts {
        map: Some(&map),
        nodes: Some(&nodes),
        rho: Some("10"),
        ..Default::default()
    };
    let (out, output) = query(tmp.path(), &index, "root", &opts);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("rho is too large"), "{stderr}");
    assert!(matches(&out).is_empty());
    let t = nodes_table(&out);
    assert_eq!(t.len(), 1, "{t:?}");
    assert_eq!(t["root"].2, 2);
    let p = props(&out);
    assert_eq!(p.keys().collect::<Vec<_>>(), vec!["root"], "{p:?}");
}
