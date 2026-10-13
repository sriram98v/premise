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
    /// `Some(rho)` runs with the penalty at that strength; `None` passes `--no-penalty`.
    rho: Option<&'a str>,
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
    match opts.rho {
        Some(rho) => args.extend(["--rho", rho, "--omega", "1e-10"].map(str::to_string)),
        None => args.push("--no-penalty".to_string()),
    }
    if let Some(m) = opts.map {
        args.push("--seq2taxid".to_string());
        args.push(m.to_str().unwrap().to_string());
    }
    if let Some(n) = opts.nodes {
        args.push("--nodes".to_string());
        args.push(n.to_str().unwrap().to_string());
    }
    if let Some(g) = opts.gamma {
        args.push("--gamma".to_string());
        args.push(g.to_string());
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

/// `.nodes` rows keyed by TaxID: (Parent, Rank, Contigs, Count, Proportion, Conditional).
fn nodes_table(prefix: &Path) -> HashMap<String, (String, String, usize, f64, f64, f64)> {
    read(prefix, "nodes")
        .lines()
        .skip(1)
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            assert_eq!(c.len(), 7, "{l}");
            (
                c[0].to_string(),
                (
                    c[1].to_string(),
                    c[2].to_string(),
                    c[3].parse().unwrap(),
                    c[4].parse().unwrap(),
                    c[5].parse().unwrap(),
                    c[6].parse().unwrap(),
                ),
            )
        })
        .collect()
}

fn props(prefix: &Path) -> HashMap<String, f64> {
    read(prefix, "props")
        .lines()
        .map(|l| {
            let (n, p) = l.split_once('\t').unwrap();
            (n.to_string(), p.parse().unwrap())
        })
        .collect()
}

/// Assigned ids of `.matches`, excluding unclassified reads.
fn matches(prefix: &Path) -> Vec<String> {
    read(prefix, "matches")
        .lines()
        .skip(1)
        .filter_map(|l| l.split('\t').nth(1))
        .filter(|r| *r != "unclassified")
        .map(str::to_string)
        .collect()
}

fn posterior_labels(prefix: &Path) -> Vec<String> {
    read(prefix, "posteriors")
        .lines()
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
fn gamma_outside_half_one_is_an_error() {
    let tmp = tempfile::tempdir().unwrap();
    let index = build_index(tmp.path());
    for g in ["0", "0.3", "0.5", "1.5"] {
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
        gamma: Some("0.51"),
        ..Default::default()
    };
    query(tmp.path(), &index, "fine", &fine);
}

#[test]
fn plain_query_writes_star_nodes_table() {
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
    assert_eq!(read(&flat_out, "aligns"), read(&tree_out, "aligns"));
    assert_eq!(matches(&flat_out), matches(&tree_out));
    assert!(read(&tree_out, "posteriors").starts_with("ReadID\tTaxID\tPosterior"));
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
    for id in posterior_labels(&out) {
        assert!(id == "ref_A" || id == "ref_B", "{id}");
    }
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

    for gamma in ["1", "0.51"] {
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
    assert!(posterior_labels(&out).iter().all(|i| i == "11"));
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
