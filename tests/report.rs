//! The report: `analyze` reproduces `query`'s `.matches` at any gamma, with and without
//! compression; `query` writes what its flags ask for; the head records the run without
//! paths.

use sha2::{Digest, Sha256};
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

/// A run configuration: `Some(rho)` with the penalty, `None` without (the default).
struct Cfg {
    tree: bool,
    rho: Option<&'static str>,
}

/// `query` arguments without `--gamma` or `--matches`.
fn query_args(index: &Path, dir: &Path, out: &Path, cfg: &Cfg) -> Vec<String> {
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
    if let Some(rho) = cfg.rho {
        args.extend(["--rho", rho, "--omega", "1e-10"].map(str::to_string));
    }
    if cfg.tree {
        let map = dir.join("same.map");
        fs::write(&map, "ref_A\t11\nref_B\t11\n").unwrap();
        args.extend([
            "--seq2taxid".to_string(),
            map.to_str().unwrap().to_string(),
            "--nodes".to_string(),
            fixtures().join("nodes.dmp").to_str().unwrap().to_string(),
        ]);
    }
    args
}

fn run_owned(args: &[String]) -> std::process::Output {
    let v: Vec<&str> = args.iter().map(String::as_str).collect();
    run(&v)
}

fn read(path: &Path) -> String {
    fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

fn with_ext(prefix: &Path, ext: &str) -> PathBuf {
    PathBuf::from(format!("{}.{ext}", prefix.display()))
}

/// `.matches` without its position columns: tied placements land on any of the
/// equally good positions, so two runs only agree up to them.
fn no_pos(table: &str) -> String {
    table
        .lines()
        .map(|l| l.split('\t').take(3).collect::<Vec<_>>().join("\t"))
        .collect::<Vec<_>>()
        .join("\n")
}

/// Split a NEXUS row into tokens, keeping `'quoted labels'` together.
fn nexus_tokens(line: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut quoted = false;
    for ch in line.chars() {
        match ch {
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

/// One block of a report, from `BEGIN <name>;` to its `END;`.
fn block(report: &str, name: &str) -> String {
    let start = report
        .find(&format!("BEGIN {name};"))
        .unwrap_or_else(|| panic!("no {name} block"));
    let end = report[start..].find("END;").expect("block end") + start;
    report[start..end].to_string()
}

const GAMMAS: [&str; 4] = ["0.5", "0.6", "0.9", "0.99"];

#[test]
fn analyze_reproduces_query_for_every_gamma() {
    let tmp = tempfile::tempdir().unwrap();
    let index = build_index(tmp.path());
    for (tag, cfg) in [
        (
            "flat",
            Cfg {
                tree: false,
                rho: None,
            },
        ),
        (
            "tree",
            Cfg {
                tree: true,
                rho: Some("3.5"),
            },
        ),
        (
            "tree0",
            Cfg {
                tree: true,
                rho: None,
            },
        ),
    ] {
        // A report without .matches, as `query` writes it by default.
        let base = tmp.path().join(format!("{tag}-base"));
        ok(&run_owned(&query_args(&index, tmp.path(), &base, &cfg)));
        assert!(
            !with_ext(&base, "matches").exists(),
            "{tag}: no --matches, no matches"
        );
        let base_report = with_ext(&base, "report");
        for g in GAMMAS {
            let direct = tmp.path().join(format!("{tag}-q{g}"));
            let mut q = query_args(&index, tmp.path(), &direct, &cfg);
            q.extend([
                "--matches".to_string(),
                "--gamma".to_string(),
                g.to_string(),
            ]);
            ok(&run_owned(&q));
            let want = read(&with_ext(&direct, "matches"));
            assert!(want.lines().count() == 6, "header + 5 reads: {want}");

            // From the run's own report: byte-identical.
            let own = tmp.path().join(format!("{tag}-own{g}"));
            ok(&run(&[
                "analyze",
                "-r",
                with_ext(&direct, "report").to_str().unwrap(),
                "--matches",
                "--gamma",
                g,
                "-o",
                own.to_str().unwrap(),
            ]));
            assert_eq!(read(&with_ext(&own, "matches")), want, "{tag} at gamma {g}");

            // From another run's report: the same assignments and posteriors.
            let gen = tmp.path().join(format!("{tag}-g{g}"));
            ok(&run(&[
                "analyze",
                "-r",
                base_report.to_str().unwrap(),
                "--matches",
                "--gamma",
                g,
                "-o",
                gen.to_str().unwrap(),
            ]));
            assert_eq!(
                no_pos(&read(&with_ext(&gen, "matches"))),
                no_pos(&want),
                "{tag} at gamma {g} from the base report"
            );
        }
    }
}

#[test]
fn query_flags_choose_the_outputs() {
    let tmp = tempfile::tempdir().unwrap();
    let index = build_index(tmp.path());
    let cfg = Cfg {
        tree: false,
        rho: None,
    };
    // --no-report with --matches: only the table.
    let only = tmp.path().join("only");
    let mut q = query_args(&index, tmp.path(), &only, &cfg);
    q.extend(["--matches", "--gamma", "0.6", "--no-report"].map(str::to_string));
    let out = run_owned(&q);
    ok(&out);
    assert!(with_ext(&only, "matches").exists());
    assert!(!with_ext(&only, "report").exists());
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        stdout.contains("only.matches") && !stdout.contains("report"),
        "{stdout}"
    );
    // --no-report without --matches is a usage error, before any work is done.
    let bad = tmp.path().join("bad");
    let mut q = query_args(&index, tmp.path(), &bad, &cfg);
    q.push("--no-report".to_string());
    let out = run_owned(&q);
    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("--matches"), "{stderr}");
    assert!(!with_ext(&bad, "matches").exists() && !with_ext(&bad, "report").exists());
    // --gamma out of [0.5, 1).
    for g in ["0.49", "1", "1.5"] {
        let mut q = query_args(&index, tmp.path(), &bad, &cfg);
        q.extend(["--matches", "--gamma", g].map(str::to_string));
        let out = run_owned(&q);
        assert!(!out.status.success(), "gamma {g}");
        assert!(String::from_utf8_lossy(&out.stderr).contains("--gamma must be in [0.5, 1)"));
    }
    // --compress conflicts with --no-report.
    let mut q = query_args(&index, tmp.path(), &bad, &cfg);
    q.extend(["--matches", "--gamma", "0.6", "--no-report", "--compress"].map(str::to_string));
    assert!(!run_owned(&q).status.success());
}

#[test]
fn compressed_report_gives_the_same_matches_and_tables() {
    let tmp = tempfile::tempdir().unwrap();
    let index = build_index(tmp.path());
    let cfg = Cfg {
        tree: true,
        rho: Some("3.5"),
    };
    let plain = tmp.path().join("plain");
    let mut q = query_args(&index, tmp.path(), &plain, &cfg);
    q.extend([
        "--matches".to_string(),
        "--gamma".to_string(),
        "0.9".to_string(),
    ]);
    ok(&run_owned(&q));

    let zst = tmp.path().join("zst");
    let mut q = query_args(&index, tmp.path(), &zst, &cfg);
    q.push("--compress".to_string());
    let out = run_owned(&q);
    ok(&out);
    assert!(String::from_utf8_lossy(&out.stdout).contains("zst.report.zst"));
    let zpath = tmp.path().join("zst.report.zst");
    let bytes = fs::read(&zpath).unwrap();
    assert_eq!(&bytes[..4], &[0x28, 0xB5, 0x2F, 0xFD], "zstd magic");
    assert!(!tmp.path().join("zst.report").exists());
    assert!(!tmp.path().join("zst.matches").exists());

    let from_zst = tmp.path().join("z");
    ok(&run(&[
        "analyze",
        "-r",
        zpath.to_str().unwrap(),
        "--matches",
        "--gamma",
        "0.9",
        "-o",
        from_zst.to_str().unwrap(),
        "--props",
        "--nodes",
    ]));
    // Head-only export needs no --matches.
    ok(&run(&[
        "analyze",
        "-r",
        tmp.path().join("plain.report").to_str().unwrap(),
        "-o",
        tmp.path().join("p").to_str().unwrap(),
        "--props",
        "--nodes",
    ]));
    assert!(!tmp.path().join("p.matches").exists());
    assert_eq!(
        no_pos(&read(&with_ext(&from_zst, "matches"))),
        no_pos(&read(&with_ext(&plain, "matches")))
    );
    assert_eq!(
        read(&tmp.path().join("z.props")),
        read(&tmp.path().join("p.props"))
    );
    assert_eq!(
        read(&tmp.path().join("z.nodes")),
        read(&tmp.path().join("p.nodes"))
    );
    let nodes = read(&tmp.path().join("p.nodes"));
    assert!(nodes.starts_with("TaxID\tParent\tRank\tContigs\tCount\tProportion\tConditional\n"));
    assert!(
        nodes.lines().any(|l| l.starts_with("11\t10\tspecies\t2\t")),
        "{nodes}"
    );
    let props = read(&tmp.path().join("p.props"));
    assert_eq!(props.lines().count(), 1, "{props}");
    assert!(props.starts_with("11\t1.0000000000e0"), "{props}");
}

#[test]
fn report_head_records_the_run_without_paths() {
    let tmp = tempfile::tempdir().unwrap();
    let index = build_index(tmp.path());
    let out = tmp.path().join("meta");
    let cfg = Cfg {
        tree: true,
        rho: Some("3.5"),
    };
    let mut q = query_args(&index, tmp.path(), &out, &cfg);
    q.extend([
        "--matches".to_string(),
        "--gamma".to_string(),
        "0.6".to_string(),
    ]);
    ok(&run_owned(&q));
    let report = read(&tmp.path().join("meta.report"));
    for needle in [
        "#NEXUS\n[PREMISE report v3]\n",
        "BEGIN PREMISE;",
        "\tVERSION 3;",
        "\tLABEL_HEADER TaxID;",
        "\tINDEX ref.fmidx;",
        "\tINDEX_SHA256 ",
        "\tR1 r1.fastq;",
        "\tR2 r2.fastq;",
        "\tSEQ2TAXID same.map;",
        "\tNODES nodes.dmp;",
        "\tTHREADS 1;",
        "rho=3.5e0",
        "init=split penalty=1;",
        "\tREADS total=5 aligned=",
        "collapsed=1",
        "BEGIN EM_TRACE;",
        "BEGIN TAXA;",
        "BEGIN TREES;",
        "TREE taxonomy = [&R] ",
        "BEGIN NODES;",
        "BEGIN REFMAP;",
        "BEGIN PROPS;",
        "BEGIN ALIGNMENT_LIKELIHOODS;",
        "BEGIN PROFILE;",
        "peak_rss_kb=",
    ] {
        assert!(report.contains(needle), "missing {needle:?} in\n{report}");
    }
    // The tree is written once per field: TREES is topology over TAXA indices, the
    // statuses live in NODES.
    assert!(!report.contains("TRANSLATE"), "{report}");
    let nodes =
        &report[report.find("BEGIN NODES;").unwrap()..report.find("BEGIN REFMAP;").unwrap()];
    for status in [" collapsed ", " pooled "] {
        assert!(nodes.contains(status), "missing {status:?} in\n{nodes}");
    }
    // No directory anywhere in the head: the temp dir's path must not appear.
    let head = &report[..report.find("BEGIN ALIGNMENT_LIKELIHOODS;").unwrap()];
    assert!(!head.contains(tmp.path().to_str().unwrap()), "{head}");
    assert!(!head.contains(fixtures().to_str().unwrap()), "{head}");
    // The hash is the file's SHA-256.
    let sha = report
        .lines()
        .find_map(|l| l.trim().strip_prefix("INDEX_SHA256 "))
        .map(|v| v.trim_end_matches(';').to_string())
        .expect("sha line");
    assert_eq!(sha.len(), 64, "{sha}");
    assert!(sha.bytes().all(|b| b.is_ascii_hexdigit()));
    let expect = format!("{:x}", Sha256::digest(fs::read(&index).unwrap()));
    assert_eq!(sha, expect);
    // NODES has 8 columns and no labels; TAXA has as many labels as NODES rows.
    let nodes = block(&report, "NODES");
    let rows: Vec<&str> = nodes
        .lines()
        .skip_while(|l| !l.contains("MATRIX"))
        .skip(1)
        .take_while(|l| l.trim() != ";")
        .collect();
    assert!(rows.iter().all(|r| nexus_tokens(r).len() == 8), "{nodes}");
    let ntax: usize = block(&report, "TAXA")
        .lines()
        .find_map(|l| l.trim().strip_prefix("DIMENSIONS NTAX="))
        .and_then(|v| v.trim_end_matches(';').parse().ok())
        .unwrap();
    assert_eq!(ntax, rows.len());
    // Every read has its alignments: 5 headers `id 0` or `id n k`, then the seqids, k
    // likelihoods, n indices when k > 1, n forward positions, 1 or n insert offsets.
    let al = block(&report, "ALIGNMENT_LIKELIHOODS");
    let lines: Vec<&str> = al
        .lines()
        .skip_while(|l| !l.contains("MATRIX"))
        .skip(1)
        .take_while(|l| l.trim() != ";")
        .collect();
    let count = |i: usize| lines[i].split_whitespace().count();
    let mut i = 0;
    let mut reads = 0;
    let mut rows_total = 0;
    while i < lines.len() {
        let hdr: Vec<usize> = lines[i]
            .split_whitespace()
            .skip(1)
            .map(|t| t.parse().unwrap())
            .collect();
        reads += 1;
        i += 1;
        if hdr == [0] {
            continue;
        }
        let [n, k] = hdr[..] else {
            panic!("bad read header {}", lines[i - 1])
        };
        assert!((1..=n).contains(&k), "{}", lines[i - 1]);
        let mut expect = vec![n, k];
        if k > 1 {
            expect.push(n);
        }
        expect.push(n);
        for (j, want) in expect.iter().enumerate() {
            assert_eq!(count(i + j), *want, "{}", lines[i + j]);
        }
        i += expect.len();
        assert!(count(i) == 1 || count(i) == n, "{}", lines[i]);
        i += 1;
        rows_total += n;
    }
    assert_eq!(reads, 5);
    assert!(rows_total >= 5, "{al}");
    let niter: usize = report
        .lines()
        .find_map(|l| l.trim().strip_prefix("DIMENSIONS NITER="))
        .and_then(|v| v.trim_end_matches(';').parse().ok())
        .expect("NITER");
    assert!(niter >= 2, "{niter}");
}

#[test]
fn no_index_hash_leaves_the_checksum_empty() {
    let tmp = tempfile::tempdir().unwrap();
    let index = build_index(tmp.path());
    let out = tmp.path().join("nohash");
    let cfg = Cfg {
        tree: false,
        rho: None,
    };
    let mut q = query_args(&index, tmp.path(), &out, &cfg);
    q.extend(["--no-index-hash", "--matches", "--gamma", "0.6"].map(str::to_string));
    ok(&run_owned(&q));
    let report = read(&with_ext(&out, "report"));
    assert!(
        report.contains("\tINDEX ref.fmidx;\n\tINDEX_SHA256 -;\n"),
        "{report}"
    );
    let gen = tmp.path().join("nohash-gen");
    ok(&run(&[
        "analyze",
        "-r",
        with_ext(&out, "report").to_str().unwrap(),
        "--matches",
        "--gamma",
        "0.6",
        "-o",
        gen.to_str().unwrap(),
    ]));
    assert_eq!(
        read(&with_ext(&gen, "matches")),
        read(&with_ext(&out, "matches"))
    );
}

#[test]
fn corrupt_reports_leave_no_partial_matches() {
    let tmp = tempfile::tempdir().unwrap();
    let index = build_index(tmp.path());
    let out = tmp.path().join("c");
    let cfg = Cfg {
        tree: false,
        rho: None,
    };
    ok(&run_owned(&query_args(&index, tmp.path(), &out, &cfg)));
    let text = read(&with_ext(&out, "report"));
    let analyze = |report_text: &str, tag: &str| -> std::process::Output {
        let path = tmp.path().join(format!("{tag}.report"));
        fs::write(&path, report_text).unwrap();
        run(&[
            "analyze",
            "-r",
            path.to_str().unwrap(),
            "--matches",
            "--gamma",
            "0.6",
            "-o",
            tmp.path().join(tag).to_str().unwrap(),
        ])
    };
    // A sequence that REFMAP does not know (a huge delta on a read's seqid line).
    let al = text.find("BEGIN ALIGNMENT_LIKELIHOODS;").unwrap();
    let matrix = text[al..].find("\tMATRIX\n").unwrap() + al + 8;
    let hdr_end = text[matrix..].find('\n').unwrap() + matrix + 1;
    let seq_end = text[hdr_end..].find('\n').unwrap() + hdr_end;
    let bad_seq = format!("{} 5000{}", &text[..seq_end], &text[seq_end..]);
    let outp = analyze(&bad_seq, "badseq");
    assert!(!outp.status.success());
    let err = String::from_utf8_lossy(&outp.stderr);
    assert!(
        err.contains("not in REFMAP") || err.contains("seqids"),
        "{err}"
    );
    assert!(!tmp.path().join("badseq.matches").exists());
    // Fewer reads than announced.
    let short = text.replacen("DIMENSIONS NREADS=5;", "DIMENSIONS NREADS=6;", 1);
    let outp = analyze(&short, "short");
    assert!(!outp.status.success());
    assert!(String::from_utf8_lossy(&outp.stderr).contains("announces 6 reads but holds 5"));
    assert!(!tmp.path().join("short.matches").exists());
    // An output naming the input is refused: a refit report under the input's prefix.
    let r = with_ext(&out, "report");
    let outp = run(&[
        "analyze",
        "-r",
        r.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--rho",
        "1",
        "--write-report",
    ]);
    assert!(!outp.status.success());
    assert!(String::from_utf8_lossy(&outp.stderr).contains("is the input report"));
    assert!(read(&r).starts_with("#NEXUS"), "the report survived");
}

#[test]
fn analyze_needs_an_output_and_a_valid_gamma() {
    let tmp = tempfile::tempdir().unwrap();
    let index = build_index(tmp.path());
    let out = tmp.path().join("g");
    let cfg = Cfg {
        tree: false,
        rho: None,
    };
    ok(&run_owned(&query_args(&index, tmp.path(), &out, &cfg)));
    let report = tmp.path().join("g.report");
    let r = report.to_str().unwrap();
    let o = tmp.path().join("o");
    let o = o.to_str().unwrap();
    // -o is required; with nothing to write it is an error before any work.
    let no_out = run(&["analyze", "-r", r, "--props"]);
    assert!(!no_out.status.success());
    assert!(String::from_utf8_lossy(&no_out.stderr).contains("--out"));
    let no_output = run(&["analyze", "-r", r, "-o", o]);
    assert!(!no_output.status.success());
    assert!(String::from_utf8_lossy(&no_output.stderr).contains("nothing to do"));
    let bad_gamma = run(&["analyze", "-r", r, "--matches", "--gamma", "1", "-o", o]);
    assert!(!bad_gamma.status.success());
    assert!(String::from_utf8_lossy(&bad_gamma.stderr).contains("--gamma must be in [0.5, 1)"));
    // --compress needs --write-report; --rho excludes --no-penalty. Each case would
    // otherwise be a valid request, so clap's rule is what fails it.
    for (extra, needle) in [
        (
            &["--compress", "--no-penalty", "--props"][..],
            "--write-report",
        ),
        (
            &["--rho", "1", "--no-penalty", "--write-report"][..],
            "--no-penalty",
        ),
    ] {
        let mut args = vec!["analyze", "-r", r, "-o", o];
        args.extend_from_slice(extra);
        let bad = run(&args);
        assert!(!bad.status.success(), "{extra:?}");
        let stderr = String::from_utf8_lossy(&bad.stderr);
        assert!(stderr.contains(needle), "{extra:?}: {stderr}");
        for ext in ["props", "matches", "report", "report.zst"] {
            assert!(
                !Path::new(&format!("{o}.{ext}")).exists(),
                "{extra:?} {ext}"
            );
        }
    }
    // Overrides are validated before any work.
    let nan = run(&["analyze", "-r", r, "-o", o, "--rho", "nan", "--props"]);
    assert!(!nan.status.success());
    assert!(String::from_utf8_lossy(&nan.stderr).contains("--rho must be"));
    assert!(!Path::new(&format!("{o}.props")).exists());
    // --write-report without an EM parameter would only copy the input.
    let copy = run(&["analyze", "-r", r, "-o", o, "--write-report"]);
    assert!(!copy.status.success());
    assert!(String::from_utf8_lossy(&copy.stderr).contains("--write-report without EM parameters"));
    assert!(!Path::new(&format!("{o}.report")).exists());
    let missing = run(&[
        "analyze",
        "-r",
        "/nonexistent.report",
        "--matches",
        "--gamma",
        "0.6",
        "-o",
        o,
    ]);
    assert!(!missing.status.success());
    // gen_matches is gone.
    assert!(!run(&["gen_matches", "-r", r]).status.success());
}

/// `label → value` of a `.props` table.
fn props_map(table: &str) -> Vec<(String, f64)> {
    table
        .lines()
        .map(|l| {
            let (label, v) = l.split_once('\t').unwrap();
            (label.to_string(), v.parse::<f64>().unwrap())
        })
        .collect()
}

/// A PREMISE command line of a report, e.g. `premise_line(r, "EM")`.
fn premise_line(report: &str, cmd: &str) -> String {
    block(report, "PREMISE")
        .lines()
        .find_map(|l| l.trim().strip_prefix(&format!("{cmd} ")))
        .unwrap_or_else(|| panic!("no {cmd} line"))
        .trim_end_matches(';')
        .to_string()
}

/// The parsed `.posteriors` rows of a read: `(label, posterior)`, `None` for `-`.
fn posterior_rows(table: &str, read: &str) -> Vec<(String, Option<f64>)> {
    table
        .lines()
        .skip(1)
        .filter_map(|l| {
            let cols: Vec<&str> = l.split('\t').collect();
            (cols[0] == read).then(|| (cols[1].to_string(), cols[2].parse::<f64>().ok()))
        })
        .collect()
}

#[test]
fn analyze_writes_every_table_from_the_stored_fit() {
    let tmp = tempfile::tempdir().unwrap();
    let index = build_index(tmp.path());
    let cfg = Cfg {
        tree: true,
        rho: Some("3.5"),
    };
    let q = tmp.path().join("q");
    let mut args = query_args(&index, tmp.path(), &q, &cfg);
    args.extend(["--matches", "--gamma", "0.6"].map(str::to_string));
    ok(&run_owned(&args));
    let t = tmp.path().join("t");
    let out = run(&[
        "analyze",
        "-r",
        with_ext(&q, "report").to_str().unwrap(),
        "-o",
        t.to_str().unwrap(),
        "--matches",
        "--gamma",
        "0.6",
        "--props",
        "--posteriors",
        "--nodes",
    ]);
    ok(&out);
    let stdout = String::from_utf8_lossy(&out.stdout);
    for ext in ["matches", "posteriors", "props", "nodes"] {
        assert!(stdout.contains(&format!("t.{ext}")), "{stdout}");
        assert!(with_ext(&t, ext).exists(), "{ext}");
    }
    assert!(!with_ext(&t, "report").exists());
    let matches = read(&with_ext(&t, "matches"));
    assert_eq!(matches, read(&with_ext(&q, "matches")));
    let props = read(&with_ext(&t, "props"));
    assert!(props.starts_with("11\t1.0000000000e0"), "{props}");
    // .posteriors: one row per EM leaf of every classified read, summing to 1; the
    // unclassified reads of .matches are unclassified here too.
    let posteriors = read(&with_ext(&t, "posteriors"));
    assert_eq!(
        posteriors.lines().next().unwrap(),
        "ReadID\tTaxID\tPosterior"
    );
    let mut classified = 0;
    for line in matches.lines().skip(1) {
        let cols: Vec<&str> = line.split('\t').collect();
        let rows = posterior_rows(&posteriors, cols[0]);
        if cols[1] == "unclassified" {
            assert_eq!(rows, vec![("unclassified".to_string(), None)], "{line}");
            continue;
        }
        classified += 1;
        assert!(!rows.is_empty(), "{line}");
        let sum: f64 = rows.iter().map(|r| r.1.unwrap()).sum();
        assert!((sum - 1.0).abs() < 1e-4, "{line}: {rows:?}");
        // The fit collapsed everything into species 11, the only leaf.
        assert_eq!(rows[0].0, "11", "{rows:?}");
    }
    assert!(classified >= 1, "{matches}");
    assert_eq!(posteriors.lines().count(), 1 + 5);
    // Without the penalty both references stay leaves: reads aligning to both get two
    // rows in label order, and the pooled view disappears.
    let two = tmp.path().join("two");
    ok(&run(&[
        "analyze",
        "-r",
        with_ext(&q, "report").to_str().unwrap(),
        "-o",
        two.to_str().unwrap(),
        "--no-penalty",
        "--matches",
        "--gamma",
        "0.51",
        "--posteriors",
    ]));
    let matches = read(&with_ext(&two, "matches"));
    let posteriors = read(&with_ext(&two, "posteriors"));
    let mut seen: Vec<String> = Vec::new();
    for line in matches.lines().skip(1) {
        let cols: Vec<&str> = line.split('\t').collect();
        let rows = posterior_rows(&posteriors, cols[0]);
        if cols[1] == "unclassified" {
            assert_eq!(rows, vec![("unclassified".to_string(), None)], "{line}");
            continue;
        }
        let labels: Vec<&str> = rows.iter().map(|r| r.0.as_str()).collect();
        let mut sorted = labels.clone();
        sorted.sort_unstable();
        assert_eq!(labels, sorted, "{line}: {rows:?}");
        assert!(
            labels.iter().all(|l| ["ref_A", "ref_B"].contains(l)),
            "{rows:?}"
        );
        assert_eq!(labels[0], cols[1], "{line}: {rows:?}");
        let sum: f64 = rows.iter().map(|r| r.1.unwrap()).sum();
        assert!((sum - 1.0).abs() < 1e-4, "{line}: {rows:?}");
        seen.extend(rows.into_iter().map(|r| r.0));
    }
    // Both leaves appear (the fixture's reads each align to one reference only, so the
    // pooled multi-row case is covered by the unit tests).
    assert!(seen.iter().any(|l| l == "ref_A") && seen.iter().any(|l| l == "ref_B"));
}

#[test]
fn refit_with_the_reports_parameters_reproduces_the_fit() {
    let tmp = tempfile::tempdir().unwrap();
    let index = build_index(tmp.path());
    let cfg = Cfg {
        tree: true,
        rho: Some("3.5"),
    };
    let q = tmp.path().join("q");
    let mut args = query_args(&index, tmp.path(), &q, &cfg);
    args.extend(["--matches", "--gamma", "0.6"].map(str::to_string));
    ok(&run_owned(&args));
    let q_report = with_ext(&q, "report");
    let base = tmp.path().join("base");
    ok(&run(&[
        "analyze",
        "-r",
        q_report.to_str().unwrap(),
        "-o",
        base.to_str().unwrap(),
        "--matches",
        "--gamma",
        "0.6",
        "--props",
        "--nodes",
        "--posteriors",
    ]));
    // The same EM parameters again, explicitly: a refit.
    let same = tmp.path().join("same");
    let out = run(&[
        "analyze",
        "-r",
        q_report.to_str().unwrap(),
        "-o",
        same.to_str().unwrap(),
        "--rho",
        "3.5",
        "--matches",
        "--gamma",
        "0.6",
        "--props",
        "--nodes",
        "--posteriors",
        "--write-report",
    ]);
    ok(&out);
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(stdout.contains("Refit of q.report"), "{stdout}");
    assert!(stdout.contains("same.report"), "{stdout}");
    // Proportions agree up to float summation order (the rows are in report order now).
    let want = props_map(&read(&with_ext(&base, "props")));
    let got = props_map(&read(&with_ext(&same, "props")));
    assert_eq!(want.len(), got.len());
    for ((la, a), (lb, b)) in want.iter().zip(&got) {
        assert_eq!(la, lb);
        assert!(
            (a - b).abs() <= 1e-9 * a.abs().max(1e-300),
            "{la}: {a} vs {b}"
        );
    }
    // Same vertices, parents, ranks and contig counts in .nodes.
    let first4 = |t: &str| -> Vec<String> {
        t.lines()
            .map(|l| l.split('\t').take(4).collect::<Vec<_>>().join("\t"))
            .collect()
    };
    assert_eq!(
        first4(&read(&with_ext(&same, "nodes"))),
        first4(&read(&with_ext(&base, "nodes")))
    );
    // Same assignments and posteriors (to the printed digits) as the stored fit.
    let two_cols = |t: &str| -> Vec<String> {
        t.lines()
            .map(|l| l.split('\t').take(2).collect::<Vec<_>>().join("\t"))
            .collect()
    };
    assert_eq!(
        two_cols(&read(&with_ext(&same, "matches"))),
        two_cols(&read(&with_ext(&base, "matches")))
    );
    assert_eq!(
        two_cols(&read(&with_ext(&same, "posteriors"))),
        two_cols(&read(&with_ext(&base, "posteriors")))
    );
    // The refit report: lineage, the same EM, the alignment rows copied through.
    let q_text = read(&q_report);
    let r_text = read(&with_ext(&same, "report"));
    assert_eq!(premise_line(&r_text, "SOURCE"), "q.report");
    assert_eq!(premise_line(&r_text, "VERSION"), "3");
    assert_eq!(
        premise_line(&r_text, "PARAMS"),
        premise_line(&q_text, "PARAMS")
    );
    assert_eq!(
        premise_line(&r_text, "READS"),
        premise_line(&q_text, "READS")
    );
    assert_eq!(
        premise_line(&r_text, "TAXONOMY"),
        premise_line(&q_text, "TAXONOMY")
    );
    assert_eq!(
        premise_line(&r_text, "INDEX"),
        premise_line(&q_text, "INDEX")
    );
    assert_eq!(
        premise_line(&r_text, "INDEX_SHA256"),
        premise_line(&q_text, "INDEX_SHA256")
    );
    let em = |t: &str| -> String {
        premise_line(t, "EM")
            .split_whitespace()
            .filter(|kv| !kv.starts_with("final_ll="))
            .collect::<Vec<_>>()
            .join(" ")
    };
    assert_eq!(em(&r_text), em(&q_text));
    assert!(premise_line(&r_text, "RUNTIME").contains("align_s=0e0"));
    assert_eq!(
        block(&r_text, "ALIGNMENT_LIKELIHOODS"),
        block(&q_text, "ALIGNMENT_LIKELIHOODS")
    );
    assert_eq!(block(&r_text, "REFMAP"), block(&q_text, "REFMAP"));
    assert!(!premise_line(&q_text, "NODES").is_empty());
    assert!(!q_text.contains("SOURCE"));
    // Analysing the refit report reproduces the refit's own tables byte for byte.
    let chain = tmp.path().join("chain");
    ok(&run(&[
        "analyze",
        "-r",
        with_ext(&same, "report").to_str().unwrap(),
        "-o",
        chain.to_str().unwrap(),
        "--matches",
        "--gamma",
        "0.6",
        "--props",
        "--posteriors",
    ]));
    for ext in ["matches", "posteriors", "props"] {
        assert_eq!(
            read(&with_ext(&chain, ext)),
            read(&with_ext(&same, ext)),
            "{ext}"
        );
    }
}

#[test]
fn refit_overrides_change_the_penalty_and_are_recorded() {
    let tmp = tempfile::tempdir().unwrap();
    let index = build_index(tmp.path());
    let cfg = Cfg {
        tree: true,
        rho: Some("3.5"),
    };
    let q = tmp.path().join("q");
    ok(&run_owned(&query_args(&index, tmp.path(), &q, &cfg)));
    let q_report = with_ext(&q, "report");
    let r = q_report.to_str().unwrap();
    let q_text = read(&q_report);
    assert!(
        premise_line(&q_text, "EM").contains("collapsed=1"),
        "{q_text}"
    );
    let refit = |tag: &str, extra: &[&str]| -> String {
        let out = tmp.path().join(tag);
        let mut args = vec!["analyze", "-r", r, "-o", out.to_str().unwrap(), "--props"];
        args.extend_from_slice(extra);
        ok(&run(&args));
        assert!(!with_ext(&out, "matches").exists());
        read(&with_ext(&out, "props"))
    };
    let report_of = |tag: &str, extra: &[&str]| -> String {
        let out = tmp.path().join(tag);
        let mut args = vec![
            "analyze",
            "-r",
            r,
            "-o",
            out.to_str().unwrap(),
            "--write-report",
        ];
        args.extend_from_slice(extra);
        ok(&run(&args));
        read(&with_ext(&out, "report"))
    };
    // No penalty: nothing pruned or collapsed, both sequences keep their own leaf.
    let text = report_of("nop", &["--no-penalty"]);
    assert!(
        premise_line(&text, "EM").contains("pruned=0 collapsed=0"),
        "{text}"
    );
    assert!(premise_line(&text, "PARAMS").ends_with("init=split penalty=0"));
    assert!(premise_line(&text, "PARAMS").contains("rho=3.5e0"));
    let props = refit("nop-props", &["--no-penalty"]);
    assert_eq!(props.lines().count(), 2, "{props}");
    assert!(
        props.contains("ref_A\t") && props.contains("ref_B\t"),
        "{props}"
    );
    // A huge rho collapses the whole sample into the root.
    let text = report_of("big", &["--rho", "1000"]);
    assert!(
        premise_line(&text, "EM").contains("root_collapsed=1"),
        "{text}"
    );
    assert!(premise_line(&text, "PARAMS").contains("rho=1e3"));
    // Every override is recorded.
    let text = report_of(
        "all",
        &[
            "--rho",
            "2",
            "--omega",
            "1e-5",
            "--iter",
            "1",
            "--em_threshold",
            "1e-3",
            "--init",
            "uniform",
        ],
    );
    let params = premise_line(&text, "PARAMS");
    for needle in [
        "iter=1 ",
        "rho=2e0",
        "omega=1e-5",
        "em_threshold=1e-3",
        "init=uniform",
        "penalty=1",
        "mem=5 ",
    ] {
        assert!(params.contains(needle), "{needle} in {params}");
    }
    assert!(
        premise_line(&text, "EM").starts_with("iterations=1 "),
        "{text}"
    );
    // A compressed refit report is read back like any other.
    let z = tmp.path().join("z");
    ok(&run(&[
        "analyze",
        "-r",
        r,
        "-o",
        z.to_str().unwrap(),
        "--no-penalty",
        "--write-report",
        "--compress",
        "--matches",
        "--gamma",
        "0.9",
    ]));
    let zpath = tmp.path().join("z.report.zst");
    assert_eq!(&fs::read(&zpath).unwrap()[..4], &[0x28, 0xB5, 0x2F, 0xFD]);
    assert!(!with_ext(&z, "report").exists());
    let again = tmp.path().join("again");
    ok(&run(&[
        "analyze",
        "-r",
        zpath.to_str().unwrap(),
        "-o",
        again.to_str().unwrap(),
        "--matches",
        "--gamma",
        "0.9",
    ]));
    assert_eq!(
        read(&with_ext(&again, "matches")),
        read(&with_ext(&z, "matches"))
    );
    // A corrupt report leaves no refit outputs behind.
    let short = q_text.replacen("DIMENSIONS NREADS=5;", "DIMENSIONS NREADS=6;", 1);
    let bad = tmp.path().join("bad.report");
    fs::write(&bad, short).unwrap();
    let b = tmp.path().join("b");
    let out = run(&[
        "analyze",
        "-r",
        bad.to_str().unwrap(),
        "-o",
        b.to_str().unwrap(),
        "--no-penalty",
        "--props",
        "--write-report",
    ]);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("announces 6 reads but holds 5"));
    assert!(!with_ext(&b, "props").exists() && !with_ext(&b, "report").exists());
    // ... but a file of an earlier run that this run never opened survives.
    fs::write(with_ext(&b, "posteriors"), "old\n").unwrap();
    let out = run(&[
        "analyze",
        "-r",
        bad.to_str().unwrap(),
        "-o",
        b.to_str().unwrap(),
        "--matches",
        "--gamma",
        "0.6",
        "--props",
    ]);
    assert!(!out.status.success());
    assert_eq!(read(&with_ext(&b, "posteriors")), "old\n");
    assert!(!with_ext(&b, "props").exists() && !with_ext(&b, "matches").exists());
}
