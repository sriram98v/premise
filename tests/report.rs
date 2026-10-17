//! The report: `analyze` reproduces `query`'s `.matches` at any gamma, with and without
//! compression; `query` writes what its flags ask for; the head records the run without
//! paths.

use sha2::{Digest, Sha256};
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

/// The rows of a `.matches` table after its leading `#` EM summary: the column header,
/// then one row per read.
fn table_lines(table: &str) -> impl Iterator<Item = &str> {
    table.lines().filter(|l| !l.starts_with('#'))
}

/// The leading `# …` EM summary line of a `.matches` table, checked for its shape and
/// for the gamma the table was cut at.
fn matches_summary(table: &str, gamma: &str) -> String {
    let first = table.lines().next().expect("a summary line");
    assert!(first.starts_with("# PREMISE "), "{first}");
    assert!(first.contains(" EM: init="), "{first}");
    assert!(first.ends_with(&format!("; gamma={gamma}")), "{first}");
    assert!(
        !table.lines().skip(1).any(|l| l.starts_with('#')),
        "a single summary line\n{table}"
    );
    first.to_string()
}

/// The `[…]` EM summary on the second comment line of a NEXUS `.props`, without brackets.
fn props_summary(props: &str) -> String {
    let line = props.lines().nth(2).expect("a summary line");
    assert!(
        props.starts_with("#NEXUS\n[PREMISE proportions v1]\n["),
        "{props}"
    );
    let inner = line
        .strip_prefix('[')
        .and_then(|l| l.strip_suffix(']'))
        .unwrap_or_else(|| panic!("{line}"));
    assert!(inner.contains("EM: init="), "{line}");
    inner.to_string()
}

/// The `.matches` summary of a table cut at `gamma` and the `.props` summary of the same
/// fit: the same text but for the `; gamma=…` that only the table has.
fn assert_same_summary(matches: &str, props: &str, gamma: &str) {
    let m = matches_summary(matches, gamma);
    let want = m
        .strip_prefix("# ")
        .and_then(|s| s.strip_suffix(&format!("; gamma={gamma}")))
        .unwrap();
    assert_eq!(props_summary(props), want);
}

/// `key=value` of an EM summary or a PARAMS line.
fn kv(text: &str, key: &str) -> String {
    text.split(|c: char| c.is_whitespace() || c == ';' || c == ',')
        .find_map(|t| t.strip_prefix(&format!("{key}=")))
        .unwrap_or_else(|| panic!("no {key}= in {text}"))
        .to_string()
}

/// The EM settings of a summary are the report's PARAMS.
fn assert_summary_params(summary: &str, report: &str) {
    let params = premise_line(report, "PARAMS");
    for key in ["init", "iter", "rho", "omega", "em_threshold", "penalty"] {
        assert_eq!(kv(summary, key), kv(&params, key), "{summary}\nvs {params}");
    }
}

/// `.matches` without its position columns: tied placements land on any of the
/// equally good positions, so two runs only agree up to them. Without the EM summary
/// too: its full-precision log-likelihood is a float sum whose order two runs need not
/// share.
fn no_pos(table: &str) -> String {
    table_lines(table)
        .map(|l| l.split('\t').take(3).collect::<Vec<_>>().join("\t"))
        .collect::<Vec<_>>()
        .join("\n")
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

/// One block of a report, from `BEGIN <name>;` to its `END;`.
fn block(report: &str, name: &str) -> String {
    let start = report
        .find(&format!("BEGIN {name};"))
        .unwrap_or_else(|| panic!("no {name} block"));
    let end = report[start..].find("END;").expect("block end") + start;
    report[start..end].to_string()
}

/// The MATRIX lines of a block (trimmed), up to its closing `;`.
fn matrix_lines(block: &str) -> Vec<String> {
    block
        .lines()
        .skip_while(|l| !l.contains("MATRIX"))
        .skip(1)
        .take_while(|l| l.trim() != ";")
        .map(|l| l.trim().to_string())
        .collect()
}

/// A delta-coded seqid line (first absolute, then differences) as absolute seqids.
fn seqids(line: &str) -> Vec<usize> {
    line.split_whitespace()
        .scan(0usize, |acc, t| {
            *acc += t.parse::<usize>().unwrap();
            Some(*acc)
        })
        .collect()
}

/// Every read of a report's ALIGNMENT_LIKELIHOODS in order: (read id, [(seqid,
/// log-likelihood)] in row order). Checks the announced read count on the way.
fn alignment_rows(report: &str) -> Vec<(String, Vec<(usize, f64)>)> {
    let blk = block(report, "ALIGNMENT_LIKELIHOODS");
    let lines = matrix_lines(&blk);
    let mut out = Vec::new();
    let mut i = 0;
    while i < lines.len() {
        let hdr = nexus_tokens(&lines[i]);
        i += 1;
        if hdr[1] == "0" {
            out.push((hdr[0].clone(), Vec::new()));
            continue;
        }
        let n: usize = hdr[1].parse().unwrap();
        let k: usize = hdr[2].parse().unwrap();
        let seqs = seqids(&lines[i]);
        let lls: Vec<f64> = lines[i + 1]
            .split_whitespace()
            .map(|t| t.parse().unwrap())
            .collect();
        let idx: Vec<usize> = if k > 1 {
            lines[i + 2]
                .split_whitespace()
                .map(|t| t.parse().unwrap())
                .collect()
        } else {
            vec![0; n]
        };
        assert_eq!((seqs.len(), lls.len(), idx.len()), (n, k, n), "{hdr:?}");
        out.push((
            hdr[0].clone(),
            seqs.into_iter()
                .zip(idx)
                .map(|(s, j)| (s, lls[j]))
                .collect(),
        ));
        // seqids, log-likelihoods, indices when k > 1, positions, insert offsets
        i += if k > 1 { 5 } else { 4 };
    }
    assert!(
        blk.contains(&format!("DIMENSIONS NREADS={};", out.len())),
        "{blk}"
    );
    out
}

/// Every read of a report's ALIGNMENT_LIKELIHOODS in order: (read id, aligned seqids).
fn alignment_seqids(report: &str) -> Vec<(String, Vec<usize>)> {
    alignment_rows(report)
        .into_iter()
        .map(|(id, rows)| (id, rows.into_iter().map(|r| r.0).collect()))
        .collect()
}

/// Every read's posterior on each of its alignment rows, recomputed from what the report
/// stores (it keeps no posteriors): `π̃_s · exp(ll_s − max ll) / Σ`, where `π̃_s` is the
/// NODES proportion of the final EM leaf owning the sequence's tip shared among its
/// contigs (0 once pruned). All zeros when none of the read's sequences has mass.
fn posteriors(report: &str) -> Vec<(String, Vec<(usize, f64)>)> {
    let tree = Tree::of(report);
    alignment_rows(report)
        .into_iter()
        .map(|(id, rows)| {
            let m = rows.iter().map(|r| r.1).fold(f64::NEG_INFINITY, f64::max);
            let w: Vec<f64> = rows
                .iter()
                .map(|&(s, ll)| tree.pi_tilde(s) * (ll - m).exp())
                .collect();
            let d: f64 = w.iter().sum();
            let post = rows
                .iter()
                .zip(&w)
                .map(|(r, w)| (r.0, if d > 0.0 { w / d } else { 0.0 }))
                .collect();
            (id, post)
        })
        .collect()
}

/// The tree of a report: TAXLABELS, NODES and REFMAP joined up.
struct Tree {
    /// label → NODES id
    id: HashMap<String, String>,
    /// NODES id → (parent id or `-`, status)
    node: HashMap<String, (String, String)>,
    /// NODES id → (contigs, proportion)
    fit: HashMap<String, (usize, f64)>,
    /// seqid → NODES id of the sequence's own vertex
    seq: HashMap<usize, String>,
}

impl Tree {
    fn of(report: &str) -> Tree {
        let labels: Vec<String> = block(report, "TAXA")
            .lines()
            .skip_while(|l| !l.contains("TAXLABELS"))
            .skip(1)
            .take_while(|l| l.trim() != ";")
            .flat_map(nexus_tokens)
            .collect();
        let rows: Vec<Vec<String>> = matrix_lines(&block(report, "NODES"))
            .iter()
            .map(|l| nexus_tokens(l))
            .collect();
        assert_eq!(labels.len(), rows.len());
        let seq = matrix_lines(&block(report, "REFMAP"))
            .iter()
            .map(|l| {
                let (s, n) = l.split_once(' ').unwrap();
                (s.parse().unwrap(), n.to_string())
            })
            .collect();
        Tree {
            id: labels
                .into_iter()
                .zip(&rows)
                .map(|(l, r)| (l, r[0].clone()))
                .collect(),
            node: rows
                .iter()
                .map(|r| (r[0].clone(), (r[1].clone(), r[3].clone())))
                .collect(),
            fit: rows
                .iter()
                .map(|r| (r[0].clone(), (r[4].parse().unwrap(), r[6].parse().unwrap())))
                .collect(),
            seq,
        }
    }

    /// The final EM leaf owning vertex `v`: `v` itself unless it was pooled (then its
    /// parent's owner, i.e. the nearest collapsed ancestor) or pruned (`None`).
    fn em_leaf<'a>(&'a self, v: &'a str) -> Option<&'a str> {
        let (parent, status) = &self.node[v];
        match status.as_str() {
            "pruned" => None,
            "pooled" => self.em_leaf(parent),
            _ => Some(v),
        }
    }

    /// `π̃_s`: the proportion of the EM leaf owning sequence `s`, over its contigs.
    fn pi_tilde(&self, s: usize) -> f64 {
        match self.em_leaf(&self.seq[&s]).map(|l| self.fit[l]) {
            Some((n, p)) if n > 0 => p / n as f64,
            _ => 0.0,
        }
    }

    /// Whether vertex `anc` is `v` or one of its ancestors.
    fn covers(&self, anc: &str, v: &str) -> bool {
        let mut cur = v;
        loop {
            if cur == anc {
                return true;
            }
            match self.node.get(cur) {
                Some((p, _)) if p != "-" => cur = p,
                _ => return false,
            }
        }
    }

    /// A read's posterior on the subtree under vertex `v` (a NODES id).
    fn mass(&self, v: &str, posts: &[(usize, f64)]) -> f64 {
        posts
            .iter()
            .filter(|(s, _)| self.covers(v, &self.seq[s]))
            .map(|(_, p)| p)
            .sum()
    }

    /// The final-tree (`active` / `collapsed`) children of vertex `v`.
    fn final_children(&self, v: &str) -> Vec<&str> {
        self.node
            .iter()
            .filter(|(_, (p, st))| p == v && (st == "active" || st == "collapsed"))
            .map(|(id, _)| id.as_str())
            .collect()
    }

    fn root(&self) -> &str {
        self.node
            .iter()
            .find(|(_, (p, _))| p == "-")
            .map(|(id, _)| id.as_str())
            .expect("a root")
    }
}

/// The report holds no PROPS or READ_POSTERIORS (the EM regenerates both from NODES and
/// ALIGNMENT_LIKELIHOODS) and its blocks come in their fixed order.
fn check_blocks(report: &str) {
    assert!(!report.contains("BEGIN PROPS;"), "{report}");
    assert!(!report.contains("READ_POSTERIORS"), "{report}");
    let order: Vec<&str> = report
        .lines()
        .filter_map(|l| l.strip_prefix("BEGIN "))
        .map(|l| l.trim_end_matches(';'))
        .collect();
    assert_eq!(
        order,
        [
            "PREMISE",
            "EM_TRACE",
            "TAXA",
            "TREES",
            "NODES",
            "REFMAP",
            "ALIGNMENT_LIKELIHOODS",
            "PROFILE"
        ]
    );
}

/// Every row of `matches` is the gamma descent over the posteriors recomputed from the
/// report ([`posteriors`]): the assigned vertex holds the printed posterior, which is
/// above gamma, and none of its final children is; an unclassified read has no child of
/// the root above gamma. Each read's posteriors sum to one, or to zero when its
/// sequences were all pruned.
fn check_posteriors(report: &str, matches: &str, gamma: f64) {
    check_blocks(report);
    let posts = posteriors(report);
    for (id, rows) in &posts {
        let sum: f64 = rows.iter().map(|r| r.1).sum();
        assert!((sum - 1.0).abs() < 1e-9 || sum == 0.0, "{id}: {rows:?}");
    }
    let tree = Tree::of(report);
    let by_read: HashMap<&str, &[(usize, f64)]> = posts
        .iter()
        .map(|(id, rows)| (id.as_str(), rows.as_slice()))
        .collect();
    let mut n = 0;
    for line in table_lines(matches).skip(1) {
        n += 1;
        let cols: Vec<&str> = line.split('\t').collect();
        let rows = by_read[cols[0]];
        let below = |v: &str| -> Vec<f64> {
            tree.final_children(v)
                .into_iter()
                .map(|c| tree.mass(c, rows))
                .collect()
        };
        if cols[1] == "unclassified" {
            let kids = below(tree.root());
            assert!(kids.iter().all(|m| *m <= gamma), "{line}: {kids:?}");
            continue;
        }
        let v = &tree.id[cols[1]];
        let mass = tree.mass(v, rows);
        let printed: f64 = cols[2].parse().unwrap();
        assert!((mass - printed).abs() < 1e-4, "{line}: {mass}");
        assert!(mass > gamma, "{line}: {mass}");
        let kids = below(v);
        assert!(kids.iter().all(|m| *m <= gamma), "{line}: {kids:?}");
    }
    assert_eq!(n, posts.len(), "one .matches row per read");
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
            assert!(
                want.lines().count() == 7,
                "summary + header + 5 reads: {want}"
            );
            // The summary is the run's EM, cut at this gamma, as in its .props.
            let direct_report = read(&with_ext(&direct, "report"));
            let summary = matches_summary(&want, g);
            assert_summary_params(&summary, &direct_report);
            assert_eq!(kv(&summary, "iter"), "20", "{summary}");
            let rho = cfg.rho.map_or("0e0".to_string(), |r| format!("{r}e0"));
            assert_eq!(kv(&summary, "rho"), rho, "{summary}");
            assert!(!summary.contains("refit of"), "{summary}");
            assert_same_summary(&want, &read(&with_ext(&direct, "props")), g);
            // Every row is the gamma descent over the posteriors the report implies.
            check_posteriors(&direct_report, &want, g.parse().unwrap());

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
    // By default: the report and the proportions.
    let dflt = tmp.path().join("dflt");
    let out = run_owned(&query_args(&index, tmp.path(), &dflt, &cfg));
    ok(&out);
    assert!(with_ext(&dflt, "report").exists() && with_ext(&dflt, "props").exists());
    assert!(!with_ext(&dflt, "matches").exists());
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        stdout.contains("dflt.report") && stdout.contains("dflt.props"),
        "{stdout}"
    );
    // --no-report with --matches: the table and the proportions, no report.
    let only = tmp.path().join("only");
    let mut q = query_args(&index, tmp.path(), &only, &cfg);
    q.extend(["--matches", "--gamma", "0.6", "--no-report"].map(str::to_string));
    let out = run_owned(&q);
    ok(&out);
    assert!(with_ext(&only, "matches").exists() && with_ext(&only, "props").exists());
    assert!(!with_ext(&only, "report").exists());
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        stdout.contains("only.matches")
            && stdout.contains("only.props")
            && !stdout.contains("report"),
        "{stdout}"
    );
    // --no-report alone is fine too: only the proportions.
    let props_only = tmp.path().join("po");
    let mut q = query_args(&index, tmp.path(), &props_only, &cfg);
    q.push("--no-report".to_string());
    let out = run_owned(&q);
    ok(&out);
    assert!(with_ext(&props_only, "props").exists());
    assert!(
        !with_ext(&props_only, "matches").exists() && !with_ext(&props_only, "report").exists()
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        stdout.contains("po.props") && !stdout.contains("report"),
        "{stdout}"
    );
    // The proportions do not depend on which other outputs were asked for.
    let leaves = props_map(&read(&with_ext(&props_only, "props")));
    assert_eq!(
        leaves.iter().map(|l| l.0.as_str()).collect::<Vec<_>>(),
        ["ref_A", "ref_B"]
    );
    for other in [&dflt, &only] {
        let theirs = props_map(&read(&with_ext(other, "props")));
        assert_eq!(theirs.len(), leaves.len());
        for ((la, a), (lb, b)) in leaves.iter().zip(&theirs) {
            assert_eq!(la, lb);
            assert!((a - b).abs() <= 1e-9, "{la}: {a} vs {b}");
        }
    }
    // --matches without --gamma uses gamma 0.6.
    let dg = tmp.path().join("dg");
    let mut q = query_args(&index, tmp.path(), &dg, &cfg);
    q.push("--matches".to_string());
    ok(&run_owned(&q));
    let at6 = tmp.path().join("at6");
    ok(&run(&[
        "analyze",
        "-r",
        with_ext(&dg, "report").to_str().unwrap(),
        "--matches",
        "--gamma",
        "0.6",
        "-o",
        at6.to_str().unwrap(),
    ]));
    assert_eq!(
        read(&with_ext(&dg, "matches")),
        read(&with_ext(&at6, "matches"))
    );
    let bad = tmp.path().join("bad");
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
    ]));
    // Head-only export needs no --matches.
    ok(&run(&[
        "analyze",
        "-r",
        tmp.path().join("plain.report").to_str().unwrap(),
        "-o",
        tmp.path().join("p").to_str().unwrap(),
        "--props",
    ]));
    assert!(!tmp.path().join("p.matches").exists());
    assert_eq!(
        no_pos(&read(&with_ext(&from_zst, "matches"))),
        no_pos(&read(&with_ext(&plain, "matches")))
    );
    // `analyze --props` reproduces each query's own .props byte for byte; the two runs
    // agree up to float summation order.
    assert_eq!(
        read(&tmp.path().join("z.props")),
        read(&tmp.path().join("zst.props"))
    );
    assert_eq!(
        read(&tmp.path().join("p.props")),
        read(&tmp.path().join("plain.props"))
    );
    assert_same_props(
        &read(&tmp.path().join("z.props")),
        &read(&tmp.path().join("p.props")),
    );
    // Both references collapsed into species 11, which pools them.
    let nodes = final_nodes(&read(&tmp.path().join("plain.report")));
    assert!(
        nodes.contains(&(
            "11".to_string(),
            "10".to_string(),
            "species".to_string(),
            "2".to_string()
        )),
        "{nodes:?}"
    );
    let props = read(&tmp.path().join("p.props"));
    let rows: Vec<(String, String)> = props_rows(&props).into_iter().map(|r| (r.0, r.1)).collect();
    let pair = |a: &str, b: &str| (a.to_string(), b.to_string());
    assert_eq!(
        rows,
        vec![
            pair("root", "-"),
            pair("1", "root"),
            pair("10", "1"),
            pair("11", "10")
        ],
        "{props}"
    );
    let leaves = props_map(&props);
    assert_eq!(leaves.len(), 1, "{props}");
    assert_eq!(leaves[0].0, "11", "{props}");
    assert!((leaves[0].1 - 1.0).abs() < 1e-12, "{props}");
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
        "#NEXUS\n[PREMISE report v4]\n",
        "BEGIN PREMISE;",
        "\tVERSION 4;",
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
        "BEGIN ALIGNMENT_LIKELIHOODS;",
        "BEGIN PROFILE;",
        "peak_rss_kb=",
    ] {
        assert!(report.contains(needle), "missing {needle:?} in\n{report}");
    }
    // No PROPS or READ_POSTERIORS (the final leaves' proportions are NODES', the
    // posteriors follow from them and ALIGNMENT_LIKELIHOODS), the blocks in order; the
    // recomputed posteriors explain .matches.
    check_blocks(&report);
    let matches = read(&tmp.path().join("meta.matches"));
    check_posteriors(&report, &matches, 0.6);
    // The .matches summary is this run's EM: the PARAMS, then the EM line's iterations
    // and final log-likelihood.
    let summary = matches_summary(&matches, "0.6");
    assert_summary_params(&summary, &report);
    assert!(
        summary.contains("rho=3.5e0") && summary.contains("penalty=1"),
        "{summary}"
    );
    let em = premise_line(&report, "EM");
    assert!(
        summary.contains(&format!("; {} iterations, ", kv(&em, "iterations"))),
        "{summary} vs {em}"
    );
    let final_ll = kv(&em, "final_ll");
    assert!(
        summary.contains(&format!("final penalized log-likelihood {final_ll}")),
        "{summary} vs {em}"
    );
    assert_same_summary(&matches, &read(&tmp.path().join("meta.props")), "0.6");
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
    assert!(String::from_utf8_lossy(&no_output.stderr)
        .contains("analyze: nothing to do; give --matches, --props or --write-report"));
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

/// The PROPORTIONS rows of a NEXUS `.props` in pre-order: (label, parent label or `-`,
/// rank, proportion). Checks that the first row is the root with all the mass and that
/// every other row's parent is an earlier row.
fn props_rows(text: &str) -> Vec<(String, String, String, f64)> {
    assert!(
        text.starts_with("#NEXUS\n[PREMISE proportions v1]\n"),
        "{text}"
    );
    let rows: Vec<(String, String, String, f64)> = matrix_lines(&block(text, "PROPORTIONS"))
        .iter()
        .map(|l| {
            let c = nexus_tokens(l);
            assert_eq!(c.len(), 4, "{l}");
            (
                c[0].clone(),
                c[1].clone(),
                c[2].clone(),
                c[3].parse().unwrap(),
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

/// `label → value` of the leaves of a `.props` tree (rows that are no other row's
/// parent), in row order: the final EM leaves with mass.
fn props_map(text: &str) -> Vec<(String, f64)> {
    let rows = props_rows(text);
    rows.iter()
        .filter(|(label, ..)| !rows.iter().any(|r| r.1 == *label))
        .map(|(label, _, _, p)| (label.clone(), *p))
        .collect()
}

/// Two `.props` trees agree row by row on labels, parents and ranks, and on the
/// proportions up to float summation order.
fn assert_same_props(a: &str, b: &str) {
    let (ra, rb) = (props_rows(a), props_rows(b));
    assert_eq!(ra.len(), rb.len(), "{a}\nvs\n{b}");
    for (x, y) in ra.iter().zip(&rb) {
        assert_eq!((&x.0, &x.1, &x.2), (&y.0, &y.1, &y.2));
        assert!(
            (x.3 - y.3).abs() <= 1e-9 * x.3.abs().max(1e-300),
            "{}: {} vs {}",
            x.0,
            x.3,
            y.3
        );
    }
}

/// The final-tree rows of a report's NODES block (`active` / `collapsed`), as (label,
/// parent label or `-`, rank, contigs).
fn final_nodes(report: &str) -> Vec<(String, String, String, String)> {
    let labels: Vec<String> = block(report, "TAXA")
        .lines()
        .skip_while(|l| !l.contains("TAXLABELS"))
        .skip(1)
        .take_while(|l| l.trim() != ";")
        .flat_map(nexus_tokens)
        .collect();
    let rows: Vec<Vec<String>> = matrix_lines(&block(report, "NODES"))
        .iter()
        .map(|l| nexus_tokens(l))
        .collect();
    let label_of: HashMap<&str, &str> = rows
        .iter()
        .zip(&labels)
        .map(|(r, l)| (r[0].as_str(), l.as_str()))
        .collect();
    rows.iter()
        .zip(&labels)
        .filter(|(r, _)| r[3] == "active" || r[3] == "collapsed")
        .map(|(r, l)| {
            let parent = if r[1] == "-" {
                "-"
            } else {
                label_of[r[1].as_str()]
            };
            (l.clone(), parent.to_string(), r[2].clone(), r[4].clone())
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
    ]);
    ok(&out);
    let stdout = String::from_utf8_lossy(&out.stdout);
    for ext in ["matches", "props"] {
        assert!(stdout.contains(&format!("t.{ext}")), "{stdout}");
        assert!(with_ext(&t, ext).exists(), "{ext}");
    }
    for ext in ["report", "posteriors", "nodes"] {
        assert!(!with_ext(&t, ext).exists(), "{ext}");
    }
    let matches = read(&with_ext(&t, "matches"));
    assert_eq!(matches, read(&with_ext(&q, "matches")));
    // .props is the query's own, the fit collapsed everything into species 11.
    let props = read(&with_ext(&t, "props"));
    assert_eq!(props, read(&with_ext(&q, "props")));
    let leaves = props_map(&props);
    assert_eq!(leaves.len(), 1, "{props}");
    assert_eq!(leaves[0].0, "11", "{props}");
    assert!((leaves[0].1 - 1.0).abs() < 1e-12, "{props}");
    // The posteriors the report implies: one entry per read, the classified reads of
    // .matches at species 11 with their posterior there.
    let report = read(&with_ext(&q, "report"));
    check_posteriors(&report, &matches, 0.6);
    let posts = posteriors(&report);
    assert_eq!(posts.len(), 5);
    let classified: Vec<&str> = table_lines(&matches)
        .skip(1)
        .filter_map(|l| l.split('\t').nth(1))
        .filter(|l| *l != "unclassified")
        .collect();
    assert!(!classified.is_empty(), "{matches}");
    assert!(classified.iter().all(|l| *l == "11"), "{matches}");
    // The flags of the old tables are gone.
    for flag in ["--nodes", "--posteriors"] {
        let gone = run(&[
            "analyze",
            "-r",
            with_ext(&q, "report").to_str().unwrap(),
            "-o",
            tmp.path().join("gone").to_str().unwrap(),
            flag,
        ]);
        assert!(!gone.status.success(), "{flag}");
    }
    // Without the penalty both references stay leaves: the refit report's posteriors put
    // each read on its own reference.
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
        "--write-report",
    ]));
    let matches = read(&with_ext(&two, "matches"));
    let report = read(&with_ext(&two, "report"));
    check_posteriors(&report, &matches, 0.51);
    // The refit's summary is its own EM (penalty off), naming the report it refits.
    let summary = matches_summary(&matches, "0.51");
    assert_summary_params(&summary, &report);
    assert!(summary.contains("penalty=0"), "{summary}");
    assert!(
        summary.contains(", refit of q.report; gamma=0.51"),
        "{summary}"
    );
    let seen: Vec<&str> = table_lines(&matches)
        .skip(1)
        .filter_map(|l| l.split('\t').nth(1))
        .filter(|l| *l != "unclassified")
        .collect();
    assert!(
        seen.iter().all(|l| ["ref_A", "ref_B"].contains(l)),
        "{matches}"
    );
    // Both leaves appear (the fixture's reads each align to one reference only, so the
    // pooled multi-row case is covered by the unit tests).
    assert!(
        seen.contains(&"ref_A") && seen.contains(&"ref_B"),
        "{matches}"
    );
}

/// The report keeps no PROPS or READ_POSTERIORS, in a query's report or a refit's: the
/// posteriors are regenerated from ALIGNMENT_LIKELIHOODS and the fit's NODES, and each
/// report's regenerated posteriors are its own fit's and explain its own `.matches`.
#[test]
fn report_omits_what_the_em_regenerates() {
    let tmp = tempfile::tempdir().unwrap();
    let index = build_index(tmp.path());
    let cfg = Cfg {
        tree: true,
        rho: None,
    };
    let q = tmp.path().join("q");
    let mut args = query_args(&index, tmp.path(), &q, &cfg);
    args.extend(["--matches", "--gamma", "0.6"].map(str::to_string));
    ok(&run_owned(&args));
    let q_text = read(&with_ext(&q, "report"));
    check_blocks(&q_text);
    let aligned: Vec<String> = alignment_seqids(&q_text).into_iter().map(|r| r.0).collect();
    assert_eq!(aligned.len(), 5);
    check_posteriors(&q_text, &read(&with_ext(&q, "matches")), 0.6);
    // Unpenalized, both species stay: every read is all on the reference it aligns to.
    for (_, rows) in posteriors(&q_text) {
        assert_eq!(rows.len(), 1, "{rows:?}");
        assert!((rows[0].1 - 1.0).abs() < 1e-12, "{rows:?}");
    }
    // A refit with the penalty writes a v4 report with the same blocks, no more.
    let r = tmp.path().join("r");
    let out = run(&[
        "analyze",
        "-r",
        with_ext(&q, "report").to_str().unwrap(),
        "-o",
        r.to_str().unwrap(),
        "--rho",
        "3.5",
        "--write-report",
        "--matches",
        "--gamma",
        "0.6",
    ]);
    ok(&out);
    let r_text = read(&with_ext(&r, "report"));
    assert!(
        r_text.starts_with("#NEXUS\n[PREMISE report v4]\n"),
        "{r_text}"
    );
    assert_eq!(premise_line(&r_text, "VERSION"), "4");
    check_blocks(&r_text);
    assert_eq!(
        block(&r_text, "ALIGNMENT_LIKELIHOODS"),
        block(&q_text, "ALIGNMENT_LIKELIHOODS")
    );
    let r_matches = read(&with_ext(&r, "matches"));
    check_posteriors(&r_text, &r_matches, 0.6);
    let summary = matches_summary(&r_matches, "0.6");
    assert_summary_params(&summary, &r_text);
    assert!(summary.contains("rho=3.5e0"), "{summary}");
    assert!(summary.contains(", refit of q.report; "), "{summary}");
    assert!(
        premise_line(&r_text, "EM").contains("collapsed=1"),
        "{r_text}"
    );
    // The refit pooled both references under species 11: the reads keep their posterior
    // on their own sequence, all of it now under 11.
    let tree = Tree::of(&r_text);
    let eleven = &tree.id["11"];
    for (id, rows) in posteriors(&r_text) {
        assert!(
            (tree.mass(eleven, &rows) - 1.0).abs() < 1e-9,
            "{id}: {rows:?}"
        );
    }
    assert_ne!(
        block(&r_text, "NODES"),
        block(&q_text, "NODES"),
        "the refit's own tree"
    );
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
        "--write-report",
    ]);
    ok(&out);
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(stdout.contains("Refit of q.report"), "{stdout}");
    assert!(stdout.contains("same.report"), "{stdout}");
    // The stored fit's .props is the query's own; the refit's agrees with it up to float
    // summation order.
    assert_eq!(
        read(&with_ext(&base, "props")),
        read(&with_ext(&q, "props"))
    );
    assert_same_props(
        &read(&with_ext(&base, "props")),
        &read(&with_ext(&same, "props")),
    );
    // Same final vertices, parents, ranks and contig counts in NODES.
    let q_text = read(&q_report);
    let r_text = read(&with_ext(&same, "report"));
    assert_eq!(final_nodes(&r_text), final_nodes(&q_text));
    // Same assignments as the stored fit (below the summaries, which name the fits).
    let two_cols = |t: &str| -> Vec<String> {
        table_lines(t)
            .map(|l| l.split('\t').take(2).collect::<Vec<_>>().join("\t"))
            .collect()
    };
    let (same_m, base_m) = (
        read(&with_ext(&same, "matches")),
        read(&with_ext(&base, "matches")),
    );
    assert_eq!(two_cols(&same_m), two_cols(&base_m));
    // The stored fit's summary is the query's own; the refit's has the same settings and
    // names its source, in .matches and in .props alike.
    assert_eq!(
        matches_summary(&base_m, "0.6"),
        matches_summary(&read(&with_ext(&q, "matches")), "0.6")
    );
    let (sb, ss) = (
        matches_summary(&base_m, "0.6"),
        matches_summary(&same_m, "0.6"),
    );
    assert!(!sb.contains("refit of"), "{sb}");
    assert!(ss.contains(", refit of q.report; gamma=0.6"), "{ss}");
    let settings = |s: &str| s.split(';').next().unwrap().to_string();
    assert_eq!(settings(&sb), settings(&ss));
    assert_same_summary(&base_m, &read(&with_ext(&base, "props")), "0.6");
    assert_same_summary(&same_m, &read(&with_ext(&same, "props")), "0.6");
    // Same read posteriors: the same sequences per read, the same values up to float
    // summation order; and each report's posteriors explain its own .matches.
    let (pq, pr) = (posteriors(&q_text), posteriors(&r_text));
    assert_eq!(pq.len(), pr.len());
    for ((ia, ra), (ib, rb)) in pq.iter().zip(&pr) {
        assert_eq!(ia, ib);
        assert_eq!(ra.len(), rb.len(), "{ia}: {ra:?} vs {rb:?}");
        for ((sa, a), (sb, b)) in ra.iter().zip(rb) {
            assert_eq!(sa, sb, "{ia}");
            assert!((a - b).abs() <= 1e-9, "{ia}: {ra:?} vs {rb:?}");
        }
    }
    check_posteriors(&q_text, &read(&with_ext(&base, "matches")), 0.6);
    check_posteriors(&r_text, &read(&with_ext(&same, "matches")), 0.6);
    // The refit report: lineage, the same EM, the alignment rows copied through.
    assert_eq!(premise_line(&r_text, "SOURCE"), "q.report");
    assert_eq!(premise_line(&r_text, "VERSION"), "4");
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
    ]));
    for ext in ["matches", "props"] {
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
    let leaves: Vec<String> = props_map(&props).into_iter().map(|l| l.0).collect();
    assert_eq!(leaves, ["ref_A", "ref_B"], "{props}");
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
    fs::write(with_ext(&b, "report"), "old\n").unwrap();
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
    assert_eq!(read(&with_ext(&b, "report")), "old\n");
    assert!(!with_ext(&b, "props").exists() && !with_ext(&b, "matches").exists());
}
