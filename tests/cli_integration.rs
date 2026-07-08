/// Integration / E2E tests for the `premise` CLI.
///
/// Each test spawns the compiled binary as a child process, drives it with
/// synthetic fixture files, and asserts on exit codes, stdout/stderr content,
/// and the shape of output files.
///
/// Run with:  cargo test --test cli_integration
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

// ─────────────────────────── helpers ────────────────────────────────────────

/// Absolute path to the compiled `premise` binary (debug build).
fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("target")
        .join("debug")
        .join("premise")
}

/// Absolute path to the fixtures directory.
fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
}

/// Run `premise <args>` and return the `std::process::Output`.
fn run(args: &[&str]) -> std::process::Output {
    Command::new(bin())
        .args(args)
        .output()
        .expect("failed to spawn premise binary")
}

/// Assert that the binary exited with code 0.
fn assert_success(output: &std::process::Output) {
    assert!(
        output.status.success(),
        "expected exit 0, got {:?}\nstdout: {}\nstderr: {}",
        output.status.code(),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
}

/// Assert that the binary exited with a non-zero code.
fn assert_failure(output: &std::process::Output) {
    assert!(
        !output.status.success(),
        "expected non-zero exit, but process succeeded\nstdout: {}\nstderr: {}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
}

/// Return stdout as a UTF-8 string.
fn stdout(output: &std::process::Output) -> String {
    String::from_utf8_lossy(&output.stdout).into_owned()
}

/// Build an FM-index from the fixture FASTA into a temp directory.
///
/// The binary writes the index to exactly the path given via `-o`, so we pass
/// `<tmpdir>/ref` and expect the file `<tmpdir>/ref` to exist afterwards.
/// Returns the path to the generated index file.
fn build_index(tmpdir: &Path) -> PathBuf {
    let fasta = fixtures().join("ref.fasta");
    let outpath = tmpdir.join("ref");
    let output = run(&[
        "build",
        "-s",
        fasta.to_str().unwrap(),
        "-o",
        outpath.to_str().unwrap(),
    ]);
    assert_success(&output);
    assert!(
        outpath.exists(),
        "expected FM-index at {} after build\nstdout: {}\nstderr: {}",
        outpath.display(),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
    outpath
}

// ─────────────────────────── global ─────────────────────────────────────────

/// `premise --help` should exit 0 and mention the available subcommands.
#[test]
fn help_exits_zero_and_lists_subcommands() {
    let output = run(&["--help"]);
    assert_success(&output);
    let out = stdout(&output);
    for cmd in &["build", "query", "server"] {
        assert!(
            out.contains(cmd),
            "help output missing expected subcommand '{}'\nfull output:\n{}",
            cmd,
            out
        );
    }
}

/// `premise --version` should exit 0 and print a version string.
#[test]
fn version_exits_zero() {
    let output = run(&["--version"]);
    assert_success(&output);
    let out = stdout(&output);
    assert!(!out.trim().is_empty(), "expected non-empty version string");
}

/// Invoking `premise` with no subcommand prints help and exits 0.
#[test]
fn no_args_exits_zero() {
    let output = run(&[]);
    assert_success(&output);
}

// ─────────────────────────── build subcommand ────────────────────────────────

/// `premise build --help` exits 0 and describes the subcommand.
#[test]
fn build_help_exits_zero() {
    let output = run(&["build", "--help"]);
    assert_success(&output);
    assert!(stdout(&output).contains("FM-index"));
}

/// `premise build` without the required `--source` flag must fail.
#[test]
fn build_missing_source_fails() {
    let output = run(&["build"]);
    assert_failure(&output);
}

/// A complete `build` run against the synthetic fixture FASTA should produce
/// an index file at the specified output path.
#[test]
fn build_produces_index_file() {
    let tmpdir = tempfile::tempdir().expect("could not create temp dir");
    build_index(tmpdir.path()); // asserts internally
}

/// Running `build` twice with the same output path should succeed (idempotent).
#[test]
fn build_is_idempotent() {
    let tmpdir = tempfile::tempdir().expect("could not create temp dir");
    build_index(tmpdir.path());
    build_index(tmpdir.path());
}

/// The `build` stdout should report where the index was written.
#[test]
fn build_stdout_reports_output_path() {
    let tmpdir = tempfile::tempdir().expect("could not create temp dir");
    let fasta = fixtures().join("ref.fasta");
    let outpath = tmpdir.path().join("myindex");
    let output = run(&[
        "build",
        "-s",
        fasta.to_str().unwrap(),
        "-o",
        outpath.to_str().unwrap(),
    ]);
    assert_success(&output);
    let out = stdout(&output);
    assert!(
        out.contains("Index written"),
        "expected 'Index written' in stdout\nactual: {}",
        out
    );
}

/// Passing a non-existent FASTA file to `build` should exit non-zero.
#[test]
fn build_nonexistent_fasta_fails() {
    let tmpdir = tempfile::tempdir().expect("could not create temp dir");
    let outpath = tmpdir.path().join("out");
    let output = run(&[
        "build",
        "-s",
        "/nonexistent/path/ref.fasta",
        "-o",
        outpath.to_str().unwrap(),
    ]);
    assert_failure(&output);
}

/// When no `-o` is given, `build` derives the output name from the source FASTA
/// (appends `.fmidx`). Write source to a temp directory so the auto-named output
/// lands there.
#[test]
fn build_default_output_name_appends_fmidx() {
    let tmpdir = tempfile::tempdir().expect("could not create temp dir");
    // Copy fixture FASTA into tmpdir so the auto-named .fmidx is created there.
    let src = fixtures().join("ref.fasta");
    let dst = tmpdir.path().join("ref.fasta");
    fs::copy(&src, &dst).expect("could not copy fixture fasta");

    let output = run(&["build", "-s", dst.to_str().unwrap()]);
    assert_success(&output);

    let expected = tmpdir.path().join("ref.fasta.fmidx");
    assert!(
        expected.exists(),
        "expected auto-named index at {}",
        expected.display()
    );
}

// ─────────────────────────── query subcommand ────────────────────────────────

/// `premise query --help` exits 0 and lists required flags.
#[test]
fn query_help_exits_zero_and_lists_flags() {
    let output = run(&["query", "--help"]);
    assert_success(&output);
    let out = stdout(&output);
    for flag in &["--source", "--r1", "--r2", "--mem"] {
        assert!(
            out.contains(flag),
            "query help missing expected flag '{}'\nfull output:\n{}",
            flag,
            out
        );
    }
}

/// `premise query` without required arguments must fail.
#[test]
fn query_missing_required_args_fails() {
    let output = run(&["query"]);
    assert_failure(&output);
}

/// Full end-to-end query: build index then run query with synthetic reads.
/// Verifies that exit code is 0 and that all four expected output files are
/// created and non-empty.
#[test]
fn query_end_to_end_produces_output_files() {
    let tmpdir = tempfile::tempdir().expect("could not create temp dir");
    let index = build_index(tmpdir.path());
    let r1 = fixtures().join("r1.fastq");
    let r2 = fixtures().join("r2.fastq");
    let outprefix = tmpdir.path().join("results");

    let output = run(&[
        "query",
        "-s",
        index.to_str().unwrap(),
        "-1",
        r1.to_str().unwrap(),
        "-2",
        r2.to_str().unwrap(),
        "-p",
        "5",
        "-o",
        outprefix.to_str().unwrap(),
        "-t",
        "1",
        "-i",
        "5", // fewer EM iterations for test speed
    ]);

    assert_success(&output);

    // All four output files must exist. The `.props` file can be empty when EM
    // collapses all proportions to zero, so we only check existence for it.
    for ext in &["matches", "posteriors", "props", "aligns"] {
        let path = tmpdir.path().join(format!("results.{}", ext));
        assert!(
            path.exists(),
            "expected output file results.{} to exist",
            ext
        );
    }

    // matches, posteriors and aligns always have at least a header row
    for ext in &["matches", "posteriors", "aligns"] {
        let path = tmpdir.path().join(format!("results.{}", ext));
        let content = fs::read_to_string(&path)
            .unwrap_or_else(|_| panic!("could not read {}", path.display()));
        assert!(
            !content.trim().is_empty(),
            "output file results.{} should not be empty",
            ext
        );
    }
}

/// The `.props` output file must exist after a successful query.
///
/// Note: the props file has no header row. Its content may be empty when the
/// EM assigns zero proportion to all references (a valid outcome for small
/// datasets under the default penalized EM). We only assert the file exists.
#[test]
fn query_props_file_exists_after_query() {
    let tmpdir = tempfile::tempdir().expect("could not create temp dir");
    let index = build_index(tmpdir.path());
    let r1 = fixtures().join("r1.fastq");
    let r2 = fixtures().join("r2.fastq");
    let outprefix = tmpdir.path().join("results");

    let output = run(&[
        "query",
        "-s",
        index.to_str().unwrap(),
        "-1",
        r1.to_str().unwrap(),
        "-2",
        r2.to_str().unwrap(),
        "-p",
        "5",
        "-o",
        outprefix.to_str().unwrap(),
        "-t",
        "1",
        "-i",
        "5",
    ]);
    assert_success(&output);

    // File must exist (may be empty — EM can collapse all proportions to 0).
    let props_path = tmpdir.path().join("results.props");
    assert!(
        props_path.exists(),
        "expected results.props to be created by query"
    );
    // Any non-empty lines that appear must contain a ref_id and a numeric proportion.
    let content = fs::read_to_string(&props_path).expect("could not read results.props");
    for line in content.lines() {
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(
            cols.len(),
            2,
            "expected 2 tab-separated columns in props row, got: {}",
            line
        );
        cols[1].parse::<f64>().unwrap_or_else(|_| {
            panic!("expected numeric proportion in props row, got: {}", cols[1])
        });
    }
}

/// The `.matches` output should contain a header row and at least one data row.
#[test]
fn query_matches_file_has_data_rows() {
    let tmpdir = tempfile::tempdir().expect("could not create temp dir");
    let index = build_index(tmpdir.path());
    let r1 = fixtures().join("r1.fastq");
    let r2 = fixtures().join("r2.fastq");
    let outprefix = tmpdir.path().join("results");

    let output = run(&[
        "query",
        "-s",
        index.to_str().unwrap(),
        "-1",
        r1.to_str().unwrap(),
        "-2",
        r2.to_str().unwrap(),
        "-p",
        "5",
        "-o",
        outprefix.to_str().unwrap(),
        "-t",
        "1",
        "-i",
        "5",
    ]);
    assert_success(&output);

    let matches_path = tmpdir.path().join("results.matches");
    let content = fs::read_to_string(&matches_path).expect("could not read results.matches");
    let lines: Vec<&str> = content.lines().collect();
    assert!(
        !lines.is_empty(),
        "expected at least a header row in results.matches"
    );
}

/// The `.posteriors` output should have one row per read plus a header.
#[test]
fn query_posteriors_file_has_data_rows() {
    let tmpdir = tempfile::tempdir().expect("could not create temp dir");
    let index = build_index(tmpdir.path());
    let r1 = fixtures().join("r1.fastq");
    let r2 = fixtures().join("r2.fastq");
    let outprefix = tmpdir.path().join("results");

    let output = run(&[
        "query",
        "-s",
        index.to_str().unwrap(),
        "-1",
        r1.to_str().unwrap(),
        "-2",
        r2.to_str().unwrap(),
        "-p",
        "5",
        "-o",
        outprefix.to_str().unwrap(),
        "-t",
        "1",
        "-i",
        "5",
    ]);
    assert_success(&output);

    let path = tmpdir.path().join("results.posteriors");
    let content = fs::read_to_string(&path).expect("could not read results.posteriors");
    assert!(
        !content.trim().is_empty(),
        "posteriors file should not be empty"
    );
}

/// The stdout of `query` should mention the output prefix.
#[test]
fn query_stdout_reports_output_path() {
    let tmpdir = tempfile::tempdir().expect("could not create temp dir");
    let index = build_index(tmpdir.path());
    let r1 = fixtures().join("r1.fastq");
    let r2 = fixtures().join("r2.fastq");
    let outprefix = tmpdir.path().join("results");

    let output = run(&[
        "query",
        "-s",
        index.to_str().unwrap(),
        "-1",
        r1.to_str().unwrap(),
        "-2",
        r2.to_str().unwrap(),
        "-p",
        "5",
        "-o",
        outprefix.to_str().unwrap(),
        "-t",
        "1",
        "-i",
        "5",
    ]);
    assert_success(&output);
    let out = stdout(&output);
    assert!(
        out.contains("Query written"),
        "expected 'Query written' in stdout\nactual: {}",
        out
    );
}

/// `query` with a non-existent index file should fail gracefully.
#[test]
fn query_nonexistent_index_fails() {
    let tmpdir = tempfile::tempdir().expect("could not create temp dir");
    let r1 = fixtures().join("r1.fastq");
    let r2 = fixtures().join("r2.fastq");
    let output = run(&[
        "query",
        "-s",
        "/nonexistent/ref.fmidx",
        "-1",
        r1.to_str().unwrap(),
        "-2",
        r2.to_str().unwrap(),
        "-p",
        "5",
        "-o",
        tmpdir.path().join("out").to_str().unwrap(),
        "-t",
        "1",
        "-i",
        "5",
    ]);
    assert_failure(&output);
}

/// `query` with a non-existent R1 file should fail gracefully.
#[test]
fn query_nonexistent_r1_fails() {
    let tmpdir = tempfile::tempdir().expect("could not create temp dir");
    let index = build_index(tmpdir.path());
    let r2 = fixtures().join("r2.fastq");
    let output = run(&[
        "query",
        "-s",
        index.to_str().unwrap(),
        "-1",
        "/nonexistent/r1.fastq",
        "-2",
        r2.to_str().unwrap(),
        "-p",
        "5",
        "-o",
        tmpdir.path().join("out").to_str().unwrap(),
        "-t",
        "1",
        "-i",
        "5",
    ]);
    assert_failure(&output);
}

/// `query` with `--no-penalty` flag should still succeed and produce output.
#[test]
fn query_no_penalty_flag_succeeds() {
    let tmpdir = tempfile::tempdir().expect("could not create temp dir");
    let index = build_index(tmpdir.path());
    let r1 = fixtures().join("r1.fastq");
    let r2 = fixtures().join("r2.fastq");
    let outprefix = tmpdir.path().join("results_nopenalty");

    let output = run(&[
        "query",
        "-s",
        index.to_str().unwrap(),
        "-1",
        r1.to_str().unwrap(),
        "-2",
        r2.to_str().unwrap(),
        "-p",
        "5",
        "-o",
        outprefix.to_str().unwrap(),
        "-t",
        "1",
        "-i",
        "5",
        "--no-penalty",
    ]);
    assert_success(&output);
    assert!(tmpdir.path().join("results_nopenalty.props").exists());
}

// ─────────────────────────── server subcommand ───────────────────────────────

/// `premise server --help` exits 0 and mentions port and ip options.
#[test]
fn server_help_exits_zero() {
    let output = run(&["server", "--help"]);
    assert_success(&output);
    let out = stdout(&output);
    assert!(
        out.contains("port") || out.contains("ip"),
        "server help missing expected options\nfull output:\n{}",
        out
    );
}

/// Starting the server on a free port should respond with HTTP 200 and an HTML
/// document on `GET /`.
#[test]
fn server_serves_html_on_root() {
    use std::net::TcpListener;
    use std::thread;
    use std::time::Duration;

    // Pick a free port.
    let listener = TcpListener::bind("127.0.0.1:0").expect("could not bind");
    let port = listener.local_addr().unwrap().port();
    drop(listener);

    let port_str = port.to_string();
    let mut child = Command::new(bin())
        .args(&["server", "--port", &port_str, "--ip", "127.0.0.1"])
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .expect("failed to spawn server");

    // Wait for the server to be ready: retry up to 10 s with 200ms intervals.
    // Under parallel test load the server may take longer to start.
    let url = format!("http://127.0.0.1:{}/", port);
    let mut response = Err(ureq::Error::from(std::io::Error::new(
        std::io::ErrorKind::ConnectionRefused,
        "not started yet",
    )));
    for _ in 0..50 {
        thread::sleep(Duration::from_millis(200));
        response = ureq::get(&url).call();
        if response.is_ok() {
            break;
        }
    }

    // Read the response body BEFORE killing the server — killing it first
    // closes the TCP connection and causes into_string() to return empty.
    let (status, body) = match response {
        Ok(resp) => {
            let status = resp.status();
            let body = resp.into_string().unwrap_or_default();
            (status, body)
        }
        Err(e) => {
            child.kill().ok();
            let _ = child.wait();
            panic!("HTTP request to server failed after retries: {}", e);
        }
    };

    child.kill().ok();
    let _ = child.wait();

    assert_eq!(status, 200, "expected HTTP 200 from GET /");
    assert!(
        body.contains("<!DOCTYPE") || body.contains("<html"),
        "expected HTML body from GET /\nactual body (first 300 chars):\n{}",
        &body[..body.len().min(300)]
    );
}

/// `GET /nonexistent` should return HTTP 404.
#[test]
fn server_returns_404_for_unknown_route() {
    use std::net::TcpListener;
    use std::thread;
    use std::time::Duration;

    let listener = TcpListener::bind("127.0.0.1:0").expect("could not bind");
    let port = listener.local_addr().unwrap().port();
    drop(listener);

    let port_str = port.to_string();
    let mut child = Command::new(bin())
        .args(&["server", "--port", &port_str, "--ip", "127.0.0.1"])
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .expect("failed to spawn server");

    // Wait for the server to accept connections (up to 10 s).
    let root_url = format!("http://127.0.0.1:{}/", port);
    for _ in 0..50 {
        thread::sleep(Duration::from_millis(200));
        if ureq::get(&root_url).call().is_ok() {
            break;
        }
    }

    let url = format!("http://127.0.0.1:{}/nonexistent-path", port);
    let response = ureq::get(&url).call();

    child.kill().ok();
    let _ = child.wait();

    match response {
        Ok(resp) => {
            assert_eq!(resp.status(), 404, "expected HTTP 404 for unknown route");
        }
        Err(ureq::Error::Status(404, _)) => {
            // ureq treats 4xx as errors — 404 is the expected outcome
        }
        Err(e) => {
            panic!("unexpected error from server: {}", e);
        }
    }
}

// ─────────────────── N-character alignment bug tests ────────────────────────
//
// These tests exercise the bug in `clean_kmer_matches` where the `enumerate()`
// call is placed AFTER the N-kmer filter, causing the read-position index to
// be wrong for any kmer that appears after a filtered (N-containing) kmer.
//
// For a 50 bp read with 5 % mismatch tolerance, kmer_size = 16.  A single N at
// position 15 filters all 16 kmers that span that position (starts 0-15), so
// the first valid kmer (actual position 16) is assigned enumerate index 0.
// This shifts `ref_pos = ref_start - read_start` from 0 (correct) to 16 (wrong).
//
// Consequence: the alignment probability at the wrong position is extremely low
// (~10^-161 with Phred-40 quality) and is dropped by the eps_1 filter, causing
// the read to appear "unclassified" even though it nearly perfectly matches ref_A.
//
// After the fix both tests below pass.

/// Full pipeline: a read with N at position 15 must still be classified as ref_A.
///
/// r1_with_n.fastq  = ref_A[0..50] with N substituted at position 15 (replacing T).
/// r2_with_n.fastq  = unchanged R2 mate from the original read1 pair.
///
/// With bug:    read1 → unclassified (alignment probability underflows to 0).
/// After fix:   read1 → ref_A, forward position 0.
#[test]
fn query_read_with_n_at_position_15_classified_as_ref_a() {
    let tmpdir = tempfile::tempdir().expect("could not create temp dir");
    let index = build_index(tmpdir.path());
    let r1 = fixtures().join("r1_with_n.fastq");
    let r2 = fixtures().join("r2_with_n.fastq");
    let outprefix = tmpdir.path().join("results_n");

    let output = run(&[
        "query",
        "-s",
        index.to_str().unwrap(),
        "-1",
        r1.to_str().unwrap(),
        "-2",
        r2.to_str().unwrap(),
        "-p",
        "5",
        "-o",
        outprefix.to_str().unwrap(),
        "-t",
        "1",
        "-i",
        "5",
    ]);
    assert_success(&output);

    let matches_path = tmpdir.path().join("results_n.matches");
    assert!(
        matches_path.exists(),
        "expected results_n.matches to be created"
    );

    let content = fs::read_to_string(&matches_path)
        .unwrap_or_else(|_| panic!("could not read results_n.matches"));

    let classified_as_ref_a = content.lines().any(|line| {
        let cols: Vec<&str> = line.split('\t').collect();
        cols.len() >= 2 && cols[0] == "read1" && cols[1] == "ref_A"
    });
    assert!(
        classified_as_ref_a,
        "BUG: read1 (ref_A[0..50] with N at position 15) was not classified as ref_A.\n\
         The N-kmer filtering bug likely caused the alignment at the wrong position \
         to have probability ≈ 0, leaving the read unclassified.\n\
         Full matches output:\n{}",
        content
    );
}

/// Full pipeline: the forward alignment position for a read with N must be 0.
///
/// Even if the read is correctly classified as ref_A, the bug causes the reported
/// forward position to be 16 instead of 0 — because ref_pos = ref_start - read_start
/// uses the wrong read_start from the shifted enumerate index.
#[test]
fn query_read_with_n_forward_position_is_zero() {
    let tmpdir = tempfile::tempdir().expect("could not create temp dir");
    let index = build_index(tmpdir.path());
    let r1 = fixtures().join("r1_with_n.fastq");
    let r2 = fixtures().join("r2_with_n.fastq");
    let outprefix = tmpdir.path().join("results_n_pos");

    let output = run(&[
        "query",
        "-s",
        index.to_str().unwrap(),
        "-1",
        r1.to_str().unwrap(),
        "-2",
        r2.to_str().unwrap(),
        "-p",
        "5",
        "-o",
        outprefix.to_str().unwrap(),
        "-t",
        "1",
        "-i",
        "5",
    ]);
    assert_success(&output);

    let matches_path = tmpdir.path().join("results_n_pos.matches");
    let content = fs::read_to_string(&matches_path)
        .unwrap_or_else(|_| panic!("could not read results_n_pos.matches"));

    let row = content.lines().find(|line| {
        let cols: Vec<&str> = line.split('\t').collect();
        cols.len() >= 2 && cols[0] == "read1" && cols[1] == "ref_A"
    });

    match row {
        None => panic!(
            "read1 was not classified as ref_A — cannot check position.\n\
             The classification test should also fail.\n\
             Full output:\n{}",
            content
        ),
        Some(line) => {
            let cols: Vec<&str> = line.split('\t').collect();
            assert!(
                cols.len() >= 4,
                "expected at least 4 columns in matches row, got: {}",
                line
            );
            assert_eq!(
                cols[3], "0",
                "BUG: expected forward position 0 (read maps to ref_A[0..50]), \
                 got position {}.\n\
                 The N-kmer enumerate() shift makes ref_pos = 16 instead of 0.",
                cols[3]
            );
        }
    }
}
