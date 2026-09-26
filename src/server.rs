use crate::io::{build_index_from_bytes, QueryWriters};
use crate::utils::{EMProb, QueryProgress};
use crate::{run_alignment, run_query};
use anyhow::Result;
use std::collections::HashMap;
use std::fs::File;
use std::sync::atomic::Ordering;
use std::sync::{Arc, Mutex};
use std::thread;
use tiny_http::{Header, Response, Server, StatusCode};

/// Serve the web interface on `addr` (`ip:port`) until the process is stopped.
pub fn serve(addr: &str) -> Result<()> {
    let html = include_str!("templates/index.html")
        .replace("{{CSS}}", include_str!("templates/styles.css"))
        .replace("{{JS}}", include_str!("templates/app.js"));

    let server =
        Server::http(addr).map_err(|e| anyhow::anyhow!("Failed to bind to {}: {}", addr, e))?;

    println!("Premise web interface running at http://{}", addr);
    println!("Press Ctrl+C to stop.");

    let mut sessions: HashMap<String, AlignSession> = HashMap::new();
    let mut query_sessions: HashMap<String, AlignSession> = HashMap::new();
    let progress_map: Arc<Mutex<HashMap<String, Arc<QueryProgress>>>> =
        Arc::new(Mutex::new(HashMap::new()));

    for request in server.incoming_requests() {
        let url = request.url().to_string();
        let path = url.split('?').next().unwrap_or("/").to_string();
        let method = request.method().clone();

        match (method, path.as_str()) {
            (tiny_http::Method::Get, "/") => {
                let ct = Header::from_bytes(b"Content-Type", b"text/html; charset=utf-8").unwrap();
                let _ = request.respond(Response::from_string(&html).with_header(ct));
            }
            (tiny_http::Method::Post, "/api/build") => {
                if std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    handle_api_build(request)
                }))
                .is_err()
                {
                    eprintln!("Error: /api/build handler panicked; request dropped");
                }
            }
            (tiny_http::Method::Post, "/api/align/upload") => {
                if std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    handle_align_upload(request, &mut sessions)
                }))
                .is_err()
                {
                    eprintln!("Error: /api/align/upload handler panicked; request dropped");
                }
            }
            (tiny_http::Method::Post, "/api/align/run") => {
                if std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    handle_align_run(request, &sessions)
                }))
                .is_err()
                {
                    eprintln!("Error: /api/align/run handler panicked; request dropped");
                }
            }
            (tiny_http::Method::Post, "/api/query/upload") => {
                if std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    handle_align_upload(request, &mut query_sessions)
                }))
                .is_err()
                {
                    eprintln!("Error: /api/query/upload handler panicked; request dropped");
                }
            }
            (tiny_http::Method::Post, "/api/query/run") => {
                if std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    handle_query_run(request, &query_sessions, progress_map.clone())
                }))
                .is_err()
                {
                    eprintln!("Error: /api/query/run handler panicked; request dropped");
                }
            }
            (tiny_http::Method::Get, "/api/query/progress") => {
                if std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    handle_query_progress(request, &progress_map)
                }))
                .is_err()
                {
                    eprintln!("Error: /api/query/progress handler panicked; request dropped");
                }
            }
            _ => {
                let _ = request
                    .respond(Response::from_string("Not Found").with_status_code(StatusCode(404)));
            }
        }
    }
    Ok(())
}

/// HTTP handler for `POST /api/build`.
fn handle_api_build(mut request: tiny_http::Request) {
    let result: Result<(Vec<u8>, String)> = (|| {
        let mut body = Vec::new();
        request.as_reader().read_to_end(&mut body)?;
        if body.is_empty() {
            return Err(anyhow::anyhow!("Empty request body"));
        }
        build_index_from_bytes(&body)
    })();

    match result {
        Ok((data, log_str)) => {
            println!("{}", log_str);
            let log_header_val = log_str.replace('\n', " | ");
            let ct = Header::from_bytes(b"Content-Type", b"application/octet-stream").unwrap();
            let cd = Header::from_bytes(
                b"Content-Disposition",
                b"attachment; filename=\"output.fmidx\"",
            )
            .unwrap();
            let lg = Header::from_bytes(b"X-Premise-Log", log_header_val.as_bytes())
                .unwrap_or_else(|_| Header::from_bytes(b"X-Premise-Log", b"").unwrap());
            let _ = request.respond(
                Response::from_data(data)
                    .with_header(ct)
                    .with_header(cd)
                    .with_header(lg),
            );
        }
        Err(e) => {
            let _ = request
                .respond(Response::from_string(e.to_string()).with_status_code(StatusCode(500)));
        }
    }
}

/// Server-side state for a single upload/run session.
struct AlignSession {
    dir: std::path::PathBuf,
    r1_ext: Option<String>,
    r2_ext: Option<String>,
}

/// Generate a unique session ID based on the current system time in nanoseconds.
fn new_session_id() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    format!(
        "{:x}",
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    )
}

/// Parse the query-string portion of a URL into a key→value map.
fn parse_qs(url: &str) -> HashMap<String, String> {
    url.split('?')
        .nth(1)
        .unwrap_or("")
        .split('&')
        .filter_map(|pair| {
            let mut it = pair.splitn(2, '=');
            let k = it.next()?.to_string();
            let v = it.next().unwrap_or("").to_string();
            Some((k, v))
        })
        .collect()
}

/// Handle a file-upload request for the alignment workflow.
fn do_align_upload(
    request: &mut tiny_http::Request,
    sessions: &mut HashMap<String, AlignSession>,
) -> Result<String> {
    let url = request.url().to_string();
    let qs = parse_qs(&url);
    let part = qs.get("part").map(|s| s.as_str()).unwrap_or("").to_string();
    let ext = qs
        .get("ext")
        .cloned()
        .unwrap_or_else(|| "fastq".to_string());

    let session_id =
        if qs.get("session").map(|s| s.as_str()) == Some("new") || !qs.contains_key("session") {
            let id = new_session_id();
            let dir = std::env::temp_dir().join(format!("premise_{}", id));
            std::fs::create_dir_all(&dir)?;
            sessions.insert(
                id.clone(),
                AlignSession {
                    dir,
                    r1_ext: None,
                    r2_ext: None,
                },
            );
            id
        } else {
            qs.get("session").unwrap().clone()
        };

    let session = sessions
        .get_mut(&session_id)
        .ok_or_else(|| anyhow::anyhow!("Session not found: {}", session_id))?;

    let file_name = match part.as_str() {
        "index" => "index.fmidx".to_string(),
        "r1" => {
            session.r1_ext = Some(ext.clone());
            format!("r1.{}", ext)
        }
        "r2" => {
            session.r2_ext = Some(ext.clone());
            format!("r2.{}", ext)
        }
        _ => return Err(anyhow::anyhow!("Unknown upload part: {}", part)),
    };

    let file_path = session.dir.join(&file_name);
    let mut out_file = File::create(&file_path)?;
    std::io::copy(request.as_reader(), &mut out_file)?;

    Ok(format!(r#"{{"session":"{}","ok":true}}"#, session_id))
}

/// HTTP handler for `POST /api/align/upload`
fn handle_align_upload(
    mut request: tiny_http::Request,
    sessions: &mut HashMap<String, AlignSession>,
) {
    match do_align_upload(&mut request, sessions) {
        Ok(json) => {
            let ct = Header::from_bytes(b"Content-Type", b"application/json").unwrap();
            let _ = request.respond(Response::from_string(json).with_header(ct));
        }
        Err(e) => {
            let _ = request
                .respond(Response::from_string(e.to_string()).with_status_code(StatusCode(500)));
        }
    }
}

/// Parse query parameters and run pairwise alignment for an uploaded session.
fn do_align_run(
    request: &tiny_http::Request,
    sessions: &HashMap<String, AlignSession>,
) -> Result<String> {
    let url = request.url().to_string();
    let qs = parse_qs(&url);
    let session_id = qs
        .get("session")
        .ok_or_else(|| anyhow::anyhow!("Missing session"))?;
    let mem_seed_length: usize = qs
        .get("mem_seed_length")
        .and_then(|s| s.parse().ok())
        .unwrap_or(22);
    let eps_2: EMProb = qs
        .get("eps_2")
        .and_then(|s| s.parse().ok())
        .unwrap_or(1e-18);
    let threads: usize = qs.get("threads").and_then(|s| s.parse().ok()).unwrap_or(0);

    let session = sessions
        .get(session_id)
        .ok_or_else(|| anyhow::anyhow!("Session not found: {}", session_id))?;

    let r1_ext = session.r1_ext.as_deref().unwrap_or("fastq");
    let r2_ext = session.r2_ext.as_deref().unwrap_or("fastq");

    let ref_path = session.dir.join("index.fmidx");
    let r1_path = session.dir.join(format!("r1.{}", r1_ext));
    let r2_path = session.dir.join(format!("r2.{}", r2_ext));

    let ref_path_str = ref_path
        .to_str()
        .ok_or_else(|| anyhow::anyhow!("session path is not valid UTF-8"))?;
    let r1_path_str = r1_path
        .to_str()
        .ok_or_else(|| anyhow::anyhow!("session path is not valid UTF-8"))?;
    let r2_path_str = r2_path
        .to_str()
        .ok_or_else(|| anyhow::anyhow!("session path is not valid UTF-8"))?;

    let mut tsv: Vec<u8> = Vec::new();
    let log_str = run_alignment(
        ref_path_str,
        r1_path_str,
        r2_path_str,
        mem_seed_length,
        eps_2,
        threads,
        &mut tsv,
    )?;
    Ok(format!(
        r#"{{"ok":true,"tsv":"{}","log":"{}"}}"#,
        json_escape(&String::from_utf8_lossy(&tsv)),
        json_escape(&log_str),
    ))
}

/// HTTP handler for `POST /api/align/run`
fn handle_align_run(request: tiny_http::Request, sessions: &HashMap<String, AlignSession>) {
    match do_align_run(&request, sessions) {
        Ok(json) => {
            let ct = Header::from_bytes(b"Content-Type", b"application/json").unwrap();
            let _ = request.respond(Response::from_string(json).with_header(ct));
        }
        Err(e) => {
            let _ = request
                .respond(Response::from_string(e.to_string()).with_status_code(StatusCode(500)));
        }
    }
}

/// Escape a string for embedding inside a JSON double-quoted value.
fn json_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c => out.push(c),
        }
    }
    out
}

/// HTTP handler for `GET /api/query/progress`
fn handle_query_progress(
    request: tiny_http::Request,
    progress_map: &Arc<Mutex<HashMap<String, Arc<QueryProgress>>>>,
) {
    let url = request.url().to_string();
    let qs = parse_qs(&url);
    let json = match qs.get("session") {
        Some(sid) => {
            let map = progress_map.lock().unwrap();
            match map.get(sid) {
                Some(p) => p.to_json(),
                None => r#"{"phase":0,"reads_done":0,"reads_total":0,"em_iter_done":0,"em_iter_total":0,"elapsed_ms":0}"#.to_string(),
            }
        }
        None => r#"{"phase":0,"reads_done":0,"reads_total":0,"em_iter_done":0,"em_iter_total":0,"elapsed_ms":0}"#.to_string(),
    };
    let ct = Header::from_bytes(b"Content-Type", b"application/json").unwrap();
    let cors = Header::from_bytes(b"Access-Control-Allow-Origin", b"*").unwrap();
    let _ = request.respond(
        Response::from_string(json)
            .with_header(ct)
            .with_header(cors),
    );
}

/// HTTP handler for `POST /api/query/run`
fn handle_query_run(
    request: tiny_http::Request,
    sessions: &HashMap<String, AlignSession>,
    progress_map: Arc<Mutex<HashMap<String, Arc<QueryProgress>>>>,
) {
    let url = request.url().to_string();
    let qs = parse_qs(&url);

    let session_id = match qs.get("session") {
        Some(s) => s.clone(),
        None => {
            let _ = request.respond(
                Response::from_string("Missing session").with_status_code(StatusCode(400)),
            );
            return;
        }
    };

    let mem_seed_length: usize = qs
        .get("mem_seed_length")
        .and_then(|s| s.parse().ok())
        .unwrap_or(22);
    let eps_1: EMProb = qs.get("eps_1").and_then(|s| s.parse().ok()).unwrap_or(0.0);
    let eps_2: EMProb = qs
        .get("eps_2")
        .and_then(|s| s.parse().ok())
        .unwrap_or(1e-18);
    let num_iter: usize = qs.get("iter").and_then(|s| s.parse().ok()).unwrap_or(100);
    let rho: EMProb = qs.get("rho").and_then(|s| s.parse().ok()).unwrap_or(150.0);
    let omega: EMProb = qs
        .get("omega")
        .and_then(|s| s.parse().ok())
        .unwrap_or(1e-10);
    let em_threshold: EMProb = qs
        .get("em_threshold")
        .and_then(|s| s.parse().ok())
        .unwrap_or(1e-6);
    let no_penalty: bool = qs.get("no_penalty").map(|s| s == "true").unwrap_or(false);
    let use_penalty = !no_penalty;
    let threads: usize = qs.get("threads").and_then(|s| s.parse().ok()).unwrap_or(0);

    let session = match sessions.get(&session_id) {
        Some(s) => s,
        None => {
            let _ = request.respond(
                Response::from_string(format!("Session not found: {}", session_id))
                    .with_status_code(StatusCode(404)),
            );
            return;
        }
    };

    let r1_ext = session.r1_ext.as_deref().unwrap_or("fastq").to_string();
    let r2_ext = session.r2_ext.as_deref().unwrap_or("fastq").to_string();

    let paths = [
        session.dir.join("index.fmidx"),
        session.dir.join(format!("r1.{}", r1_ext)),
        session.dir.join(format!("r2.{}", r2_ext)),
    ];
    let path_strs: Option<Vec<String>> = paths
        .iter()
        .map(|p| p.to_str().map(|s| s.to_string()))
        .collect();
    let path_strs = match path_strs {
        Some(v) => v,
        None => {
            eprintln!("Error: session path is not valid UTF-8");
            let _ = request.respond(
                Response::from_string("session path is not valid UTF-8")
                    .with_status_code(StatusCode(500)),
            );
            return;
        }
    };
    let ref_path = path_strs[0].clone();
    let r1_path = path_strs[1].clone();
    let r2_path = path_strs[2].clone();

    // Register fresh progress entry for this session
    let progress = Arc::new(QueryProgress::new());
    progress_map
        .lock()
        .unwrap()
        .insert(session_id.clone(), progress.clone());

    // Run query in background thread; main server loop can serve progress polls
    thread::spawn(move || {
        let mut tables = QueryWriters {
            matches: Vec::new(),
            posteriors: Vec::new(),
            props: Vec::new(),
            aligns: Vec::new(),
        };
        let result = run_query(
            &ref_path,
            &r1_path,
            &r2_path,
            mem_seed_length,
            eps_1,
            eps_2,
            num_iter,
            rho,
            omega,
            em_threshold,
            use_penalty,
            threads,
            Some(progress.clone()),
            &mut tables,
        );

        let json = match result {
            Ok((em_likelihoods, log_str)) => {
                progress.phase.store(3, Ordering::Relaxed);
                let matches_tsv = String::from_utf8_lossy(&tables.matches);
                let posteriors_tsv = String::from_utf8_lossy(&tables.posteriors);
                let props_tsv = String::from_utf8_lossy(&tables.props);
                let aligns_tsv = String::from_utf8_lossy(&tables.aligns);
                let convergence_json = format!(
                    "[{}]",
                    em_likelihoods
                        .iter()
                        .map(|v| if v.is_finite() {
                            format!("{:.10e}", v)
                        } else {
                            "null".to_string()
                        })
                        .collect::<Vec<_>>()
                        .join(",")
                );
                format!(
                    r#"{{"ok":true,"matches":"{}","posteriors":"{}","props":"{}","aligns":"{}","convergence":{},"log":"{}"}}"#,
                    json_escape(&matches_tsv),
                    json_escape(&posteriors_tsv),
                    json_escape(&props_tsv),
                    json_escape(&aligns_tsv),
                    convergence_json,
                    json_escape(&log_str),
                )
            }
            Err(e) => {
                format!(
                    r#"{{"ok":false,"error":"{}"}}"#,
                    json_escape(&e.to_string())
                )
            }
        };

        let ct = Header::from_bytes(b"Content-Type", b"application/json").unwrap();
        let _ = request.respond(Response::from_string(json).with_header(ct));
    });
}
