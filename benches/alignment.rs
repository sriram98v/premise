//! Criterion benchmarks for the PREMISE alignment stages (seed -> align).
//!
//! EM benches are intentionally out of scope (future work). All inputs are
//! synthetic and generated deterministically in-process (no disk, no RNG that
//! would vary run-to-run), so baselines are reproducible.
//!
//! Stages covered:
//!   scorer/compute_match_log_prob  pure per-base ungapped rescore (LUT win lands here)
//!   scorer/error_prob              isolates the Phred->error-prob transform (the powf)
//!   seed/clean_mem_matches         SMEM seeding + locate against the FM-index
//!   rescore/query_read             seed -> project -> ungapped rescore for one read
//!   merge/merge_read_pairs         FR position-map cross-product
//!   phase/process_read_pairs       full rayon-parallel align phase over many pairs
//!
//! Profiling: run with `--profile-time=<secs>` to emit an in-process pprof
//! profile (no `perf`/sudo needed). `PPROF_OUT=pb` writes a pprof protobuf
//! (`profile.pb`); otherwise a flamegraph SVG is written under
//! `target/criterion/<group>/<bench>/profile/`.

use bio::io::fastq;
use bio::stats::LogProb;
use criterion::{
    black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput,
};
use haystackfm::alphabet::encode_char;
use pprof::criterion::{Output, PProfProfiler};
use premise::{
    build_index_from_bytes, build_sparse_array, clean_mem_matches, merge_read_pairs,
    process_read_pairs, query_read, IOFMIndex, ReadIdx, ReadPair,
};
use std::collections::HashMap;

// ─── deterministic synthetic data ────────────────────────────────────────────

const READ_LEN: usize = 150;
const SEED_LEN: usize = 20; // representative --mem for 150 bp reads
const N_REFS: usize = 4;
const REF_LEN: usize = 2000;
const N_PAIRS: usize = 200;
const INSERT: usize = 200; // gap between R1 start and R2 start on the reference

/// Tiny LCG (Numerical Recipes constants) — deterministic, no external RNG.
fn lcg(state: &mut u64) -> u64 {
    *state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
    *state
}

fn random_dna(len: usize, state: &mut u64) -> Vec<u8> {
    const B: [u8; 4] = [b'A', b'C', b'G', b'T'];
    (0..len).map(|_| B[(lcg(state) >> 33) as usize & 3]).collect()
}

/// A `>refN`-headed FASTA of `N_REFS` pseudo-random references that all share a
/// single conserved block, so some seeds locate to multiple references (the
/// realistic metagenomic case that drives the locate cost).
fn synthetic_fasta() -> (Vec<u8>, Vec<Vec<u8>>) {
    let mut state = 0x9E3779B97F4A7C15u64; // fixed seed
    let conserved = random_dna(300, &mut state);
    let mut fasta = Vec::new();
    let mut refs = Vec::with_capacity(N_REFS);
    for i in 0..N_REFS {
        let mut seq = random_dna(REF_LEN, &mut state);
        // splice the shared conserved block into each reference at a fixed offset
        seq[500..800].copy_from_slice(&conserved);
        fasta.extend_from_slice(format!(">ref{i}\n").as_bytes());
        fasta.extend_from_slice(&seq);
        fasta.push(b'\n');
        refs.push(seq);
    }
    (fasta, refs)
}

fn rec(id: &str, seq: &[u8]) -> fastq::Record {
    fastq::Record::with_attrs(id, None, seq, &vec![b'I'; seq.len()])
}

fn main_index() -> IOFMIndex {
    let (fasta, _) = synthetic_fasta();
    let (bytes, _) = build_index_from_bytes(&fasta).expect("build index");
    IOFMIndex::from_bytes(&bytes).expect("load index")
}

/// `N_PAIRS` FR read pairs sampled from ref0: R1 forward at offset o, R2 the
/// reverse-complement of the window `INSERT` bp downstream.
fn synthetic_pairs(ref0: &[u8]) -> Vec<ReadPair> {
    let max_off = REF_LEN - INSERT - READ_LEN;
    (0..N_PAIRS)
        .map(|i| {
            let o = (i * 7) % max_off; // spread deterministically across the genome
            let r1 = ref0[o..o + READ_LEN].to_vec();
            let r2 = bio::alphabets::dna::revcomp(&ref0[o + INSERT..o + INSERT + READ_LEN]);
            ReadPair::new(&format!("r{i}"), rec("a", &r1), rec("b", &r2))
        })
        .collect()
}

// ─── benches ──────────────────────────────────────────────────────────────────

fn bench_alignment(c: &mut Criterion) {
    let idx = main_index();
    let refs = &idx.idx_to_seq;
    let h2r = &idx.header_to_idx;
    let ref0 = refs.values().next().unwrap().clone();

    // A representative 150 bp read taken verbatim from a reference (with 2 SNPs),
    // plus its raw + encoded forms and quality string.
    let mut read = ref0[100..100 + READ_LEN].to_vec();
    read[40] ^= 0b0000_0110; // flip a couple of bases so it isn't a perfect match
    read[110] ^= 0b0000_0110;
    let qual = vec![b'I'; READ_LEN];
    let ref_seg = ref0[100..100 + READ_LEN].to_vec();
    let q_enc: Vec<u8> = read.iter().filter_map(|&b| encode_char(b as char)).collect();
    let record = rec("read", &read);

    // scorer/compute_match_log_prob ───────────────────────────────────────────
    let mut g = c.benchmark_group("scorer");
    g.throughput(Throughput::Bytes(READ_LEN as u64));
    g.bench_function("compute_match_log_prob", |b| {
        b.iter(|| {
            premise::utils::compute_match_log_prob(
                black_box(&read),
                black_box(&qual),
                black_box(&ref_seg),
            )
        })
    });
    // scorer/error_prob — the powf, summed over a full read's quality string
    g.bench_function("error_prob", |b| {
        b.iter(|| {
            let mut acc = 0.0f64;
            for &q in black_box(&qual) {
                acc += *premise::utils::error_prob(q);
            }
            acc
        })
    });
    g.finish();

    // seed/clean_mem_matches ────────────────────────────────────────────────────
    let mut g = c.benchmark_group("seed");
    g.bench_function("clean_mem_matches", |b| {
        b.iter(|| {
            clean_mem_matches(
                black_box(&idx.fmidx),
                black_box(h2r),
                black_box(&q_enc),
                black_box(SEED_LEN),
            )
        })
    });
    g.finish();

    // rescore/query_read ────────────────────────────────────────────────────────
    let mut g = c.benchmark_group("rescore");
    g.bench_function("query_read", |b| {
        b.iter(|| {
            query_read(
                black_box(&idx.fmidx),
                black_box(h2r),
                black_box(refs),
                black_box(&record),
                black_box(SEED_LEN),
                black_box(false),
            )
            .unwrap()
        })
    });
    g.finish();

    // merge/merge_read_pairs ─────────────────────────────────────────────────────
    // Representative position maps (a handful of diagonals per mate).
    let fwd: HashMap<usize, LogProb> =
        (0..4).map(|k| (100 + k * 10, LogProb(-(k as f64) - 1.0))).collect();
    let rev: HashMap<usize, LogProb> =
        (0..4).map(|k| (250 + k * 10, LogProb(-(k as f64) - 1.0))).collect();
    let mut g = c.benchmark_group("merge");
    g.bench_function("merge_read_pairs", |b| {
        b.iter(|| merge_read_pairs(black_box(&fwd), black_box(&rev), black_box(READ_LEN)))
    });
    g.finish();

    // phase/process_read_pairs ───────────────────────────────────────────────────
    let pairs = synthetic_pairs(&ref0);
    let eps_1 = LogProb((1e-64f64).ln());
    let eps_2 = LogProb((1e-10f64).ln());
    let mut g = c.benchmark_group("phase");
    g.throughput(Throughput::Elements(N_PAIRS as u64));
    g.sample_size(10); // full parallel phase — keep wall-clock bounded
    g.bench_with_input(
        BenchmarkId::new("process_read_pairs", N_PAIRS),
        &pairs,
        |b, pairs| {
            b.iter(|| {
                process_read_pairs(
                    black_box(&idx.fmidx),
                    black_box(h2r),
                    black_box(refs),
                    black_box(pairs),
                    black_box(SEED_LEN),
                    black_box(eps_1),
                    black_box(eps_2),
                    None,
                )
                .unwrap()
            })
        },
    );
    g.finish();

    // build/sparse_from_alignments ──────────────────────────────────────────────
    // The DoK SparseArray build that sits between the parallel align and parallel
    // EM phases. Uses a larger read set (where the serial->parallel win matters).
    const BUILD_PAIRS: usize = 2000;
    let big_pairs: Vec<ReadPair> = (0..BUILD_PAIRS)
        .map(|i| {
            let max_off = REF_LEN - INSERT - READ_LEN;
            let o = (i * 7) % max_off;
            let r1 = ref0[o..o + READ_LEN].to_vec();
            let r2 = bio::alphabets::dna::revcomp(&ref0[o + INSERT..o + INSERT + READ_LEN]);
            ReadPair::new(&format!("r{i}"), rec("a", &r1), rec("b", &r2))
        })
        .collect();
    let (big_aligns, _) =
        process_read_pairs(&idx.fmidx, h2r, refs, &big_pairs, SEED_LEN, eps_1, eps_2, None).unwrap();
    let mut read_ids = HashMap::new();
    for (n, val) in big_aligns.iter().enumerate() {
        read_ids.insert(*val.key(), ReadIdx::new(n));
    }
    let mut g = c.benchmark_group("build");
    g.throughput(Throughput::Elements(BUILD_PAIRS as u64));
    g.bench_function("sparse_from_alignments", |b| {
        b.iter(|| build_sparse_array(black_box(&big_aligns), black_box(&read_ids)))
    });
    g.finish();
}

// pprof profiler: Output chosen via PPROF_OUT (svg default, "pb" = protobuf).
fn profiled() -> Criterion {
    let output = match std::env::var("PPROF_OUT").as_deref() {
        Ok("pb") => Output::Protobuf,
        _ => Output::Flamegraph(None),
    };
    Criterion::default().with_profiler(PProfProfiler::new(100, output))
}

criterion_group! {
    name = benches;
    config = profiled();
    targets = bench_alignment
}
criterion_main!(benches);
