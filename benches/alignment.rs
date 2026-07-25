use bio::io::fastq;
use bio::stats::LogProb;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use haystackfm::alphabet;
use haystackfm::alphabet::encode_byte;
#[cfg(unix)]
use pprof::criterion::{Output, PProfProfiler};
use premise::{
    build_index_from_bytes, build_sparse_array, clean_mem_matches, load_index, merge_read_pairs,
    process_read_pairs, query_read, ReadIdx, ReadPair, RefIndex,
};
use std::collections::HashMap;


const READ_LEN: usize = 150;
const SEED_LEN: usize = 20;
const N_REFS: usize = 4;
const REF_LEN: usize = 2000;
const N_PAIRS: usize = 200;
const INSERT: usize = 200;

fn lcg(state: &mut u64) -> u64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    *state
}

fn random_dna(len: usize, state: &mut u64) -> Vec<u8> {
    const B: [u8; 4] = [b'A', b'C', b'G', b'T'];
    (0..len)
        .map(|_| B[(lcg(state) >> 33) as usize & 3])
        .collect()
}

fn synthetic_fasta() -> (Vec<u8>, Vec<Vec<u8>>) {
    let mut state = 0x9E3779B97F4A7C15u64; // fixed seed
    let conserved = random_dna(300, &mut state);
    let mut fasta = Vec::new();
    let mut refs = Vec::with_capacity(N_REFS);
    for i in 0..N_REFS {
        let mut seq = random_dna(REF_LEN, &mut state);

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

fn main_index() -> (RefIndex, Vec<Vec<u8>>) {
    let (fasta, refs) = synthetic_fasta();
    let (bytes, _) = build_index_from_bytes(&fasta).expect("build index");
    (load_index(&bytes, "<bench>").expect("load index"), refs)
}

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

fn bench_alignment(c: &mut Criterion) {
    let (idx, refs) = main_index();
    let ref0 = refs[0].clone();

    let mut read = ref0[100..100 + READ_LEN].to_vec();
    read[40] ^= 0b0000_0110; // flip a couple of bases so it isn't a perfect match
    read[110] ^= 0b0000_0110;
    let qual = vec![b'I'; READ_LEN];
    let ref_seg: Vec<u8> = ref0[100..100 + READ_LEN]
        .iter()
        .map(|&b| encode_byte(b).unwrap_or(alphabet::N))
        .collect();
    let q_enc: Vec<u8> = read
        .iter()
        .map(|&b| encode_byte(b).unwrap_or(alphabet::N))
        .collect();
    let record = rec("read", &read);

    // scorer/compute_match_log_prob ───────────────────────────────────────────
    let mut g = c.benchmark_group("scorer");
    g.throughput(Throughput::Bytes(READ_LEN as u64));
    g.bench_function("compute_match_log_prob", |b| {
        b.iter(|| {
            premise::utils::compute_match_log_prob(
                black_box(&q_enc),
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

    let mut g = c.benchmark_group("seed");
    g.bench_function("clean_mem_matches", |b| {
        b.iter(|| clean_mem_matches(black_box(&idx), black_box(&q_enc), black_box(SEED_LEN)))
    });
    g.finish();

    let mut g = c.benchmark_group("rescore");
    g.bench_function("query_read", |b| {
        b.iter(|| {
            query_read(
                black_box(&idx),
                black_box(&record),
                black_box(SEED_LEN),
                black_box(false),
            )
            .unwrap()
        })
    });
    g.finish();

    let fwd: HashMap<usize, LogProb> = (0..4)
        .map(|k| (100 + k * 10, LogProb(-(k as f64) - 1.0)))
        .collect();
    let rev: HashMap<usize, LogProb> = (0..4)
        .map(|k| (250 + k * 10, LogProb(-(k as f64) - 1.0)))
        .collect();
    let mut g = c.benchmark_group("merge");
    g.bench_function("merge_read_pairs", |b| {
        b.iter(|| merge_read_pairs(black_box(&fwd), black_box(&rev), black_box(READ_LEN)))
    });
    g.finish();

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
                    black_box(&idx),
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
        process_read_pairs(&idx, &big_pairs, SEED_LEN, eps_1, eps_2, None).unwrap();
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

#[cfg(unix)]
fn profiled() -> Criterion {
    let output = match std::env::var("PPROF_OUT").as_deref() {
        Ok("pb") => Output::Protobuf,
        _ => Output::Flamegraph(None),
    };
    Criterion::default().with_profiler(PProfProfiler::new(100, output))
}

#[cfg(not(unix))]
fn profiled() -> Criterion {
    Criterion::default()
}

criterion_group! {
    name = benches;
    config = profiled();
    targets = bench_alignment
}
criterion_main!(benches);
