/// Benchmark: genedex FM-index vs webgpu-fmidx (CPU and GPU) on a synthetic dataset.
///
/// Run:
///   cargo bench --bench fmidx_comparison
///
/// For GPU path:
///   cargo bench --bench fmidx_comparison -- gpu
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};

// ── genedex backend ──────────────────────────────────────────────────────────
use genedex::text_with_rank_support::{Block64, FlatTextWithRankSupport};
use genedex::{alphabet as genedex_alphabet, FmIndexConfig, FmIndexFlat64};

// ── webgpu-fmidx backend ─────────────────────────────────────────────────────
use webgpu_fmidx::alphabet::encode_char;
use webgpu_fmidx::{DnaSequence, FmIndex as WgpuFmIndex, FmIndexConfig as WgpuConfig};

// ── GPU path (async, blocked with pollster) ──────────────────────────────────
#[cfg(feature = "gpu")]
use webgpu_fmidx::FmIndexConfig as WgpuGpuConfig;

// ── test data ────────────────────────────────────────────────────────────────
/// A simple synthetic FASTA with 3 reference sequences of 500 bp each.
const FASTA: &[(&str, &str)] = &[
    ("ref_A", "AGCTAGCTAGCTAGCTTACGATCGATCGAATCGAATCGATCGATCGATCGAGCTAGCTAGCTAGCTTACGATCGATCGAATCGAATCGATCGATCGATCGAGCTAGCTAGCTAGCTTACGATCGATCGAATCGAATCGATCGATCGATCGAGCTAGCTAGCTAGCTTACGATCGATCGAATCGAATCGATCGATCGATCGAGCTAGCTAGCTAGCTTACGATCGATCGAATCGAATCGATCGATCGATCGAGCTAGCTAGCTAGCTTACGATCGATCGAATCGAATCGATCGATCGATCGAGCTAGCTAGCTAGCTTACGATCGATCGAATCGAATCGATCGATCGATCGAGCTAGCTAGCTAGCTTACGATCGATCGAATCGAATCGATCGATCGATCGAGCTAGCTAGCTAGCTTACGATCGATCGAATCGAATCGATCGATCGATCGAGCTAGCTAGCTAGCTTACGATCGATCGAATCGAATCGATCGATCGATCG"),
    ("ref_B", "TGCATGCATGCATGCAGTCAGTCAGTCAGTCAAATCGATCGATCGATCGTGCATGCATGCATGCAGTCAGTCAGTCAGTCAAATCGATCGATCGATCGTGCATGCATGCATGCAGTCAGTCAGTCAGTCAAATCGATCGATCGATCGTGCATGCATGCATGCAGTCAGTCAGTCAGTCAAATCGATCGATCGATCGTGCATGCATGCATGCAGTCAGTCAGTCAGTCAAATCGATCGATCGATCGTGCATGCATGCATGCAGTCAGTCAGTCAGTCAAATCGATCGATCGATCGTGCATGCATGCATGCAGTCAGTCAGTCAGTCAAATCGATCGATCGATCGTGCATGCATGCATGCAGTCAGTCAGTCAGTCAAATCGATCGATCGATCGTGCATGCATGCATGCAGTCAGTCAGTCAGTCAAATCGATCGATCGATCGTGCATGCATGCATGCAGTCAGTCAGTCAGTCAAATCGATCGATCGATCG"),
    ("ref_C", "CCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCC"),
];

/// Query patterns (16-mers from ref_A).
const PATTERNS: &[&str] = &[
    "AGCTAGCTAGCTAGCT",
    "TACGATCGATCGAATC",
    "GAATCGATCGATCGAT",
    "CGAGCTAGCTAGCTAG",
];

// ─── genedex helpers ──────────────────────────────────────────────────────────

fn build_genedex_index() -> FmIndexFlat64<i64> {
    let texts: Vec<&str> = FASTA.iter().map(|(_, seq)| *seq).collect();
    let dna_alphabet = genedex_alphabet::ascii_dna_iupac_as_dna_with_n();
    FmIndexConfig::<i64, FlatTextWithRankSupport<i64, Block64>>::new()
        .suffix_array_sampling_rate(1)
        .lookup_table_depth(13)
        .construct_index(texts, dna_alphabet)
}

fn query_genedex(fmidx: &FmIndexFlat64<i64>, patterns: &[&str]) -> usize {
    let mut total = 0usize;
    for pattern in patterns {
        let encoded: Vec<_> = fmidx
            .locate_many(std::iter::once(pattern.as_bytes()))
            .flatten()
            .collect();
        total += encoded.len();
    }
    total
}

// ─── webgpu-fmidx CPU helpers ─────────────────────────────────────────────────

fn build_wgpu_cpu_index() -> WgpuFmIndex {
    let sequences: Vec<DnaSequence> = FASTA
        .iter()
        .map(|(header, seq)| DnaSequence::from_str_with_header(seq, header).unwrap())
        .collect();
    let config = WgpuConfig {
        sa_sample_rate: 1,
        use_gpu: false,
    };
    WgpuFmIndex::build_cpu(&sequences, &config).unwrap()
}

fn encode_pattern(s: &str) -> Vec<u8> {
    s.chars().filter_map(|c| encode_char(c)).collect()
}

fn query_wgpu(fmidx: &WgpuFmIndex, patterns: &[&str]) -> usize {
    let mut total = 0usize;
    for pattern in patterns {
        total += fmidx.locate(&encode_pattern(pattern)).len();
    }
    total
}

// ─── webgpu-fmidx GPU helpers ─────────────────────────────────────────────────

#[cfg(feature = "gpu")]
fn build_wgpu_gpu_index() -> WgpuFmIndex {
    let sequences: Vec<DnaSequence> = FASTA
        .iter()
        .map(|(header, seq)| DnaSequence::from_str_with_header(seq, header).unwrap())
        .collect();
    let config = WgpuGpuConfig {
        sa_sample_rate: 1,
        use_gpu: true,
    };
    pollster::block_on(WgpuFmIndex::build(&sequences, &config)).unwrap()
}

// ─── benchmarks ───────────────────────────────────────────────────────────────

fn bench_index_construction(c: &mut Criterion) {
    let mut group = c.benchmark_group("index_construction");

    group.bench_function(BenchmarkId::new("genedex", "3x500bp"), |b| {
        b.iter(|| build_genedex_index())
    });

    group.bench_function(BenchmarkId::new("webgpu_cpu", "3x500bp"), |b| {
        b.iter(|| build_wgpu_cpu_index())
    });

    #[cfg(feature = "gpu")]
    group.bench_function(BenchmarkId::new("webgpu_gpu", "3x500bp"), |b| {
        b.iter(|| build_wgpu_gpu_index())
    });

    group.finish();
}

fn bench_query(c: &mut Criterion) {
    let genedex_idx = build_genedex_index();
    let wgpu_cpu_idx = build_wgpu_cpu_index();

    #[cfg(feature = "gpu")]
    let wgpu_gpu_idx = build_wgpu_gpu_index();

    let mut group = c.benchmark_group("query_4x16mer");

    group.bench_function(BenchmarkId::new("genedex", "locate"), |b| {
        b.iter(|| query_genedex(&genedex_idx, PATTERNS))
    });

    group.bench_function(BenchmarkId::new("webgpu_cpu", "locate"), |b| {
        b.iter(|| query_wgpu(&wgpu_cpu_idx, PATTERNS))
    });

    #[cfg(feature = "gpu")]
    group.bench_function(BenchmarkId::new("webgpu_gpu", "locate"), |b| {
        b.iter(|| query_wgpu(&wgpu_gpu_idx, PATTERNS))
    });

    group.finish();
}

criterion_group!{
    name = benches;
    config = Criterion::default().sample_size(100);
    targets = bench_index_construction, bench_query
}
criterion_main!(benches);
