# Alignment benchmark baseline

**Baseline id:** `before` (Criterion `--save-baseline before`)

| Field | Value |
|---|---|
| CPU | Intel Core i7-10700K @ 3.80 GHz (8C/16T) |
| rustc | 1.96.1 (2026-06-26) |
| criterion | 0.5.1 · pprof 0.13 |
| profile | `bench` (opt-level 3, **no LTO**, codegen-units=16) — pre-optimization |
| commit | `7dd104e` (feat/code-optim) |
| date | 2026-07-24 |
| data | deterministic synthetic: 4 refs × 2000 bp w/ shared 300 bp conserved block; 150 bp reads; seed len 20 |

## Results

Median is the headline number (robust to outliers); mean±CI is the 95% bootstrap confidence interval Criterion reports.

| bench | what it measures | median | mean [95% CI] | std dev |
|---|---|--:|--:|--:|
| `scorer/compute_match_log_prob` | pure per-base ungapped rescore (150 bp) | 3.37 µs | 3.49 µs [3.41 µs – 3.57 µs] | 415.8 ns |
| `scorer/error_prob` | Phred→error-prob transform ×150 (the powf) | 1.96 µs | 2.03 µs [1.99 µs – 2.08 µs] | 237.0 ns |
| `seed/clean_mem_matches` | SMEM seeding + locate, one read | 26.75 µs | 29.29 µs [27.42 µs – 32.30 µs] | 13.13 µs |
| `rescore/query_read` | seed→project→rescore, one read | 37.96 µs | 40.27 µs [38.44 µs – 43.32 µs] | 13.42 µs |
| `merge/merge_read_pairs` | FR position-map cross-product (4×4) | 114.2 ns | 118.8 ns [115.8 ns – 122.4 ns] | 17.1 ns |
| `phase/process_read_pairs/200` | full rayon align phase, 200 pairs | 4.857 ms | 4.839 ms [4.742 ms – 4.933 ms] | 163.32 µs |

## Reading the numbers

- **`scorer/error_prob` ≈ 2.03 µs** for 150 calls → **~14 ns/call**, almost all of it the `10f64.powf(...)`. This alone is **~58% of `compute_match_log_prob`** (3.49 µs) — the single clearest target for the LUT win.
- **`compute_match_log_prob` ≈ 3.49 µs** for a 150 bp read → **~23.2 ns/base**. In the real pipeline this runs `read_len × #diagonals × 4 orientations × #pairs` times.
- **`seed/clean_mem_matches` ≈ 29.29 µs** and **`rescore/query_read` ≈ 40.27 µs** (query_read = seed + rescore): seeding/locate is the larger share of a single read's alignment, matching the flamegraph (`find_smems`/`locate` frames dominate the phase). This is what **LTO** should help by inlining haystackfm FM primitives.
- **`phase/process_read_pairs/200` ≈ 4.839 ms** for 200 pairs across 16 threads → **~24.2 µs/pair** wall-clock (parallel).
- **`merge/merge_read_pairs` ≈ 118.8 ns**: negligible, as predicted (tiny diagonal counts).

## Artifacts (industry-standard formats)

- **Criterion HTML** (interactive, per-bench PDF/violin + regression plots): open `target/criterion/report/index.html`.
- **Criterion JSON** (machine-readable): `target/criterion/<bench>/new/estimates.json` (+ `sample.json`, `benchmark.json`).
- **Flamegraph** (pprof, in-process): `reports/flamegraph-baseline.svg` — open in a browser; frame width = share of samples. `find_smems`/`locate` and `compute_match_log_prob`→`error_prob`→`powf` are visible.
- **pprof protobuf**: `reports/profile-baseline.pb` — `pprof -http=: reports/profile-baseline.pb`.

See `reports/README.md` for what each format is and how to read it.
