# Alignment benchmarks — optimized vs baseline

**Comparison:** optimized run vs Criterion baseline `before` (`cargo bench --bench alignment -- --baseline before`).
Same machine / data / seed as the baseline (see `baseline.md`).

**Changes applied** (commit `e4d9b1d`):
1. **LUT** — `error_prob` serves `10^(-((q-33)/10))` from a 256-entry `LazyLock` table instead of a per-base `powf` (`src/utils.rs`).
2. **LTO** — `lto = "fat"` + `codegen-units = 1` on `[profile.release]` and `[profile.bench]` (`Cargo.toml`), so haystackfm's FM primitives inline into the seeding loop.

## Results (median; % change is Criterion's, with its t-test verdict)

| bench | driver | before | after | change | significance |
|---|---|--:|--:|--:|---|
| `scorer/error_prob` | LUT | 1.96 µs | 121 ns | **-93.9%** | p<0.05 ✓ |
| `scorer/compute_match_log_prob` | LUT | 3.37 µs | 745 ns | **-77.9%** | p<0.05 ✓ |
| `rescore/query_read` | LUT+LTO | 37.96 µs | 23.90 µs | **-38.9%** | p<0.05 ✓ |
| `seed/clean_mem_matches` | LTO | 26.75 µs | 20.92 µs | **-26.2%** | p<0.05 ✓ |
| `phase/process_read_pairs/200` | LUT+LTO | 4.857 ms | 4.241 ms | **-11.1%** | p<0.05 ✓ |
| `merge/merge_read_pairs` | — | 114 ns | 114 ns | **-2.0%** | p=0.27 (noise) |

## Attribution

- **LUT → scorer.** `error_prob` dropped **−94%** (2.03 µs → 124 ns for 150 calls); `compute_match_log_prob` **−78%** (3.49 µs → 774 ns). The `powf` is gone: it appeared in **84** flamegraph frames at baseline and **0** after (`error_prob`: 42 → 0). This is the single biggest per-call win.
- **LTO → seeding + phase.** `clean_mem_matches` **−26%** and `query_read` **−39%** (the latter also gets the LUT via its rescore). Cross-crate inlining of haystackfm's `find_smems`/`locate`/`lf_mapping` into PREMISE is the driver. The full parallel phase improved **−11%** wall-clock (LUT + LTO combined), i.e. ~24.2 → ~21.5 µs/pair.
- **merge unchanged** (−2%, p=0.27): as predicted — no `powf`, trivially small diagonal counts, nothing for LTO to inline meaningfully. A clean negative control that the harness detects *no* change where none exists.

## Visual confirmation

Compare `flamegraph-baseline.svg` vs `flamegraph-optimized.svg` (Ctrl-F `powf`): the
`compute_match_log_prob → error_prob → powf` stack visible at baseline is absent
after the LUT. `find_smems`/`locate` frames shrink under LTO.

pprof diff of the two profiles:
```bash
pprof -diff_base reports/profile-baseline.pb reports/profile-optimized.pb
```

## Caveats

- Synthetic in-memory data (4 refs, 300 bp shared block). Real metagenomic DBs have far higher seed multiplicity, so the **LTO** locate gains should be *larger* in production, while the **LUT** gain is input-independent (always removes the per-base `powf`).
- The phase bench runs 200 pairs over 16 threads; absolute per-pair time is parallel wall-clock, not single-core cost.
- EM stages not measured (future work).

---

# Follow-up edits: diagonal dedup (#2), borrow forward reads (#3), DoK build (#4)

Same machine/data/seed. Columns are **cumulative** (each adds to the one on its
left). Medians. `≈` = statistically unchanged (edit doesn't touch that path).

| bench | before | LUT+LTO | +#2 dedup | +#3 Cow | verdict |
|---|--:|--:|--:|--:|---|
| `scorer/error_prob` | 1.96 µs | 121 ns | ≈ 124 ns | ≈ 124 ns | untouched |
| `scorer/compute_match_log_prob` | 3.37 µs | 745 ns | ≈ 795 ns | ≈ 795 ns | untouched |
| `rescore/query_read` | 37.96 µs | 23.90 µs | 24.06 µs | 24.50 µs | flat (noise) |
| `seed/clean_mem_matches` | 26.75 µs | 20.92 µs | ≈ 22.0 µs | ≈ 22.0 µs | untouched |
| `merge/merge_read_pairs` | 114 ns | 114 ns | ≈ 117 ns | ≈ 117 ns | untouched |
| `phase/process_read_pairs/200` | 4.857 ms | 4.241 ms | **4.10 ms** | **4.07 ms** | **−3.3% then −0.7%** |

### #4 — parallel DoK build (measured separately; not cumulative)

New bench `build/sparse_from_alignments` (2000 reads):

| variant | median | outcome |
|---|--:|---|
| serial (`iter`) | **627 µs** | kept |
| parallel (`par_iter` over shared `DashMap`) | 729 µs | **reverted — +16% regression** |

## Interpretation

- **#2 (dedup diagonals)** is the only follow-up with a real signal: **phase −3.3%**
  (4.24 → 4.10 ms). It fires only when several SMEMs on one reference project to
  the same diagonal; the synthetic data has low multiplicity, so this is a floor —
  on real metagenomic DBs with high seed multiplicity the gain should be larger.
  `query_read` is flat because post-LUT the per-read rescore is already cheap.
- **#3 (Cow / borrow forward reads)** is within noise here (saves two ~150-byte
  `Vec` clones on the forward orientations). Kept anyway: strictly fewer
  allocations, lower allocator pressure at scale, cleaner code — no downside.
- **#4 (parallelize DoK build)** was implemented and benchmarked, and it **regressed
  ~16%**: rayon-over-a-shared-`DashMap` trades a cheap serial insert for shard-lock
  contention at this read count. Kept **serial**. The extraction + bench remain so a
  future lock-free design (thread-local buffers + merge) can be measured directly.

**Bottom line:** the large wins are still LUT (input-independent) + LTO. Of the
follow-ups only diagonal dedup pays off measurably here, and it should pay off more
on real data. #4's takeaway is a negative result worth keeping: naive shared-map
parallelism loses to the serial build at these sizes.
