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

---

# Follow-up #1 — SeqId ids + index-owned references

`haystackfm` 0.3.0 reports match locations as an integer `SeqId` and serves the bases
of any indexed sequence. That removes two lookups from the alignment inner loop:

- **Seeding.** `clean_mem_matches` used `header_to_ref[&j.0]` — a `String` hash + hashmap
  probe **per located occurrence** (per read × per SMEM × per hit), against a key the
  index had already heap-allocated to hand over. The id now comes straight off the match.
- **Rescoring.** `query_read` probed a `HashMap<RefIdx, Vec<u8>>` for the reference bases;
  it now slices them out of the index directly.

Plus a narrower key everywhere: `SeqId(u32)` replaced `RefIdx(usize)` through alignment,
EM and the sparse array, halving the reference half of every `(ReadIdx, RefIdx)` map key.

Baseline `opt234` (the `+#3 Cow` column above), same machine/data/seed.

| bench | opt234 | +SeqId | change | p | verdict |
|---|--:|--:|--:|--:|---|
| `phase/process_read_pairs/200` | 4.07 ms | **3.77 ms** | **−8.2%** | 0.00 | **improved** |
| `build/sparse_from_alignments` | 627 µs | **613 µs** | **−15.8%** | 0.00 | **improved** |
| `rescore/query_read` | 24.50 µs | **23.60 µs** | **−4.1%** | 0.01 | **improved** |
| `seed/clean_mem_matches` | 22.0 µs | 22.27 µs | +0.8% | 0.84 | no change |
| `scorer/compute_match_log_prob` | 795 ns | 776 ns | −1.9% | 0.32 | no change (control) |
| `scorer/error_prob` | 124 ns | 124 ns | +5.6% | 0.48 | no change (control) |
| `merge/merge_read_pairs` | 117 ns | 121 ns | +3.1% | 0.06 | no change (control) |

Flamegraph (`reports/flamegraph-seqid.svg` vs `reports/flamegraph-optimized.svg`):
`String`-related frames fall **37 → 9**.

## Interpretation

- **The phase-level −8.2% is the headline** and the number that matters: it is the whole
  parallel align stage, and it is the aggregate of the two removed lookups plus the
  narrower map keys.
- **`build/sparse_from_alignments` −15.8% was not predicted.** Nothing in that function
  changed except the key type: `(ReadIdx, SeqId)` is 4 bytes narrower than
  `(ReadIdx, RefIdx)`, which is enough to show up in a hash-heavy DoK build. A free win
  from the narrower id.
- **`seed/clean_mem_matches` is flat (p = 0.84), and that is a measurement artifact, not
  a null result.** This is the bench that should have improved most — it is where the
  `String` hash lived. But the synthetic fixture has only **4 references with 4-character
  headers** (`ref0`…`ref3`), so the removed hash was over a tiny key and the removed
  allocation was a tiny `String`. Real reference DBs carry thousands of long accession
  headers, where both costs scale up. **This bench understates the seeding win; do not
  read the flat result as "no gain on real data."** The phase-level number already
  reflects part of it.
- **The three controls stayed flat**, as intended. `compute_match_log_prob` now compares
  alphabet codes rather than ASCII bytes — same work per base, and the measurement
  confirms it.

## Cumulative, baseline → now

`phase/process_read_pairs/200`: **4.857 ms → 3.77 ms**, a **−22.4%** total reduction
across LUT+LTO, diagonal dedup, and SeqId.

## Correctness

Outputs are **byte-identical** across `feat/code-optim` → `feat/seqid-ids` for all four
TSVs (`.matches`, `.posteriors`, `.props`, `.aligns`) on `tests/fixtures/`, which is pure
uppercase ACGTN. 12 lib + 23 integration tests pass.

Two deliberate behavior changes apply only to references that are **not** pure uppercase
ACGTN, both corrections — scoring now uses the index's normalized text rather than raw
FASTA bytes:

1. Lowercase (soft-masked) reference bases no longer score as mismatches. The scorer
   compares bytes, so `'a'` against a read `'A'` was previously counted as a mismatch.
2. IUPAC ambiguity codes are ambiguous on both sides, instead of being indexed as `N` but
   scored literally.

A real-DB diff quantifying how many bases this touches is still outstanding.

---

# Real-data runtime: SeqId premise vs pre-refactor premise

The micro-benches above run against a synthetic 4-reference, 8 kb index. This section
measures the actual binary on a real workload, which tells a materially different story.

**Setup.** Reference DB `sequences.fasta` — 4400 influenza references, 7.5 Mb, mean header
131 chars (max 194). Reads: 400,000 real pairs from `SRR31013463` (109 bp). 8 threads,
`-p 20`. Binaries: `premise-old` = `feat/code-optim` (`d4a5d2d`), `premise-new` =
`feat/seqid-ids` (`7458e6d`). 3 alternating repetitions each.

| stage | old (mean) | new (mean) | change |
|---|--:|--:|--:|
| `query` (align + EM) | 89.71 s | **88.10 s** | **−1.8%** |
| `align` only | 88.20 s | **86.81 s** | **−1.6%** |
| index build | 1.14 s | 1.15 s | ~0% |

Run-to-run spread was tight — sd 0.08 s (old) and 0.17 s (new) on ~88 s, i.e. ~0.2% — so
the difference is well outside noise, just small. Scaling from a 10,000-pair run gives a
fixed overhead of ≈0 s and a per-pair cost of **220.4 µs → 216.8 µs**, confirming the gain
is genuinely per-pair rather than a startup artifact.

Index footprint: **85,827,340 → 85,646,978 bytes** (180 KB smaller — the reference bases
moved into the index instead of being stored a second time in premise's own blob). Peak
RSS during build rose 344 MB → 351 MB (+1.9%), the expected cost of retaining the text.

## Why this is far below the micro-bench's −8%

**A prediction in the previous section was wrong, and the real data corrects it.** That
section argued the synthetic bench *understates* the seeding win because its headers are
4 characters (`ref0`) against ~131 for real accessions, so the eliminated `String` hash and
allocation were unrealistically cheap. The header-length effect is real, but it is swamped
by an index-size effect running the other way:

- Per-pair cost is **220 µs on the real DB vs ~20 µs in the synthetic bench** — 11× more
  expensive per pair against an index ~1000× larger.
- That extra cost is FM-index backward search and SA resolution over an 85 MB structure —
  cache-miss-bound work that this change does not touch.
- The id bookkeeping that was removed is a *fixed* cost per occurrence. As search cost
  grows with index size, that bookkeeping shrinks as a fraction of the total.

So the synthetic bench **overstates** the end-to-end gain, because its trivially small index
makes bookkeeping a large share of a small total. Both effects are real; the index-size one
dominates.

**The remaining bottleneck is SMEM search itself, not reference-id handling.** Further work
on ids, map keys, or allocation around seeding has little left to recover on realistic
inputs.

## What the change is still worth

−1.6% per pair, reproducibly, plus a smaller index, one fewer side table, and the deletion
of `IOFMIndex`. It is free and correct. It is not the ~8% the micro-bench implied.

## Correctness on real data

193 of 405,784 alignment rows differ (**0.048%**), and the cause is confirmed rather than
assumed:

- **105 of 105 gained alignments (100%) fall on references carrying IUPAC ambiguity codes.**
  Top: `KT225475.1` (+57 rows, 42 ambiguity codes), `GU646028.1` (+37, 2 codes).
- 88 rows lost, half on ambiguity-bearing references — the knock-on from reads whose best
  hit moved.
- **6 reads went unclassified → classified.**
- Abundances (`.props`) list the same 4 references and agree to 5 significant figures;
  ordering differs because the map key type changed, so float summation order changed.

This is the predicted ambiguity-code correction, measured: previously an `R` in a reference
scored as a hard mismatch against every read base, often pushing an alignment below `eps_2`;
it is now skipped as ambiguous. Blast radius on this DB is **345 ambiguity codes in
7,495,851 bases (0.0046%)** and **zero lowercase bases**, so the lowercase fix has no effect
here.

## Criterion re-run on the final (crates.io) build

| bench | change vs opt234 | p |
|---|--:|--:|
| `build/sparse_from_alignments` | −16.2% | 0.00 |
| `phase/process_read_pairs/200` | −5.3% | 0.00 |
| `rescore/query_read` | −2.7% | 0.35 (ns) |
| `seed/clean_mem_matches` | +0.1% | 0.96 (ns) |
| `merge/merge_read_pairs` | +7.3% | 0.00 |

Two caveats. The phase figure was −8.2% on the earlier run and −5.3% here, so these
micro-benches carry several points of run-to-run variance — the real-data measurement above
is the more trustworthy signal. And `merge/merge_read_pairs` reports a *significant* 7.3%
regression despite being **completely unmodified** (it never touches reference ids); at
125 ns this is almost certainly a code-layout/LTO artifact rather than a real regression,
but it is recorded here rather than omitted.

---

# Cumulative: `origin/main` vs `feat/seqid-ids`

Every comparison above measures one increment. This one measures the whole effort end to end:
**`origin/main` (`35a3dcc`) vs `feat/seqid-ids` (`7458e6d`)** — LUT + fat LTO + diagonal dedup +
`Cow` borrowing + SeqId, together.

This matters because the 16-sample sweep's "old" binary was `feat/code-optim`, which *already*
contained the LUT and dedup work. That sweep therefore measured only the SeqId increment (−0.60%)
and badly understated the total.

The baseline is genuinely unoptimised: `origin/main` has `[profile.release] debug = true` and
nothing else — no LTO, default `codegen-units = 16` — and `haystackfm 0.1.0`. Its binary is
66.5 MB against 31.1 MB for the LTO'd build.

Same controls as the sweep: harness-matching filtered inputs and parameters (`-p 22`,
`-i 200`/`100`, same eps/rho/omega), `-t 16`, page cache warmed before each pair, run order
alternating, separate index per binary (the formats differ).

| sample | pairs | `origin/main` | `feat/seqid-ids` | change |
|---|--:|--:|--:|--:|
| real/isolate/SRR31013467 | 421,150 | 86.87 s | 75.98 s | **−12.54%** |
| real/mixed/SRR3360140 | 792,685 | 186.51 s | 164.58 s | **−11.76%** |
| synthetic/isolate/Dataset-1 | 500,000 | 155.05 s | 135.40 s | **−12.67%** |
| synthetic/mixed/Dataset-1 | 500,000 | 182.23 s | 147.52 s | **−19.05%** |
| **total** | | **610.66 s** | **523.48 s** | **−14.28%** |

**1.167× faster overall** — 87 s saved on 611 s of baseline work — and consistent across every
category, unlike the SeqId increment which was category-dependent.

## Where the gain actually came from

Splicing in each sample's `feat/code-optim` time from the 16-sample sweep splits the total into
its two phases:

| sample | `origin/main` → `feat/code-optim` | `feat/code-optim` → `feat/seqid-ids` |
|---|--:|--:|
| real/isolate/SRR31013467 | −10.4% | −2.4% |
| real/mixed/SRR3360140 | −11.5% | −0.3% |
| synthetic/isolate/Dataset-1 | −12.7% | +0.0% |
| synthetic/mixed/Dataset-1 | −17.8% | −1.5% |

**Roughly 90% of the improvement came from the first phase** — memoising `error_prob` into a
256-entry LUT, fat LTO, and diagonal dedup. The SeqId refactor contributed the remainder.

*Caveat:* the `feat/code-optim` column is taken from a different run session, so this split is
indicative rather than a controlled paired measurement. The `origin/main` vs `feat/seqid-ids`
columns above **are** controlled — paired, same session, alternating order, warm cache.

That split matches the profile. Before the LUT, `error_prob`'s `powf` ran once per base of every
diagonal of every read across four orientations, and scoring dominated. Afterwards `src/utils.rs`
measures **0.09%** of runtime and FM-index search accounts for ~81% — which is why the later,
seeding-side work had so much less headroom to recover.
