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
