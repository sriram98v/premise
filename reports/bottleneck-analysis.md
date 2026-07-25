# PREMISE runtime bottleneck analysis

Profiling of `premise-new` (`feat/seqid-ids`, `7458e6d`) on the synthetic benchmark samples,
against the real 4400-reference DB (`sequences.fasta`, 7.5 Mbp).

**Method.** `perf record` on the release binary (debug symbols on, fat LTO). Two flat profiles
at 499 Hz — `synthetic/isolate/Dataset-1` (383,538 samples) and `synthetic/mixed/Dataset-1`
(415,143 samples), 200k read pairs each, `-p 22 -i 200 -t 16`. Rankings agree closely between
the two, so the numbers below are the isolate profile unless noted.

`perf` works here despite `perf_event_paranoid = 2` — that setting blocks *kernel* profiling for
unprivileged users, not user-space sampling of your own processes. An earlier plan in this repo
assumed it was fully blocked and fell back to in-process `pprof`; that was unnecessary.

---

## Where the time goes

| origin | % runtime |
|---|--:|
| `occ/mod.rs` + `fm_index/mod.rs` | **56.8%** |
| `uint_macros.rs` (integer ops inlined into `rank`) | 14.2% |
| `fm_index/bidir.rs` (interval extension) | 7.2% |
| `malloc.c` / allocator | 6.5% |
| `range.rs` | 4.7% |
| `fm_index/smem.rs` | 2.6% |
| **premise's own scoring (`src/utils.rs`)** | **0.09%** |

Roughly **81% of runtime is inside haystackfm's FM-index search**. Premise's own alignment and
scoring code is statistical noise.

Two independent confirmations that the cost is *backward search*, not locate or rescoring:

1. **Seed-length sweep.** Raising `-p` from 22 to 70 collapses the number of qualifying SMEMs —
   and therefore locate calls and rescored diagonals — but runtime barely moves:

   | `-p` | 22 | 40 | 70 |
   |---|--:|--:|--:|
   | runtime | 12.87 s | 12.76 s | 12.71 s |

   A 1.2% spread across a 3× change in seed length. If locate or rescoring mattered, this would
   fall sharply.

2. **`sa_sample_rate = 1`.** Premise builds with a fully sampled suffix array, so `resolve_sa`
   is a direct lookup with no LF-walk. Locate is therefore close to free by construction.

**Consequence: further premise-side micro-optimisation is pointless.** The LUT, diagonal dedup,
`Cow` borrowing and `SeqId` work have driven premise's own cost to ~0.1%. Everything that remains
is in the index library or in configuration handed to it.

---

## Bottleneck 1 — IUPAC ambiguity expansion in SMEM search

**The largest single opportunity, and it is almost entirely wasted work.**

`premise` builds via `RefIndex::build_cpu`, which is `build_cpu_with::<IupacDna>`. The stored
alphabet drives `extend_multi_right` (`haystackfm/src/fm_index/smem.rs`):

```rust
fn extend_multi_right(ivs: &[BidirInterval], c: u8, rev: &FmIndex) -> Vec<BidirInterval> {
    let bases = (rev.alphabet_fns.compatible_fn)(c);   // IupacDna: 8 symbols for a plain base
    let mut result = Vec::new();
    for &base in bases {
        for iv in ivs {
            if let Some(ext) = iv.extend_right(base, rev) { result.push(ext); }
        }
    }
```

With `IupacDna` (`alphabet.rs:174`):

```
A => [A, N, R, W, M, D, H, V]      // 8
N => [A, C, G, T, N, R, Y, S, W, K, M, B, D, H, V]   // 15
```

So every extension step issues **8 interval extensions per interval**, each doing `rank()` calls.
`extend_multi_right` alone is 19.9% self time, and it drives much of the `occ` cost attributed to
the inlined `query_read`.

**But premise's indexed text contains only A, C, G, T, N.** `build_index_from_bytes` normalises
every reference base:

```rust
.map(|c| match c.to_ascii_uppercase() { 'A'|'C'|'G'|'T'|'N' => .., _ => 'N' })
```

and reads are encoded the same way (`encode_byte(b).unwrap_or(alphabet::N)`). Measured composition
of `sequences.fasta`: 7,495,851 bases, of which 345 are ambiguity codes — all collapsed to `N` at
build time. **R, W, M, D, H, V never occur in the text.** Six of the eight symbol searches per
extension step are guaranteed to return nothing, and still pay full `rank()` cost.

### Recommendations

| option | effect | cost / risk |
|---|---|---|
| **A. Build with `ExactDna`** (`build_cpu_with::<ExactDna>`) | 8 symbols → 1. Potentially the largest available win. | `N` in the reference stops matching, so reads overlapping those 345 positions lose seeds there. Same blast radius (0.0046%) as the ambiguity change already made — and arguably in the same direction. Requires haystackfm to expose a `BidirFmIndex::build_cpu_with` path premise can call (it exists). |
| **B. An ACGTN-only alphabet** upstream | `A => [A, N]`, `N => [A,C,G,T,N]`. 8 → 2 for ordinary bases, a 4× cut, while preserving current N semantics exactly. | New alphabet variant in haystackfm; a small, additive change. My preferred option — it keeps behaviour identical to today. |
| **C. Skip provably-absent symbols** upstream | Build a per-index symbol-presence bitmap and have `extend_multi_right` skip symbols with zero occurrences. | Fully automatic and behaviour-preserving for *any* input; helps every haystackfm user, not just premise. Slightly more invasive. |

Expected magnitude: `rank()` and its callers are >70% of runtime, and this cuts the number of
`rank()` calls per extension step by 4–8×. A **2–3× end-to-end speedup** is plausible. This should
be measured, not assumed — the interval set is not uniform and some extensions terminate early.

---

## Bottleneck 2 — Occ bitplane reconstruction (27.6%)

The single hottest source line in the whole profile is `occ/mod.rs:291`, inside `rank()`:

```rust
let mut mask = u64::MAX;
for p in 0..num_planes {
    let plane_val = self.word_at(base, p);
    mask &= if (lane >> p) & 1 == 1 { plane_val } else { !plane_val };   // line 291
}
```

With a 16-symbol alphabet this loop runs 4 iterations on **every** `rank()` call, reconstructing a
one-hot mask from bitplanes.

haystackfm already offers the alternative: `OccEncoding::OneHot` stores the one-hot vectors
directly and, per its own documentation, *"skips the reconstruction loop entirely"*. But
`OccEncoding::Bitplane` is `#[default]`, and premise builds with `..Default::default()`, so it
silently inherits the slower encoding.

### Recommendation

Set `occ_encoding: OccEncoding::OneHot` in `build_index_from_bytes`'s `RefIndexConfig`. One line.

**Caveat worth testing rather than assuming:** OneHot stores 16 `u64` per block instead of 4, so
the Occ table grows ~4×. Since `rank()` is memory-latency-bound, a larger table could give back
some of the win through cache misses. Measure both encodings end-to-end. This interacts with
Bottleneck 3 — freeing index memory there may make the larger Occ table affordable.

---

## Bottleneck 3 — Full suffix array buying nothing (index size)

Premise builds with `sa_sample_rate: 1`, i.e. a *completely* unsampled suffix array. The index is
**81.7 MB for 7.5 Mbp — 11.4 bytes per base**, of which roughly 60 MB is the forward and reverse
suffix arrays.

The point of `sa_sample_rate: 1` is O(1) locate with no LF-walk. But the seed-length sweep above
shows locate is a negligible share of runtime, so premise is paying ~60 MB of index — and the
attendant cache pressure on the *hot* Occ table — to accelerate something that isn't a bottleneck.

### Recommendation

Sweep `sa_sample_rate` over 1 / 8 / 16 / 32 and measure runtime and index size together. Expect
the index to shrink several-fold with little runtime cost, and possibly a runtime *gain* from
improved cache locality. This is the change most likely to make Bottleneck 2's OneHot encoding
a net win.

---

## Bottleneck 4 — Allocator churn (6.5%)

`malloc` / `realloc` / `free` / `RawVecInner::finish_grow` together account for ~6.5%.
`find_smems` returns `Vec<Mem>`, and each `Mem` owns a `Vec<(SeqId, u32)>` of positions;
`extend_multi_right` allocates a fresh `Vec<BidirInterval>` on **every extension step of every
read**. At 500k read pairs that is a very large number of short-lived small allocations.

### Recommendations

1. **Cheap and non-invasive:** switch the global allocator to `mimalloc` or `jemalloc`. A
   dependency plus three lines in `lib.rs`, no logic change. Typically recovers a good share of
   allocator overhead in this pattern.
2. **Upstream:** have `extend_multi_right` write into a caller-supplied scratch buffer instead of
   returning a fresh `Vec`, and offer a `find_smems_into(&mut Vec<Mem>)` variant so premise can
   reuse per-thread buffers across reads.

---

## Non-finding: the prefix lookup table does not apply

`FmIndexConfig::lookup_depth` defaults to 0, and a depth-k prefix table would ordinarily be an
obvious lever — it seeds backward search and skips the first k `rank()` calls.

**It does not help here.** The lookup table is only consulted in `query.rs:225`, on the
`count`/`locate` path. `smem.rs` never references it, so SMEM seeding — which is where all of
premise's time goes — cannot use it. Raising `lookup_depth` would grow the index for no gain.

Recorded so this is not re-proposed later.

---

## Summary and suggested order

| # | bottleneck | share | fix location | expected |
|---|---|--:|---|---|
| 1 | IUPAC expansion searching absent symbols | drives most of the >70% in `rank` | premise config, or upstream alphabet | **2–3×**, highest confidence of large gain |
| 2 | Occ bitplane reconstruction loop | 27.6% | premise config one-liner | large, but must verify cache trade |
| 3 | Full SA inflating the index | indirect | premise config | index shrink; enables #2 |
| 4 | Allocator churn | 6.5% | Cargo + upstream API | modest, cheap to try |
| — | premise's own scoring | 0.09% | — | nothing left |

Suggested sequence: **3 → 2 → 1 → 4.** Shrink the index first so the memory budget for OneHot
exists, then the encoding, then the alphabet change (largest but also the one that touches
matching semantics, so it wants a clean baseline and a correctness diff), then the allocator.

All four are configuration or upstream changes. **None of them are in premise's alignment code**,
which is now effectively free.
