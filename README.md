# PREMISE:  A probabilistic framework for source assignment of Illumina reads

[![Build](https://github.com/sriram98v/premise/actions/workflows/ci.yml/badge.svg)](https://github.com/sriram98v/premise/actions)
[![License](https://img.shields.io/github/license/sriram98v/premise)](https://github.com/sriram98v/premise/blob/main/LICENSE)
[![GitHub release](https://img.shields.io/github/v/release/sriram98v/premise)](https://github.com/sriram98v/premise/releases)
[![Stars](https://img.shields.io/github/stars/sriram98v/premise?style=social)](https://github.com/sriram98v/premise/stargazers)
[![Rust](https://img.shields.io/badge/built_with-Rust-orange?logo=rust)](https://www.rust-lang.org)
[![MSRV](https://img.shields.io/badge/MSRV-1.70-blue)](https://github.com/sriram98v/premise)
[![status: pre-release](https://img.shields.io/badge/status-pre--release-yellow)](https://github.com/sriram98v/premise)

**Authors:** Sriram Vijendran, Karin Dorman, Tavis Anderson, Oliver Eulenstein

## Introduction

PREMISE is an EM-based metagenomic classifier for paired-end Illumina short reads. Given a set of reference sequences and a paired-end FASTQ dataset, PREMISE:

1. Builds a run-length compressed FM-index over the reference sequences.
2. Seeds and extends alignments for each read pair against all references.
3. Runs a penalized Expectation-Maximization algorithm to estimate the posterior probability of each read originating from each reference (assignments) and the relative abundance of each reference in the sample (proportions).

The tool is implemented in Rust and ships with an optional browser-based GUI (served locally) for interactive use.

## Requirements

- [Rust / Cargo](https://rustup.rs/) ≥ 1.70
- Paired-end Illumina reads in FASTQ or FASTQ.gz format
- Reference sequences in FASTA format

## Installation

Install from crates.io / GitHub:

```bash
# via Cargo (recommended)
cargo install --git https://github.com/sriram98v/premise

# or build from source
git clone https://github.com/sriram98v/premise
cd premise
cargo install --path .
```

Or with [Nix](https://nixos.org/download/) (flakes enabled), which pins the whole toolchain:

```bash
# run without installing
nix run github:sriram98v/premise -- --help

# install into your profile
nix profile install github:sriram98v/premise

# or a dev shell with the pinned Rust toolchain (from a clone)
nix develop
```

## Usage

### Step 1 — Build an FM-index

```bash
premise build -s <reference.fasta>
```

Produces `<reference>.fmidx`, a run-length bidirectional FM-index ([haystackfm](https://github.com/sriram98v/haystackfm) `BidirRlFmIndex`). Its size grows with the number of BWT runs, i.e. with the variation in the reference collection rather than its length, so a collection of near-identical strains is far smaller than a plain FM-index (CAMI II strain madness, 1.48 Gbp at `--sa_sample_rate 16`: 1.8 GB vs 6.8 GB). `--sa_sample_rate R` (≥ 1) keeps one suffix-array sample in every R positions; `--lookup_depth K` sets the k-mer tables that seed the SMEM search. Indexes written by earlier versions (plain FM-index) are not readable; rebuild them.

### Step 2 — Classify reads (CLI)

```bash
premise query \
  -s <reference.fmidx>     \  # FM-index built in Step 1
  -1 <R1.fastq.gz> \
  -2 <R2.fastq.gz> \
  -m <min_seed_length>     \  # minimum MEM seed length (default 22)
  --eps_1 <float>          \  # alignment likelihood cutoff (default 0: keep every alignment)
  --eps_2 <float>          \  # minimum match log-probability (default 1e-18)
  --rho   <float>          \  # EM penalty weight ρ; > 0 turns the penalty on (default 0: plain EM, no pruning)
  --omega <float>          \  # EM penalty weight ω, used when ρ > 0 (default 1e-10)
  --iter  <int>            \  # EM iterations (default 100)
  --em_threshold <float>   \  # EM convergence: relative change in the penalized log-likelihood (default 1e-6)
  --matches                \  # also write <output>.matches, the vertex each read is assigned to
  --gamma <float>          \  # with --seq2taxid: descend to a vertex while its posterior exceeds gamma, in [0.5, 1) (default 0.5)
  --no-report              \  # with --matches: write only <output>.matches
  --compress               \  # zstd-compress the report (<output>.report.zst)
  --no-index-hash          \  # skip the SHA-256 of the index in the report (a full read of the index, run alongside the alignment)
  --init <split|uniform>   \  # initial conditionals: split-tie argmax counts (default) or 1/|children|
  -t <threads>             \  # default 2; 0 = all available cores
  --seq2taxid <map>        \  # optional: the taxid of every sequence, a vertex above it (see below)
  --nodes <nodes.dmp>      \  # optional, needs --seq2taxid: the taxonomy over those taxids (see below)
  -o <output_prefix>
```

`--seq2taxid` takes a kraken/centrifuge-style `seqid2taxid.map` (`seqid<whitespace>taxid`, extra columns ignored, blank and `#` lines skipped; a seqid is matched on its full FASTA header or its first token) and makes the *taxid* the unit the EM estimates. Alignment stays per sequence, P(read | taxid) = Σ P(read | sequence) over the taxid's sequences, and `.props`, `.posteriors` and `.matches` name taxids. Every taxid is an independent column. `--nodes` adds the lineages above them from an NCBI `nodes.dmp` (`taxid | parent | rank | …`) or plain `taxid parent [rank]` lines. Sequences absent from the map hang under the root, named by their header, and a warning counts them. Each read is then assigned by descending from the root to the child with the largest subtree posterior while that posterior exceeds `--gamma`. Without `--seq2taxid` there is no taxonomy to descend: each read goes to the sequence with its largest posterior (ties to the first in the index), the strain-level classification of earlier versions, and `--gamma` has no effect.

Outputs:
| File | Contents |
|------|----------|
| `<output>.report` | The run: taxonomy, every alignment log-likelihood of every read, proportions, EM trace and metadata (NEXUS, see [Output Format](#output-format)) |
| `<output>.matches` | Per-read assignments, only with `--matches` (TSV) |

Run `premise query -h` for the full option list. The report holds every alignment of every read, so on a large sample it is about the size of an alignment dump.

### Step 2b — Tables at any gamma, and EM refits, from the report

Everything but the read assignment is independent of `--gamma`, and the posteriors a read needs are recomputed from its stored log-likelihoods and the fitted proportions, so `analyze` turns a report into the legacy tables for any `gamma` without the index or the reads:

```bash
premise analyze -r <output>.report[.zst] -o <prefix> --matches --gamma 0.9 --posteriors --props --nodes
#   <prefix>.matches (--matches), .posteriors (--posteriors), .props, .nodes; --gamma defaults to 0.5
```

The report also holds every alignment likelihood the EM saw and the contracted tree it ran on, so the EM can be re-run with other penalty parameters, again without the index or the reads. Any EM option makes `analyze` a refit; whatever is not given keeps the report's value:

```bash
premise analyze -r <output>.report -o <prefix> --rho 8 --matches --props            # tables of the refit
premise analyze -r <output>.report -o <prefix> --no-penalty --write-report [--compress]  # + <prefix>.report[.zst]
```

| Option | Effect |
|--------|--------|
| `--rho R` | Penalty weight; R > 0 turns the penalty on, `--rho 0` turns it off |
| `--omega W`, `--iter N`, `--em_threshold T`, `--init split\|uniform` | As for `query` |
| `--no-penalty` | Plain EM (excludes `--rho`) |
| `--write-report` | Write the refit as a report (`SOURCE` names the input); needs an EM option |
| `-t N` | Threads for the EM |

The alignment parameters (`-m`, `--eps_1`, `--eps_2`) and the taxonomy (`--seq2taxid`, `--nodes`) are fixed by the report. `analyze` streams the report twice for a refit (once into the likelihood matrix, once for the tables), so a compressed report of any size works.

### Step 3 — Interactive GUI (optional)

```bash
premise server
```

Opens a browser UI at `http://localhost:8080` with drag-and-drop file upload, interactive results tables, pie chart, and EM convergence plot. Supports light/dark mode.

## Project Structure

```
premise/
├── src/
│   ├── main.rs          # binary entry point
│   ├── lib.rs           # the query pipeline (align -> EM -> tables) and the CLI
│   ├── align.rs         # SMEM seeding and ungapped extension: query_smems, query_read, query_fastq
│   ├── em.rs            # sparse likelihood matrix and the penalized EM iteration
│   ├── tree.rs          # the taxonomy (a `phylo` tree): construction, contraction, per-vertex M-step, pruning, read assignment
│   ├── group.rs         # seqid -> taxid grouping of index sequences
│   ├── io.rs            # index build/load, FASTQ loading, output tables
│   ├── server.rs        # local HTTP server for the browser GUI
│   ├── utils.rs         # shared types and the per-base likelihood model
│   └── templates/       # browser GUI (embedded at compile time)
├── tests/               # integration tests (CLI + local server)
├── benches/             # Criterion benchmarks
├── flake.nix            # the pinned PREMISE package and dev shell
├── Cargo.toml
└── README.md
```

## Output Format

### `.matches` (TSV)
| Column | Description |
|--------|-------------|
| `read_id` | Read identifier |
| `ref_id` | The vertex the read stopped at: a reference sequence ID, a taxid (with `--seq2taxid`), or `unclassified` |
| `posterior` | Posterior probability of the read's subtree at the assigned vertex |
| `forward/reverse position` | Position of the read's best-scoring sequence under that vertex |

### `.report` (NEXUS)

A `#NEXUS` text file (zstd-compressed by `query --compress`; `analyze` detects either). The standard `TAXA` and `TREES` blocks let tree viewers open it; the other blocks are private
and skipped by conforming readers. Labels are NEXUS-quoted (`'...'`, `''` for a quote). This is format version 3; PREMISE reads only its own version, so reports written by earlier versions must be regenerated.

| Block | Contents |
|-------|----------|
| `PREMISE` | Format version, tool version, timestamp, input file names (basenames only, never a path) and the SHA-256 of the index (`-` with `--no-index-hash`), `SOURCE` (the report a refit started from; absent for `query`), every parameter, thread count, read counts, taxonomy sizes, EM summary (iterations, pruned, collapsed, final penalized log-likelihood), load / align / EM seconds (`align_s` 0 for a refit) |
| `EM_TRACE` | Per EM iteration: the penalized log-likelihood and `(leaf, π_ℓ)` over the leaves of the tree at that iteration (entry 0 is the initial state, the last the final proportions) |
| `TAXA` | One label per vertex of the stored tree, in `NODES` order (the label of `NODES` row *i* is `TAXLABELS` *i*) |
| `TREES` | The tree topology as Newick over 1-based `TAXA` indices (row *i* of `NODES` is `i+1`), so viewers show the labels; rank and status are in `NODES` |
| `NODES` | Per vertex, in pre-order: `id parent rank status contigs count proportion conditional`. The tree is the contracted taxonomy down to the index sequences; `status` is `active` (in the final tree), `collapsed` (a final leaf pooling the sequences beneath it), `pooled` (beneath a collapsed vertex; `proportion` is its share `π̃_s`) or `pruned` |
| `REFMAP` | `seqid node`: the tip of every index sequence with reads |
| `PROPS` | `node proportion` for the final leaves with mass (what `.props` lists) |
| `ALIGNMENT_LIKELIHOODS` | Per read, `'read_id' 0` (no alignment) or `'read_id' n_rows n_lls` followed by: the `n_rows` sequence ids (the first absolute, then the difference to the previous one); the `n_lls` distinct alignment log-likelihoods in order of first appearance; each row's index into them (line omitted when `n_lls` is 1); the `n_rows` forward positions; the insert offsets `reverse − forward` (one value when every row has the same offset, else `n_rows`). Every finite alignment is kept, pruned references included |
| `PROFILE` | Write and total seconds, peak resident set size (KiB); written last |


### `.posteriors` (TSV, from `analyze --posteriors [--gamma G]`)
| Column | Description |
|--------|-------------|
| `read_id` | Read identifier |
| `ref_id` | An EM leaf the read aligned to: a reference sequence ID or the taxid of a collapsed vertex (its pooled sequences summed) |
| `posterior` | The read's posterior on that leaf |

One row per leaf with mass, sorted by label; the rows of a read sum to one (to the printed digits).

### `.props` (TSV, from `analyze --props`)
| Column | Description |
|--------|-------------|
| `ref_id` | Reference sequence ID, or the taxid of a collapsed vertex |
| `proportion` | Estimated relative abundance `π_ℓ` |

The leaves partition the mass, so the rows sum to one. The proportion of every taxid (the sum
over its subtree) is in `.nodes`.

### `.nodes` (TSV, from `analyze --nodes`)
| Column | Description |
|--------|-------------|
| `TaxID` | Vertex label (taxid, sequence header, or `root`) |
| `Parent` | Parent's label, `-` at the root |
| `Rank` | Rank from the nodes file (`sequence` for an index record, `root` for the synthetic root, `-` unknown) |
| `Contigs` | Sequences a leaf holds: 1 for a sequence, the pooled set for a collapsed vertex, 0 for an internal vertex |
| `Count` | Expected reads in the vertex's subtree |
| `Proportion` | Subtree proportion `π_v` (nested: a parent's includes its children's) |
| `Conditional` | Share of the parent's mass, `q_{v|parent}` |

Rows are in pre-order; pruned vertices are absent (the report's `NODES` block keeps them, with their status). Without `--seq2taxid` the table is the root plus one row per sequence.

## Citation

If you use PREMISE in your research, please cite:

> Vijendran S. *PREMISE: Probabilistic Read-level Expectation Maximization for Integrated Source Estimation.* (manuscript in preparation)

## License

See [LICENSE](LICENSE).
