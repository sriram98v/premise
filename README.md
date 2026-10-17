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
  --em_stop <none|loglik|props> \  # stopping criterion. none (default): run all --iter iterations; loglik: relative change in the penalized log-likelihood; props: Euclidean distance between consecutive leaf proportion vectors in inverse-softmax (centred log-ratio) space, over leaves with at least 1/N of the mass
  --em_threshold <float>   \  # with --em_stop: stop once the change falls to this, checked after 50 iterations (default 1e-6)
  --matches                \  # also write <output>.matches, the vertex each read is assigned to at --gamma
  --gamma <float>          \  # descend to a vertex while its posterior exceeds gamma (default 0.6)
  --no-report              \  # skip the report: write only <output>.props
  --compress               \  # zstd-compress the report (<output>.report.zst)
  --no-index-hash          \  # skip the SHA-256 of the index in the report (a full read of the index, run alongside the alignment)
  --init <split|uniform>   \  # initial conditionals: split-tie argmax counts (default) or 1/|children|
  -t <threads>             \  # default 2; 0 = all available cores
  --seq2taxid <map>        \  # optional: the taxid of every sequence, a vertex above it (see below)
  --nodes <nodes.dmp>      \  # optional, needs --seq2taxid: the taxonomy over those taxids (see below)
  -o <output_prefix>
```

`--seq2taxid` takes a kraken/centrifuge-style `seqid2taxid.map` (`seqid<whitespace>taxid`, extra columns ignored, blank and `#` lines skipped; a seqid is matched on its full FASTA header or its first token) and makes the *taxid* the unit the EM estimates. Alignment stays per sequence, P(read | taxid) = Σ P(read | sequence) over the taxid's sequences, and `.props` and `.matches` name taxids. Every taxid is an independent column. `--nodes` adds the lineages above them from an NCBI `nodes.dmp` (`taxid | parent | rank | …`) or plain `taxid parent [rank]` lines. Sequences absent from the map hang under the root, named by their header, and a warning counts them. Each read in `.matches` is assigned by descending from the root to the child with the largest subtree posterior while that posterior exceeds `--gamma` (default 0.6): the lowest vertex the read is confidently placed under. Without `--seq2taxid` the tree is the root over the sequences, so a read is assigned a sequence only when its posterior on it exceeds `--gamma`, and is `unclassified` otherwise (e.g. a read split evenly between near-identical references).

Outputs:
| File | Contents |
|------|----------|
| `<output>.report` | The run: taxonomy with the fitted proportions, every alignment log-likelihood of every read, EM trace and metadata (NEXUS, see [Output Format](#output-format)) |
| `<output>.props` | The fitted taxonomy with the proportion of every vertex that has any (NEXUS) |
| `<output>.matches` | Per-read assignments at `--gamma`, only with `--matches` (TSV) |

Run `premise query -h` for the full option list. The report holds every alignment of every read, so on a large sample it is about the size of an alignment dump.

### Step 2b — Tables at any gamma, and EM refits, from the report

Only the read assignment depends on `--gamma`, and the posteriors a read needs follow from its stored log-likelihoods and the fitted proportions, so `analyze` turns a report into `.matches` at any `gamma`, and `.props`, without the index or the reads:

```bash
premise analyze -r <output>.report[.zst] -o <prefix> --matches --gamma 0.9 --props
#   <prefix>.matches (--matches), <prefix>.props (--props); --gamma defaults to 0.6
```

The report also holds every alignment likelihood the EM saw and the contracted tree it ran on, so the EM can be re-run with other penalty parameters, again without the index or the reads. Any EM option makes `analyze` a refit; whatever is not given keeps the report's value:

```bash
premise analyze -r <output>.report -o <prefix> --rho 8 --matches --props            # tables of the refit
premise analyze -r <output>.report -o <prefix> --no-penalty --write-report [--compress]  # + <prefix>.report[.zst]
```

| Option | Effect |
|--------|--------|
| `--rho R` | Penalty weight; R > 0 turns the penalty on, `--rho 0` turns it off |
| `--omega W`, `--iter N`, `--em_threshold T`, `--em_stop none\|loglik\|props`, `--init split\|uniform` | As for `query` |
| `--no-penalty` | Plain EM (excludes `--rho`) |
| `--write-report` | Write the refit as a report (`SOURCE` names the input); needs an EM option |
| `--em-props` | Write `<prefix>.em_props.tsv`: the EM's trace, one row per iteration (`0` is the initial state, the last row the final fit) and one column per reference (REFMAP order, named by sequence) holding its proportion π̃. Runs the EM, with the report's parameters unless overridden, and the other tables of the call come from that run |
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

### `.matches` (TSV, with `--matches`)

The first line is a `#` comment naming the EM that produced the proportions and the gamma, e.g.

```
# PREMISE 0.6.0 EM: init=split iter=100 rho=0e0 omega=1e-10 em_threshold=1e-6 em_stop=none penalty=0; 100 iterations, final penalized log-likelihood -1.2e4; gamma=0.6
```

(a refit adds `, refit of <report>` before the gamma). The column header follows:

| Column | Description |
|--------|-------------|
| `read_id` | Read identifier |
| `ref_id` | The deepest vertex whose subtree posterior exceeds `--gamma` along the most probable path from the root: a taxid (with `--seq2taxid`), a reference sequence ID, or `unclassified` |
| `posterior` | Posterior probability of the read's subtree at the assigned vertex |
| `forward/reverse position` | Position of the read's best-scoring sequence under that vertex |

### `.props` (NEXUS)

The fitted taxonomy restricted to the vertices of the final tree with a non-zero proportion, which form a tree under the root:

```
#NEXUS
[PREMISE proportions v1]
[PREMISE 0.6.0 EM: init=split iter=100 rho=0e0 omega=1e-10 em_threshold=1e-6 em_stop=none penalty=0; 100 iterations, final penalized log-likelihood -1.2e4]

BEGIN TAXA;
	DIMENSIONS NTAX=5;
	TAXLABELS
		root
		...
	;
END;

BEGIN TREES;
	TREE proportions = [&R] ((3,(5)4)2)1;
END;

BEGIN PROPORTIONS;
	DIMENSIONS N=5;
	[label parent rank proportion; ...]
	MATRIX
	root - root 1e0
	561 root genus 1e0
	562 561 species 4e-1
	...
	;
END;
```

| Block | Contents |
|-------|----------|
| `TAXA` | One label per vertex (taxid, sequence header, or `root`), NEXUS-quoted |
| `TREES` | The topology as Newick over 1-based `TAXA` indices |
| `PROPORTIONS` | Per vertex, in pre-order: `label parent rank proportion`, with `parent` the parent's label (`-` at the root) and `proportion` the end-of-EM `π_v`, the mass of the subtree under the vertex (1 at the root). The leaves partition the mass, so their proportions sum to one |

Without `--seq2taxid` the tree is the root plus one row per sequence with mass.

### `.report` (NEXUS)

A `#NEXUS` text file (zstd-compressed by `query --compress`; `analyze` detects either). The standard `TAXA` and `TREES` blocks let tree viewers open it; the other blocks are private
and skipped by conforming readers. Labels are NEXUS-quoted (`'...'`, `''` for a quote). This is format version 4; PREMISE reads only its own version, so reports written by earlier versions must be regenerated.

| Block | Contents |
|-------|----------|
| `PREMISE` | Format version, tool version, timestamp, input file names (basenames only, never a path) and the SHA-256 of the index (`-` with `--no-index-hash`), `SOURCE` (the report a refit started from; absent for `query`), every parameter, thread count, read counts, taxonomy sizes, EM summary (iterations, pruned, collapsed, final penalized log-likelihood), load / align / EM seconds (`align_s` 0 for a refit) |
| `EM_TRACE` | Per EM iteration: the penalized log-likelihood and `(leaf, π_ℓ)` over the leaves of the tree at that iteration (entry 0 is the initial state, the last the final proportions) |
| `TAXA` | One label per vertex of the stored tree, in `NODES` order (the label of `NODES` row *i* is `TAXLABELS` *i*) |
| `TREES` | The tree topology as Newick over 1-based `TAXA` indices (row *i* of `NODES` is `i+1`), so viewers show the labels; rank and status are in `NODES` |
| `NODES` | Per vertex, in pre-order: `id parent rank status contigs count proportion conditional`. The tree is the contracted taxonomy down to the index sequences; `status` is `active` (in the final tree), `collapsed` (a final leaf pooling the sequences beneath it), `pooled` (beneath a collapsed vertex; `proportion` is its share `π̃_s`) or `pruned` |
| `REFMAP` | `seqid node`: the tip of every index sequence with reads |
| `ALIGNMENT_LIKELIHOODS` | Per read, `'read_id' 0` (no alignment) or `'read_id' n_rows n_lls` followed by: the `n_rows` sequence ids (the first absolute, then the difference to the previous one); the `n_lls` distinct alignment log-likelihoods in order of first appearance; each row's index into them (line omitted when `n_lls` is 1); the `n_rows` forward positions; the insert offsets `reverse − forward` (one value when every row has the same offset, else `n_rows`). Every finite alignment is kept, pruned references included |
| `PROFILE` | Write and total seconds, peak resident set size (KiB); written last |

The report holds only what the EM cannot regenerate. The per-vertex counts, contigs, conditionals and proportions of the fit are in the `NODES` block; the final leaves' proportions are those of its active or collapsed vertices with no such vertex below them. Each read's posterior on a reference, `π̃_s · L(read | s)` normalised over the read's alignments, follows from `ALIGNMENT_LIKELIHOODS` and those proportions, which is how `analyze` writes `.matches` at any gamma. Reports written by earlier 0.6.0 builds, which also held `PROPS` and `READ_POSTERIORS` blocks, are still read; the extra blocks are skipped.

## Citation

If you use PREMISE in your research, please cite:

> Vijendran S. *PREMISE: Probabilistic Read-level Expectation Maximization for Integrated Source Estimation.* (manuscript in preparation)

## License

See [LICENSE](LICENSE).
