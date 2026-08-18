# PREMISE:  A probabilistic framework for source assignment of viral Illumina reads

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

1. Builds a compressed FM-index over the reference sequences.
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

**Step 1:** Build a bidirectional FM-index

```bash
premise build -s <reference.fasta>
```

Produces `<reference>.fmidx`. This index is required for both the CLI and GUI query steps.

**Step 2:** Analyze the reads in sample

```bash
premise query \
  -s <reference.fmidx>     \  # FM-index built in Step 1
  -1 <R1.fastq.gz> \
  -2 <R2.fastq.gz> \
  -m <min_seed_length>     \  # minimum MEM seed length (default 22)
  --eps_1 <float>          \  # alignment likelihood cutoff (default 0 = disabled)
  --eps_2 <float>          \  # minimum match log-probability (default 1e-18)
  --rho   <float>          \  # EM penalty weight $\rho$ (default 150)
  --omega <float>          \  # EM penalty weight $\omega$ (default 1e-10)
  --iter  <int>            \  # EM iterations (default 100)
  --em_threshold <float>   \  # EM convergence threshold (default 1e-6)
  --no-penalty             \  # disable the L1 penalty (plain EM)
  -t <threads>             \  # default 2; 0 = all available cores
  -o <output_prefix>
```

Outputs:
| File | Contents |
|------|----------|
| `<output>.matches` | Per-read assignments (TSV) |
| `<output>.posteriors` | Per-read posterior probabilities (TSV) |
| `<output>.props` | Reference abundance proportions (TSV) |
| `<output>.aligns` | Raw per-alignment likelihoods, one row per (read, reference) (TSV) |

Run `premise query -h` for the full option list.

### Interactive GUI (optional)

```bash
premise server
```

Opens a browser UI at `http://localhost:8080` with drag-and-drop file upload, interactive results tables, pie chart, and EM convergence plot.

## Algorithm

PREMISE seeds each read with **Super-Maximal Exact Matches (SMEMs)** found via the reference FM-index. Each seed is projected onto a reference diagonal, and the full read is then scored **ungapped** against that offset --- there is no chaining and no gapped extension. Read-level alignment log-likelihoods are computed from base quality scores (Phred-scaled error probabilities in natural log space); `-m`/`--mem` sets the minimum seed length.

The EM step solves a penalized likelihood maximization:

- $\rho$ controls an L1-style sparsity penalty on the proportion vector.
- $\omega$ is a small regularization floor.
- Convergence is tracked by the total data log-likelihood across iterations.

Parameters $\eps_1$ and **$\eps_2$** control alignment filtering: $\eps_1$ is a minimum alignment likelihood threshold (linear space); $\eps_2$ is a minimum match log-probability per read. The default $\eps_1 = 0$ disables the cutoff entirely.

## Project Structure

```
premise/
├── src/
│   ├── main.rs          # CLI, HTTP server, EM algorithm, alignment logic
│   ├── utils.rs         # Quality-score utilities, match log-probability
│   └── templates/
│       ├── index.html   # Browser GUI markup (embedded at compile time)
│       ├── styles.css   # Pico.css overrides (embedded at compile time)
│       └── app.js       # Frontend logic — D3 charts, dropzones, dark mode
├── tests/               # Integration tests (CLI + local server)
├── Cargo.toml
└── README.md
```

## Output Format
PREMISE produced four files: The final per-read assignments in ```*.matches```, the read-reference alignment probabilities for all read-reference matches founf prior to EM in ```*.aligns```, the posterior probabilities of all read matches found after EM in ```*.posteriors```, and the proportions of detected references after EM in ```*.props```. Below are the schemas of each output.

### `.matches` (TSV)
| Column | Description |
|--------|-------------|
| `read_id` | Read identifier |
| `ref_id` | Assigned reference sequence ID |
| `posterior` | Posterior probability of assignment |


### `.aligns` (TSV)
| Column | Description |
|--------|-------------|
| `read_id` | Read identifier |
| `ref_id` | Assigned reference sequence ID |
| `probability` | Probability of read-reference alignment|

### `.posteriors` (TSV)
| Column | Description |
|--------|-------------|
| `read_id` | Read identifier |
| `ref_id` | Assigned reference sequence ID |
| `posterior` | Posterior probability of assignment |


### `.props` (TSV)
| Column | Description |
|--------|-------------|
| `ref_id` | Reference sequence ID |
| `proportion` | Estimated relative abundance |

The reported proportions are those inferred by EM, not normalizing the read counts in `.matches`. After reads whose reference was pruned by EM are reclassified as unclassified, the likelihood matrix is restricted to the remaining reads and one last E-step and M-step are run; `.props` reports that result.

## Comparative Analysis

PREMISE is evaluated against six other read-classification tools --- [Centrifuger](https://github.com/mourisl/centrifuger), [Ganon](https://github.com/pirovc/ganon), [Karp](https://github.com/mreppell/Karp), [KMCP](https://github.com/shenwei356/kmcp), [MORA](https://github.com/AlgoLab/MORA), and [Sylph](https://github.com/bluenote-1577/sylph) --- across five dataset splits: three of which are simulated (synthetic isolate, synthetic mixed, synthetic same-subtype mixed) and two are real samples acquired from NCBI SRA (real isolate and real mixed-infection samples). Every method is evaluated on runtime, peak memory, per-read precision, coverage, Ruzicka distance, and Jaccard distance. A separate abaltion study that sweeps over the parameters `-m`, `--eps_1`, `--eps_2`, `--rho` and `--omega` is also implemented.

The toolchain for every competing method is pinned in `flake.nix`:

```bash
nix develop .#benchmark              # creates a temporary shell environment with all binaries needed for the comparative analysis on PATH
cd comparative-analysis
python3 run-analysis.py --threads <n>
```

See [comparative-analysis/README.md](comparative-analysis/README.md) for data preparation, configuration (`params.toml`), directory layout, and the full ablation procedure.

## Citation

If you use PREMISE in your research, please cite:

> Vijendran S. *PREMISE: Probabilistic Read-level Expectation Maximization for Integrated Source Estimation.* (manuscript in preparation)

## License

See [LICENSE](LICENSE).
