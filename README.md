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

## Usage

### Step 1 — Build an FM-index

```bash
premise build -s <reference.fasta>
```

Produces `<reference>.fmidx`. This index is required for both the CLI and GUI query steps.

### Step 2 — Classify reads (CLI)

```bash
premise query \
  -s <reference.fasta> \
  -1 <R1.fastq.gz> \
  -2 <R2.fastq.gz> \
  -p <percent_mismatch>   \  # e.g. 5
  --eps_1 <float>          \  # alignment likelihood cutoff (default 1e-4)
  --eps_2 <float>          \  # minimum match log-probability (default 1e-18)
  --rho   <float>          \  # EM penalty weight ρ (default 20)
  --omega <float>          \  # EM penalty weight ω (default 1e-20)
  --iter  <int>            \  # EM iterations (default 100)
  -t <threads>             \  # 0 = all available cores
  -o <output_prefix>
```

Outputs:
| File | Contents |
|------|----------|
| `<output>.matches` | Per-read assignments (TSV) |
| `<output>.posteriors` | Per-read posterior probabilities (TSV) |
| `<output>.props` | Reference abundance proportions (TSV) |

Run `premise query -h` for the full option list.

### Step 3 — Interactive GUI (optional)

```bash
premise server
```

Opens a browser UI at `http://localhost:8080` with drag-and-drop file upload, interactive results tables, pie chart, and EM convergence plot. Supports light/dark mode.

## Algorithm

PREMISE uses a seeded alignment strategy based on exact FM-index lookups, followed by seed extension with a configurable mismatch tolerance. Read-level alignment log-likelihoods are computed using base quality scores (Phred-scaled error probabilities in natural log space).

The EM step solves a penalized likelihood maximization:

- **ρ** controls an L1-style sparsity penalty on the proportion vector.
- **ω** is a small regularization floor.
- Convergence is tracked by the total data log-likelihood across iterations.

Parameters **ε₁** and **ε₂** control alignment filtering: ε₁ is a minimum alignment likelihood threshold (linear space); ε₂ is a minimum match log-probability per read.

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
├── pkg/                 # WASM build artefacts (experimental)
├── eval_tool/           # Evaluation scripts and notebooks
├── Cargo.toml
└── README.md
```

## Output Format

### `.matches` (TSV)
| Column | Description |
|--------|-------------|
| `read_id` | Read identifier |
| `ref_id` | Assigned reference sequence ID |
| `posterior` | Posterior probability of assignment |

### `.posteriors` (TSV)
Full posterior probability matrix: one row per read, one column per reference.

### `.props` (TSV)
| Column | Description |
|--------|-------------|
| `ref_id` | Reference sequence ID |
| `proportion` | Estimated relative abundance |

## Citation

If you use PREMISE in your research, please cite:

> Vijendran S. *PREMISE: Probabilistic Read-level Expectation Maximization for Integrated Source Estimation.* (manuscript in preparation)

## License

See [LICENSE](LICENSE).
