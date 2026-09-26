# PREMISE Comparative Analysis

This directory holds the **code** for the comparative evaluation of PREMISE against [Centrifuger](https://github.com/mourisl/centrifuger), [Ganon](https://github.com/pirovc/ganon), [Karp](https://github.com/mreppell/Karp), [KMCP](https://github.com/shenwei356/kmcp), [MORA](https://github.com/AlgoLab/MORA), [Salmon](https://github.com/COMBINE-lab/salmon), and [Sylph](https://github.com/bluenote-1577/sylph), and for the parameter ablation of PREMISE itself.

```bash
nix develop ..#benchmark                 # every tool and Python dependency, pinned by flake.nix
python3 -m premise_bench --help          # the commands
python3 -m premise_bench <command> --help
```

| Command | Does |
|---|---|
| `run` | build every index, classify every split, score them (the end-to-end driver) |
| `analyze` | score existing results into `results/tables/comparative-<split>.csv` |
| `ablation` | one-at-a-time sweep of PREMISE's parameters |
| `fetch-real` | download the real samples from NCBI SRA |
| `prepare-real` | trim, filter and derive truth for the real samples; add the decoy references |
| `read-truth` | per-read truth from a bwa-mem2 stream (the pipe stage `prepare-real` uses) |
| `synth` | regenerate the synthetic datasets with InSilicoSeq |
| `clean-db` | drop duplicate and contained reference records |
| `normalize` | reshape one tool's native output by hand |
| `params` | print `params.toml` as shell assignments |

`./run-analysis.py` is kept as a shortcut for `python3 -m premise_bench run`.

Method codes: `pre` PREMISE, `cen` Centrifuger, `gan` Ganon, `kap` Karp, `kmp` KMCP, `mor` MORA, `sal` Salmon, `syl` Sylph.

> **Note:** `run` executes the full comparative analysis and ablation study.

## Data

Every data path resolves under `[run] bench_data` in `params.toml` (or `$BENCH_DATA`, or `run --bench-data`); leaving it empty defaults to this directory.

```toml
[run]
bench_data = "/data/"
```

The reference database and the simulated samples are archived on Zenodo. The real samples are on NCBI SRA and are fetched by accession.

**Step 1:** Download the archive from Zenodo:

> **Zenodo:** [10.5281/zenodo.19026876](https://doi.org/10.5281/zenodo.19026876)

then unpack it and point `bench_data` to it.

```bash
mkdir -p "$BENCH_DATA"
tar -xf premise-benchmark-data.tar.gz -C "$BENCH_DATA"
```

**Step 2:** Preprocess the real samples

```bash
python3 -m premise_bench fetch-real      # fetch the 8 runs from NCBI SRA by accession
python3 -m premise_bench prepare-real    # trim + filter both splits, derive their truth, add the decoys
```

The two real splits are not classified as they arrive from SRA. `prepare-real` handles both, and emits the reads **and** the per-read truth in one pass. It aligns each sample against its own `true_sources.fasta`.

Every path is relative to the sample's own directory, `samples/real/{isolate,mixed}/<SRR>/`:

| Output | Contents |
|---|---|
| `truth_assignments.tsv` | `read_id` \t `reference_accession` |
| `unclassified_reads.txt` | `read_id` of pairs `read-truth` could not resolve |

`prepare-real` then adds the decoy references to `indexes/` (see [Decoys](#decoys-near-twins-of-the-real-sources)).


## Running the analysis

```bash
python3 -m premise_bench run --threads <n>       # full run: build indexes + all methods + all splits
python3 -m premise_bench run --dry-run           # print every command it would run, run nothing
python3 -m premise_bench run --build-only        # stop after the indexes
```

> **Warning:** by default this deletes and rebuilds every index under `indexes/<method>/` and overwrites `db_build.csv`.
> Pass `--skip-build` to reuse existing indexes.

A failed index build aborts the run (unless `--keep-going`), and each build begins by deleting its own index directory. Every classification runs under the `classify_timeout` wall-clock cap, with its time and peak memory recorded in the dataset's `time-mem` file.

The benchmark writes results for every split into `results/tables/comparative-<split>.csv`: one tidy row per (method, dataset, metric) with a `status` column marking timed-out and missing values.

## Splits

The benchmark is organised into five splits, all of which run by default. Each split maps to one directory under `samples/`, and results mirror that layout under `results/<method>/`.

| Split | `samples/` path | Contents | Reads fed |
|---|---|---|---|
| `syn-iso` | `synthetic/isolate` | 4 single-genome simulated samples (Dataset-1…4) | `reads_R{1,2}.fastq` |
| `syn-mix-sub` | `synthetic/mixed` | 4 three-genome simulated mixtures | `reads_R{1,2}.fastq` |
| `syn-mix-strain` | `synthetic/mixed-subtype` | 4 same-subtype simulated mixtures | `reads_R{1,2}.fastq` |
| `real-iso` | `real/isolate` | SRR31013463, SRR31013465, SRR31013467, SRR31013473 | `<SRR>_{1,2}-filtered.ca.fastq` |
| `real-mix` | `real/mixed` | SRR3360139, SRR3360140, SRR3360145, SRR3360146 | `<SRR>_{1,2}-filtered.ca.fastq` |

Read truth assignments comes from the InSilicoSeq read names on the synthetic splits, and from the sample's own `truth_assignments.tsv` on both real splits.

> KMCP is skipped on the three synthetic splits by design: `kmcp search` runs at roughly 4.6k queries/min on the 500k-pair simulated samples (~3.6 h each), which exceeds the classification timeout.
> Override with `KMCP_ALLOW_SYNTHETIC=1`.

## Code layout

```
comparative-analysis/
├── params.toml                 # all configuration
├── run-analysis.py             # shortcut for `python3 -m premise_bench run`
└── premise_bench/
    ├── __main__.py             # the command table
    ├── config.py               # params.toml, the code root and $BENCH_DATA paths
    ├── splits.py               # the five splits: folders, datasets, read-file names
    ├── utils.py                # shared helpers (gzip reading, label and read-id conventions, time-mem parsing, FASTA/k-mers)
    ├── metrics.py              # Ruzicka/Jaccard distances, reference-level FP/FN, per-read precision
    ├── runner.py               # child processes, the wall-clock cap and resource accounting
    ├── pipeline.py             # `run`: build -> classify -> score
    ├── methods/                # one module per tool: build, classify, normalize (if needed), load
    │   └── __init__.py         #   the registry: codes, binaries, index files, scoring flags
    ├── evaluate/               # truth.py, score.py, analyze.py
    ├── data/                   # fetch_real, prepare_real (+ decoys), read_truth, synthetic, clean_db
    └── ablation/sweep.py       # `ablation`
```

## Data layout

The data archive holds only the **inputs** both studies need.

```
$BENCH_DATA/
├── db_build.csv                     ← run: per-method index build time + on-disk size
├── indexes/
│   ├── sequences.fasta              # combined reference FASTA, as obtained (input of clean-db)
│   ├── sequences-cleaned.fasta      # non-redundant subset — every index is built from it
│   │                                #   ← prepare-real / run append the DECOY_* records
│   ├── seqid2taxid.map              # sequence ID → NCBI taxon ID (← decoy rows appended too)
│   ├── names.dmp, nodes.dmp         # NCBI taxonomy
│   ├── decoys.{fasta,tsv}, decoy-targets.tsv   ← prepare-real / run: the decoys and how they were drawn
│   ├── sequences.tax, sequences-cleaned.fasta.fai   ← run (Karp build)
│   └── <method>/                    ← run: per-tool index
├── samples/
│   ├── real/{isolate,mixed}/<SRR>/
│   │   ├── <SRR>_{1,2}.fastq                 # raw reads (from SRA; input of prepare-real)
│   │   ├── <SRR>_{1,2}-filtered.ca.fastq     # analysis-ready reads (what is classified)
│   │   ├── truth_assignments.tsv             # read_id -> true segment
│   │   └── true_sources.fasta                # the sample's reference set (truth + decoys are derived from it)
│   │                                         #   prepare-real --force also writes *.ca.fastq, unclassified_reads.txt, logs
│   └── synthetic/{isolate,mixed,mixed-subtype}/Dataset-<N>/
│       ├── reads_R{1,2}.fastq                # simulated reads; truth is in the read names
│       └── genomes.tsv                       # mixed-subtype only: its strains (the decoys' similarity targets)
└── results/                         ← run and ablation
    ├── params.used.toml             # the params.toml of the last run
    ├── tables/comparative-<split>.csv
    ├── ablation/<split>/ablation.csv (+ ablation/indexes/, the per-rate indexes)
    └── <method>/{real,synthetic}/{isolate,mixed,mixed-subtype}/<sample>/
        ├── <sample>.*               # method output
        ├── <sample>.log             # method stderr/stdout
        └── time-mem                 # wall time, peak memory, exit status
```

## Synthetic samples

All simulated with [InSilicoSeq](https://github.com/HadrienG/InSilicoSeq) 2.0.1 using the `MiSeq` error model at 500,000 read pairs of 301 bp per sample. The **`syn-iso`** and **`syn-mix-sub`** samples are reassortant-like chimeras: one sequence per influenza segment (1–8), each drawn from a *different* strain. `syn-iso` has one genome per sample, `syn-mix-sub` three. **`syn-mix-strain`** instead uses three **real** same-subtype strains per dataset, each contributed as its own complete 8-segment set.

```bash
python3 -m premise_bench synth syn-iso samples/synthetic/isolate
python3 -m premise_bench synth syn-mix-sub samples/synthetic/mixed
python3 -m premise_bench synth syn-mix-strain samples/synthetic/mixed-subtype

# one dataset using the corresponding random seed
python3 -m premise_bench synth syn-iso samples/synthetic/isolate --only 3
```

## Reference database

The reference set is 6,524 influenza records: **6,508 GenBank/RefSeq accessions**, listed one per line in `indexes/accessions.txt`, plus **16 local `PR8_*`/`WSN33_*` segments** that carry no accession under those names.

### Cleaning

The sequences are cleaned by removing duplicate entries with different accessions, entries with long consecutive runs of ambiguous characters (N's).

```bash
python3 -m premise_bench clean-db     # ~1 min
```

This reads `indexes/sequences.fasta` and writes `indexes/sequences-cleaned.fasta` (5,268 records) plus `indexes/sequences-dropped.tsv`. Indexes are then built from `indexes/sequences-cleaned.fasta`.

### Decoys: near-twins of the real sources

In the real sample, every true segment's nearest database record is far away (31-mer containment ~0.2–0.4). `prepare-real` therefore derives, from each real true-source segment, `per_source` decoy sequences `DECOY_<accession>_<n>` at the similarity the three co-occurring strains of a `syn-mix-strain` sample have to each other, and appends them to `indexes/sequences-cleaned.fasta` and `indexes/seqid2taxid.map`.

The parameters live in `[decoy]` of `params.toml` (`seed`, `per_source`, `ti_frac`, `k`). Every draw comes from `Random(f"{seed}:{accession}")`. `prepare-real` regenerate them on every call and rewrite the reference set only when they differ from the ones it holds.

```bash
python3 -m premise_bench prepare-real --decoys-only                  # just the decoys
python3 -m premise_bench prepare-real --decoys-only --verify-decoys  # + each real source's nearest neighbours
python3 -m premise_bench prepare-real --no-decoys                    # reads and truth only
```

Alongside the reference set, `indexes/` gets `decoys.fasta`, `decoys.tsv` and `decoy-targets.tsv`.

## Parameter ablation

A one-at-a-time sweep of PREMISE's tunable parameters, measuring how each trades off runtime, peak memory, precision, coverage, Ruzicka distance and Jaccard distance.

```bash
python3 -m premise_bench ablation --dry-run                   # print the plan
python3 -m premise_bench ablation                             # every split, 20 datasets
python3 -m premise_bench ablation --split syn-mix-sub         # one split
python3 -m premise_bench ablation --split syn-mix-sub --dataset Dataset-3
python3 -m premise_bench ablation --only mem                  # one parameter
python3 -m premise_bench ablation --only sa_sample_rate       # just the SA sampling-rate sweep
python3 -m premise_bench ablation --force                     # ignore cached points
python3 -m premise_bench ablation --keep-large                # keep .aligns/.posteriors/.matches
```

Results land in `results/ablation/<split>/ablation.csv`, one row per (dataset, parameter, value); the ablation figures are drawn from these CSVs.
