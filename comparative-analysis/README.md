# PREMISE Comparative Analysis

This directory holds the **code** for the comparative evaluation of PREMISE against [Centrifuger](https://github.com/mourisl/centrifuger), [Ganon](https://github.com/pirovc/ganon), [Karp](https://github.com/mreppell/Karp), [KMCP](https://github.com/shenwei356/kmcp), [MORA](https://github.com/AlgoLab/MORA), and [Sylph](https://github.com/bluenote-1577/sylph).

## Dependencies

All dependencies are managed via ```flake.nix```.

```bash
nix develop .#benchmark              # creates a temporary shell environment with all binaries needed for the comparative analysis on PATH
```

## Data

Every data path in the benchmark resolves under `[run] bench_data` in `params.toml`; leaving it empty defaults to the `premise/comparative-analysis/`. All scripts read it from the environment, so the driver exports it for them. Point it at a data disk to keep the checkout small:

```toml
[run]
bench_data = "/data/"
```

The reference database and the simulated samples are archived on Zenodo and are **downloaded by hand**; the real samples are on NCBI SRA and are to be fetched by accession.

**Step 1** Download the archive from Zenodo:

<!-- TODO at publication: add the Zenodo link, and the checksums for the archive alongside it. -->
> **Zenodo:** _<link to be added>_
>
> _Checksums to be published with the deposit._

and unpack it and point `bench_data` to it.

```bash
mkdir -p "$BENCH_DATA"
tar -xf premise-benchmark-data.tar.gz -C "$BENCH_DATA"
```

Verify the layout before going further — every later step assumes it:

```
$BENCH_DATA/indexes/sequences.fasta              # all reference sequences
$BENCH_DATA/indexes/{names.dmp,nodes.dmp,seqid2taxid.map,sequences.tax}   # taxonomy, complete
$BENCH_DATA/indexes/refs/
$BENCH_DATA/samples/synthetic/{isolate,mixed,mixed-subtype}/Dataset-<N>/
$BENCH_DATA/samples/real/{isolate,mixed}/<SRR>/true_sources.fasta      # the sample's 16 reference segments
```

**Step 2:** Preprocess the real samples

```bash
./scripts/fetch_real_samples.py     # fetch the 8 runs from NCBI SRA by accession
./scripts/prepare_real_samples.py   # trim + filter both splits, and derive their truth
./scripts/read_truth.py             # Assign truth sources to the preprocessed reads using true_sources.fasta as the reference.

```
The two real splits are not classified as they arrive from SRA. `scripts/prepare_real_samples.py` handles both, and emits the reads **and** the per-read truth in one pass, because the reads a method is scored on and the reads truth covers must be the same set. ```scripts/prepare_real_samples.py``` also includes additional options for split/sample specific steps; see ```scripts/prepare_real_samples.py --help``` for the full list of options.

Every read path is relative to the sample's own directory, `samples/real/{isolate,mixed}/<SRR>/`

| Output | Contents |
|---|---|
| `truth_assignments.tsv` | `read_id` \t `reference_accession` |
| `unclassified_reads.txt` | `read_id` of pairs `read_truth.py` could not resolve |

Both are written by `scripts/read_truth.py` in a single pass, piped straight from the aligner.

> **Note:** `run-analysis.py` runs this itself, for **both** real splits, whenever a sample's filtered reads or truth are missing. Running it by hand is only needed to force a rebuild (`--force`) or to prepare ahead of time. The `true_sources.fasta` for each file holds the references present in the sample; preprocessing scripts index`true_sources.fasta` for each sample separately to ensure a read is only ever assigned within the genome the sample actually contains.

## Running the analysis

```bash
./run-analysis.py --threads <n>       # full run: build indexes + all methods + all splits
./run-analysis.py --dry-run           # print every command it would run, run nothing
```

> **Warning:** by default this deletes and rebuilds every index under `indexes/<method>/` and overwrites `db_build.csv`.
> Pass `--skip-build` to reuse existing indexes.

Method codes: `pre` PREMISE, `cen` Centrifuger, `gan` Ganon, `kap` Karp, `kmp` KMCP, `mor` MORA, `syl` Sylph.

A failed index build aborts the run, and each build begins by deleting its own index directory.

## Configuration

Everything lives in `params.toml`: the `[run]` section holds the driver-specific paramters, one `[<code>]` section for method-specific parameters, and `[iss]` holds the simulation parameters. `scripts/ablation.py` reads the same file, so its sweep baseline cannot drift from the configuration published here.

| Section | Holds | Read by |
|---|---|---|
| `[run]` | data root, thread count, classification timeout, which splits/methods run, kmcp scope | `run-analysis.py` |
| `[pre]` … `[syl]` | per-method classification and index-build parameters (seed length, EM settings, k-mer sizes) | `run-analysis.py`, `scripts/ablation.py` |
| `[iss]` | InSilicoSeq cpu count and the per-mode RNG seeds | `scripts/make_synthetic_datasets.py` |


```bash
python3 run-analysis.py --methods "pre cen syl" --splits syn-iso --max-ds 1 --skip-build
RUN_METHODS="pre cen syl" RUN_SPLITS=syn-iso MAX_DS_PER_SPLIT=1 SKIP_BUILD=1 python3 run-analysis.py
```

## Splits

The benchmark is organised into five splits, all of which `RUN_SPLITS` runs by default. Each split maps to one directory under `samples/`, and results mirror that layout under `results/<method>/`.

| Split | `samples/` path | Contents |
|---|---|---|
| `syn-iso` | `synthetic/isolate` | 4 single-genome simulated samples |
| `syn-mix` | `synthetic/mixed` | 4 three-genome simulated mixtures |
| `syn-mix-subtype` | `synthetic/mixed-subtype` | 4 same-subtype simulated mixtures |
| `real-iso` | `real/isolate` | 4 real SRA runs |
| `real-mix` | `real/mixed` | 4 real mixed-infection runs |


> KMCP is skipped on the three synthetic splits by design: `kmcp search` runs at roughly 4.6k queries/min on the 500k-pair simulated samples (~3.6 h each), which exceeds the classification timeout.
> Override with `KMCP_ALLOW_SYNTHETIC=1`.

## Directory structure

All scripts and parameter files are git tracked, while data files are not. The directory tree after the full comparative analysis pipeline is written below.

```
comparative-analysis/
├── run-analysis.py         # load params.toml, then call scripts/ in order
├── params.toml             # all configuration
├── scripts/
│                           #   fetch_real_samples.py      real reads, from SRA
│                           #   prepare_real_samples.py    trim/filter both real splits + truth
│                           #   read_truth.py              SAM -> truth_assignments.tsv + unclassified_reads.txt
│                           #   make_synthetic_datasets.py syn-iso | syn-mix | syn-mix-subtype generators
│                           #   clean_sequences.py         drop duplicate/substring records
│                           #   analyze.py + bench_metrics.py   scoring and the LaTeX tables
│                           #   normalize_*.py             per-method output -> common schema
│                           #   load_params.py             params.toml -> Python dict / one value
│                           #   utils.py                   code root & $BENCH_DATA, resolved lazily
│                           #   clean_sequences.py         sequences.fasta -> sequences-cleaned
│                           #   ablation{,_report,_figures}.py  PREMISE parameter sweep harness
├── db_build.csv            # Per-method index build time + on-disk size (written by the driver)
├── ablation/figs/          # Ablation figures (see Parameter ablation, below)
│
├── indexes/                # Reference data and tool-specific indexes
│   ├── sequences.fasta          # Combined reference FASTA, as obtained
│   ├── sequences-cleaned.fasta  # Non-redundant subset — THIS is what every index is built from
│   ├── sequences-dropped.tsv    # Every dropped record, its reason, and what superseded it
│   ├── sequences.tax       # Per-sequence taxonomy (Karp)
│   ├── names.dmp           # NCBI taxonomy names table
│   ├── nodes.dmp           # NCBI taxonomy nodes table
│   ├── seqid2taxid.map     # Sequence ID → NCBI taxon ID (Centrifuger)
│   ├── refs/               # Per-accession FASTA files
│   └── <method>/           # Per-tool index, all generated by the driver
│
├── samples/                # Input reads
│   ├── real/                       # each <SRR>/ also holds, for both splits:
│   │   │                           #   true_sources.fasta       the sample's reference set
│   │   │                           #   truth_assignments.tsv    read_id -> segment (classified)
│   │   │                           #   unclassified_reads.txt   read_ids with no truth
│   │   │                           #   <SRR>.{cutadapt,bwa,truth}.log
│   │   ├── isolate/<SRR>/          # <SRR>_{1,2}.fastq raw, *_{1,2}.ca.fastq (primer-trimmed),
│   │   │                           #   *_{1,2}-filtered.ca.fastq (analysis-ready reads)
│   │   └── mixed/<SRR>/            # <SRR>_{1,2}.fastq raw, *_{1,2}.ca.fastq (adapter-trimmed),
│   │                               #   *-filtered.ca.fastq (analysis-ready reads)
│   └── synthetic/
│       └── {isolate,mixed,mixed-subtype}/
│           ├── seeds.tsv               # conatains random seeds for simulation reproducibility
│           └── Dataset-<N>/
│               ├── src-abundance.txt   # True per-source abundance (ref_id \t proportion)
│               └── all-abundances.txt  # True abundance for all references (incl. zeros)
│
└── results/<method>/{real,synthetic}/{isolate,mixed,mixed-subtype}/<sample>/   # Note: mixed-subtype will appear only in synthetic 
    ├── <sample>.*          # Method output (see below)
    └── <sample>.log        # Method stderr/stdout
```


## Samples

### Real (SRA)

| Split | Accessions |
|---|---|
| `real-iso` | SRR31013463, SRR31013465, SRR31013467, SRR31013473 |
| `real-mix` | SRR3360139, SRR3360140, SRR3360145, SRR3360146 |


### Synthetic

All simulated with [InSilicoSeq](https://github.com/HadrienG/InSilicoSeq) 2.0.1 using the `MiSeq` error model at 500,000 read pairs of 301 bp per sample. The simulated **`syn-iso`** and **`syn-mix`** samples are composed of reassortant-like chimeras: one sequence per influenza segment (1–8), each drawn from a *different* strain. The `syn-iso` samples have one genome per sample while `syn-mix` samples have three genomes each.

Both ship in the Zenodo archive; `scripts/make_synthetic_datasets.py` can be used to rebuild them. It has three modes and sources are always drawn from `indexes/sequences-cleaned.fasta`, so no simulated read can come from a sequence the indexes lack and truth is parsed from the ISS read names by `analyze.py` rather than written to a file:
```bash
./scripts/make_synthetic_datasets.py syn-iso samples/synthetic/isolate
./scripts/make_synthetic_datasets.py syn-mix samples/synthetic/mixed
./scripts/make_synthetic_datasets.py syn-mix-subtype samples/synthetic/mixed-subtype

# one dataset using the corresponding random seed
./scripts/make_synthetic_datasets.py syn-iso samples/synthetic/isolate --only 3
```

For `syn-mix-subtype` samples, instead of chimeras it uses three **real** same-subtype strains per dataset, each contributed as its own complete 8-segment set.

> Reproducibility: all parameters used in simulation live in `params.toml` under `[iss]` — `cpus`, and the parent seeds `seed_syn_iso` / `seed_syn_mix` / `seed_syn_mix_subtype`. Each dataset's RNG is derived from its parent seed and its index alone, so a whole split is a pure function of that seed, `--cpus`, and `indexes/sequences-cleaned.fasta`; nothing depends on what is already in the output tree The simulation seed each dataset used is recorded in `seeds.tsv` at the split root, so the archive can be audited without re-running anything. The `--cpus` and `--seed` flags override the defaults for one-off regenerations. ISS is deterministic only at a **fixed `--cpus`**: both the work partition and the per-worker seeding depend on it, and synthetic truth is parsed from ISS read names, so changing it changes the truth as well as the reads. It is not a performance knob.

### Datasets and truth

All four splits, four samples each:

| Split | `sub` | Samples | Reads fed |
|---|---|---|---|
| `syn-isolate` | synthetic/isolate | Dataset-1…4 | `reads_R{1,2}.fastq` |
| `syn-mixed` | synthetic/mixed | Dataset-1…4 | `reads_R{1,2}.fastq` |
| `syn-mixed-subtype` | synthetic/mixed-subtype | Dataset-1…4 | `reads_R{1,2}.fastq` |
| `real-isolate` | real/isolate | SRR31013463/65/67/73 | `<SRR>_{1,2}-filtered.ca.fastq` |
| `real-mix` | real/mixed | SRR3360139/40/45/46 | `<SRR>_{1,2}-filtered.ca.fastq` |

Truth comes from ISS read names on the synthetic splits (exact, no aligner in the loop) and from the sample's own `truth_assignments.tsv` on both real splits. There the per-read truth and the abundance truth are derived from that one file in a single pass, so they cannot disagree.

## Reference database

The reference set is 6,524 influenza records: **6,508 GenBank/RefSeq accessions**, listed one per line in `indexes/accessions.txt`, plus **16 local `PR8_*`/`WSN33_*` segments** that carry no accession under those names; the strains spiked into the real-mixed samples, on which `analyze.py::_strain_of` keys real-mixed precision.

### Cleaning

The sequences are cleaned by removing duplicate entries with different accessions, entries with long consecutive subsequences of ambiguous characters (N's), and entries that are contained in another longer entry.

```bash
./scripts/clean_sequences.py     # ~1 min
```

This reads `indexes/sequences.fasta` and writes `indexes/sequences-cleaned.fasta` (5,264 records) plus `indexes/sequences-dropped.tsv`, a manifest naming every dropped record, why it was dropped, and the record that superseded it. Indexes are then built from `indexes/sequences-cleaned.fasta`.

## Parameter ablation

A separate study from the method comparison above: a one-at-a-time sweep of PREMISE's five tunable parameters, measuring how each trades off runtime, peak memory, precision, coverage, Ruzicka distance and Jaccard distance. It shares this directory's `scripts/`, `params.toml` and `$BENCH_DATA`. To run the ablation study:

```bash
./scripts/ablation.py --dry-run                   # print the plan
./scripts/ablation.py                             # all 4 splits, 16 datasets
./scripts/ablation.py --split syn-mix            # one split
./scripts/ablation.py --split syn-mix --dataset Dataset-3
./scripts/ablation.py --only mem                  # one parameter
./scripts/ablation.py --force                     # ignore cached points
./scripts/ablation.py --keep-large                # keep .aligns/.posteriors/.matches
```

To aggregate and visualize the results of the ablation study:
```bash
./scripts/ablation_report.py                      # text + results/ablation-tables.generated.tex
./scripts/ablation_figures.py                     # ablation/figs/<split>/*.{pdf,png}
./scripts/ablation_figures.py --split syn-mix --only mem,rho --formats pdf,png,svg
./scripts/ablation_figures.py --no-relative
```
