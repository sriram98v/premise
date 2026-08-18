#!/usr/bin/env python3
"""Regenerate the synthetic datasets from scratch with sources drawn from indexes/sequences-cleaned.fasta

    syn-iso          1 chimeric pseudo-genome, 8 segments, uniform 1/8 abundance
    syn-mix          3 chimeric pseudo-genomes, 24 sequences, Dirichlet(8,8,8) genome shares
    syn-mix-subtype  3 REAL same-subtype strains, mutually similar, Dirichlet(8,8,8) shares

Outputs, matching the existing dataset layout exactly:
    src.fasta            the chosen sequences, original headers, database order
    src-abundance.txt    ref_id \t proportion
    all-abundances.txt   every reference in the database
    reads_R{1,2}.fastq   written by InSilicoSeq, called in process (iss.app.generate_reads)
    genomes.tsv          syn-mix-subtype mode only: strain -> segment -> accession
    seeds.tsv            at the split root: the simulation seed each dataset used

Usage:
    python3 scripts/make_synthetic_datasets.py syn-iso samples/synthetic/isolate
    python3 scripts/make_synthetic_datasets.py syn-mix samples/synthetic/mixed
    python3 scripts/make_synthetic_datasets.py syn-mix-subtype samples/synthetic/mixed-subtype
    python3 scripts/make_synthetic_datasets.py syn-iso samples/synthetic/isolate --only 3
"""
from __future__ import annotations

import argparse
import logging
import random
import re
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import load_params
import utils
from iss.app import generate_reads


def db() -> Path:
    return utils.indexes() / "sequences-cleaned.fasta"

logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(name)s:%(message)s")

SEGMENTS = list(range(1, 9))
N_GENOMES = 3
N_DATASETS = 4
K = 21

SEG_RE = re.compile(r"segment (\d)")
STRAIN_RE = re.compile(r"Influenza A virus \((.+?)\) segment")
SUB_RE = re.compile(r"\(H(\d+)N(\d+)\)")

ALL_ABUND_NAME = {"syn-iso": "other-props.txt", "syn-mix": "all-abundances.txt",
                  "syn-mix-subtype": "all-abundances.txt"}

SEEDS_NAME = "seeds.tsv"

def _iss_int(key: str) -> int:
    """One [iss] value as an int. Every value in params.toml is a quoted string by design.
    """
    params = load_params.section("iss")
    try:
        raw = params[key]
    except KeyError:
        raise SystemExit(f"params.toml: [iss] has no '{key}' "
                         f"(have: {', '.join(sorted(params)) or 'none'})") from None
    try:
        return int(raw)
    except ValueError:
        raise SystemExit(f"params.toml: [iss].{key} must be an integer, got {raw!r}") from None


def read_headers(fasta: Path) -> dict[str, tuple[int | None, str]]:
    """Returns {accession: (segment_or_None, header_line)} for every record, in file order."""
    out: dict[str, tuple[int | None, str]] = {}
    with open(fasta) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            head = line[1:].rstrip("\n")
            acc = head.split()[0]
            m = SEG_RE.search(head)
            out[acc] = (int(m.group(1)) if m else None, head)
    return out


def dataset_rng(seed: int, index: int) -> random.Random:
    """The RNG for Dataset-<index>, independent of every other dataset in the split.
    """
    return random.Random(f"{seed}:{index}")


def indices(only: int | None) -> list[int]:
    """Which dataset numbers this run writes."""
    if only is None:
        return list(range(1, N_DATASETS + 1))
    if not 1 <= only <= N_DATASETS:
        sys.exit(f"ERROR: --only {only} is outside the {N_DATASETS} datasets in a split")
    return [only]


def write_seeds(out_root: Path, parent: int, rows: dict[int, int]) -> None:
    """Merge {index: iss_seed} into <out_root>/seeds.tsv, so the archive records its own provenance.
    """
    path = out_root / SEEDS_NAME
    merged: dict[int, tuple[int, int]] = {}
    if path.exists():
        with open(path) as f:
            next(f, None)
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) == 3 and parts[0].startswith("Dataset-"):
                    merged[int(parts[0].split("-", 1)[1])] = (int(parts[1]), int(parts[2]))
    for idx, iss_seed in rows.items():
        merged[idx] = (parent, iss_seed)
    with open(path, "w") as f:
        f.write("dataset\tparent_seed\tiss_seed\n")
        for idx in sorted(merged):
            p, s = merged[idx]
            f.write(f"Dataset-{idx}\t{p}\t{s}\n")
    print(f"seeds -> {path}")


def extract(fasta: Path, wanted: set[str], out: Path) -> None:
    """Copy the requested records out of a FASTA, preserving original headers and file order."""
    keep = False
    written = set()
    with open(fasta) as f, open(out, "w") as w:
        for line in f:
            if line.startswith(">"):
                acc = line[1:].split()[0]
                keep = acc in wanted
                if keep:
                    written.add(acc)
            if keep:
                w.write(line)
    missing = wanted - written
    if missing:
        sys.exit(f"ERROR: {len(missing)} sequences not found in {fasta}: {sorted(missing)[:5]}")


def write_abundances(out_dir: Path, mode: str, ordered_accs, per_seq: dict[str, float],
                     rows: list[str]) -> None:
    """Write src-abundance.txt (chosen only) and the all-references table.
    """
    with open(out_dir / "src-abundance.txt", "w") as f:
        f.writelines(f"{acc}\t{per_seq[acc]}\n" for acc in rows)
    with open(out_dir / ALL_ABUND_NAME[mode], "w") as f:
        f.writelines(f"{acc} {per_seq.get(acc, 0)}\n" for acc in ordered_accs)


def dirichlet_shares(rng: random.Random, n: int) -> list[float]:
    """n genome shares ~ Dirichlet(8,…,8) — the spread the existing mixtures sit in (0.22–0.44)."""
    gammas = [rng.gammavariate(8.0, 1.0) for _ in range(n)]
    total = sum(gammas)
    return [g / total for g in gammas]


def iss_args(out_dir: Path, reads: str, cpus: int, seed: int) -> argparse.Namespace:
    """The Namespace `iss generate` would have built for these parameters.

    generate_reads() reads its arguments off an argparse Namespace, so every option has to be
    present even when it is left at its default — the values below are iss's own parser defaults
    (iss/app.py, parser_gen). `abundance` stays at "lognormal" exactly as on the command line: it
    is ignored whenever abundance_file is set.
    """
    return argparse.Namespace(
        mode="kde", model="MiSeq", seed=seed, cpus=cpus,
        genomes=[str(out_dir / "src.fasta")],
        abundance_file=str(out_dir / "src-abundance.txt"),
        n_reads=reads, output=str(out_dir / "reads"),
        draft=None, ncbi=None, n_genomes_ncbi=None, n_genomes=None,
        abundance="lognormal", coverage=None, coverage_file=None,
        readcount_file=None, sequence_type="metagenomics",
        gc_bias=False, compress=False, store_mutations=False,
        fragment_length=None, fragment_length_sd=None,
    )


def run_iss(out_dir: Path, reads: str, cpus: int, seed: int) -> None:
    """Simulate reads into out_dir, then drop the intermediates. ISS 2.x writes reads_R{1,2}.fastq."""
    args = iss_args(out_dir, reads, cpus, seed)
    print(f"\niss generate  genomes={args.genomes[0]}  abundance={args.abundance_file}\n"
          f"              n_reads={reads}  cpus={cpus}  seed={seed}  model=MiSeq\n"
          f"              output={args.output}", flush=True)
    try:
        generate_reads(args)
    except SystemExit as e:
        raise SystemExit(f"iss generate failed for {out_dir}: {e}") from e
    for tmp in out_dir.glob("reads.iss.tmp*"):
        tmp.unlink()


def segment_pools(headers, need: int) -> dict[int, list[str]]:
    """Returns {segment: sorted candidate accessions}, sorted so the draw is seed-reproducible."""
    by_seg: dict[int, list[str]] = {s: [] for s in SEGMENTS}
    for acc, (seg, _) in headers.items():
        if seg in by_seg:
            by_seg[seg].append(acc)
    for s in SEGMENTS:
        by_seg[s].sort()
        if len(by_seg[s]) < need:
            sys.exit(f"ERROR: only {len(by_seg[s])} candidates for segment {s}, need {need}")
    return by_seg

def write_isolate_dataset(out_dir: Path, rng: random.Random, headers, by_seg,
                          reads: str, cpus: int) -> int:
    """One syn-iso dataset. Returns the simulation seed it used."""
    out_dir.mkdir(parents=True, exist_ok=True)

    genome = [rng.choice(by_seg[s]) for s in SEGMENTS]
    if len(set(genome)) != len(SEGMENTS):
        sys.exit("ERROR: drew the same accession for two segments")

    share = 1.0 / len(SEGMENTS)
    per_seq = {acc: share for acc in genome}

    extract(db(), set(per_seq), out_dir / "src.fasta")
    write_abundances(out_dir, "syn-iso", headers.keys(), per_seq, genome)

    # Drawn after the sources, so the simulation seed is a function of (parent seed, index) too.
    iss_seed = rng.randint(1, 2**31 - 1)

    print(f"{out_dir.name}: 8 segments @ {share}")
    for s, acc in zip(SEGMENTS, genome):
        print(f"  segment {s}: {acc}  {headers[acc][1][:90]}")

    run_iss(out_dir, reads, cpus, iss_seed)
    return iss_seed


def write_mixed_dataset(out_dir: Path, rng: random.Random, headers, by_seg,
                        reads: str, cpus: int) -> int:
    """One syn-mix dataset. Returns the simulation seed it used."""
    out_dir.mkdir(parents=True, exist_ok=True)

    chosen: list[list[str]] = []
    taken: set[str] = set()
    for _ in range(N_GENOMES):
        genome = []
        for s in SEGMENTS:
            pick = rng.choice([a for a in by_seg[s] if a not in taken])
            taken.add(pick)
            genome.append(pick)
        chosen.append(genome)

    shares = dirichlet_shares(rng, N_GENOMES)
    per_seq: dict[str, float] = {}
    for genome, share in zip(chosen, shares):
        for acc in genome:
            per_seq[acc] = share / len(SEGMENTS)

    extract(db(), set(per_seq), out_dir / "src.fasta")
    write_abundances(out_dir, "syn-mix", headers.keys(), per_seq,
                     [acc for genome in chosen for acc in genome])

    iss_seed = rng.randint(1, 2**31 - 1)

    print(f"{out_dir.name}: genome shares " + ", ".join(f"{s:.4f}" for s in shares))
    for i, genome in enumerate(chosen):
        print(f"  genome {i + 1} ({shares[i]:.4f}): {' '.join(genome)}")

    run_iss(out_dir, reads, cpus, iss_seed)
    return iss_seed


def isolate(out_root: Path, seed: int, reads: str, cpus: int,
            only: int | None = None) -> None:
    """syn-iso: one chimeric genome of eight segments at uniform abundance, per dataset."""
    out_root = Path(out_root)
    headers = read_headers(db())
    by_seg = segment_pools(headers, 1)

    seeds: dict[int, int] = {}
    for i in indices(only):
        seeds[i] = write_isolate_dataset(
            out_root / f"Dataset-{i}", dataset_rng(seed, i), headers, by_seg, reads, cpus)

    write_seeds(out_root, seed, seeds)
    print("done:", out_root)


def mixed(out_root: Path, seed: int, reads: str, cpus: int,
          only: int | None = None) -> None:
    """syn-mix: three chimeric genomes with Dirichlet-drawn shares, per dataset."""
    out_root = Path(out_root)
    headers = read_headers(db())
    by_seg = segment_pools(headers, N_GENOMES)

    seeds: dict[int, int] = {}
    for i in indices(only):
        seeds[i] = write_mixed_dataset(
            out_root / f"Dataset-{i}", dataset_rng(seed, i), headers, by_seg, reads, cpus)

    write_seeds(out_root, seed, seeds)
    print("done:", out_root)


def load_db():
    """Returns (strain -> {segment: accession}, accession -> sequence, accessions in file order)."""
    strains: dict[str, dict[int, str]] = defaultdict(dict)
    seqs: dict[str, list[str]] = {}
    order: list[str] = []
    cur = None
    with open(db()) as f:
        for line in f:
            if line.startswith(">"):
                head = line[1:].rstrip("\n")
                cur = head.split()[0]
                seqs[cur] = []
                order.append(cur)
                sm, gm = STRAIN_RE.search(head), SEG_RE.search(head)
                if sm and gm:
                    # setdefault: keep the first record for a (strain, segment) pair
                    strains[sm.group(1).replace(" (", "(")].setdefault(int(gm.group(1)), cur)
            elif cur:
                seqs[cur].append(line.strip())
    return strains, {a: "".join(v).upper() for a, v in seqs.items()}, order


def kmerize(seq: str) -> set[str]:
    return {seq[i:i + K] for i in range(len(seq) - K + 1)}


def candidate_groups(strains, seqs, min_jaccard: float):
    """Same-subtype cliques of N_GENOMES complete strains, all pairwise >= min_jaccard.

    Returns ([(subtype, min_pairwise, [strain, ...])] sorted by min_pairwise desc, complete_strains)
    """
    complete = {s: v for s, v in strains.items() if len(v) == len(SEGMENTS)}
    by_sub: dict[str, list[str]] = defaultdict(list)
    for s in complete:
        m = SUB_RE.search(s)
        if m:
            by_sub[f"H{m.group(1)}N{m.group(2)}"].append(s)

    gk = {s: kmerize("".join(seqs[segmap[i]] for i in SEGMENTS)) for s, segmap in complete.items()}

    def jac(a, b):
        ka, kb = gk[a], gk[b]
        return len(ka & kb) / len(ka | kb)

    out = []
    for sub, members in by_sub.items():
        if len(members) < N_GENOMES:
            continue
        members = sorted(members)
        best = None
        for combo in combinations(members, N_GENOMES):
            lo = min(jac(a, b) for a, b in combinations(combo, 2))
            if lo >= min_jaccard and (best is None or lo > best[0]):
                best = (lo, list(combo))
        if best:
            out.append((sub, best[0], best[1]))
    out.sort(key=lambda t: -t[1])
    return out, complete


def write_subtype_dataset(out_dir: Path, group, complete, order, rng,
                          reads: str, cpus: int, dry_run: bool) -> int:
    sub, lo, members = group
    out_dir.mkdir(parents=True, exist_ok=True)

    shares = dirichlet_shares(rng, N_GENOMES)
    per_seq: dict[str, float] = {}
    rows: list[str] = []
    table: list[tuple[str, int, str, float]] = []
    for strain, share in zip(members, shares):
        for s in SEGMENTS:
            acc = complete[strain][s]
            per_seq[acc] = share / len(SEGMENTS)
            rows.append(acc)
            table.append((strain, s, acc, share))

    if len(set(rows)) != N_GENOMES * len(SEGMENTS):
        sys.exit(f"ERROR: {out_dir.name}: {len(set(rows))} distinct accessions, "
                 f"expected {N_GENOMES * len(SEGMENTS)} (shared segment between strains?)")

    extract(db(), set(per_seq), out_dir / "src.fasta")
    write_abundances(out_dir, "syn-mix-subtype", order, per_seq, rows)

    with open(out_dir / "genomes.tsv", "w") as f:
        f.write("subtype\tmin_pairwise_jaccard\tstrain\tgenome_share\tsegment\taccession\n")
        for strain, s, acc, share in table:
            f.write(f"{sub}\t{lo:.4f}\t{strain}\t{share:.6f}\t{s}\t{acc}\n")

    iss_seed = rng.randint(1, 2**31 - 1)

    print(f"{out_dir.name}: {sub}  min pairwise Jaccard {lo:.3f}")
    for strain, share in zip(members, shares):
        print(f"    {share:.4f}  {strain}")

    if not dry_run:
        run_iss(out_dir, reads, cpus, iss_seed)
    return iss_seed


def subtype(out_root: Path, seed: int, reads: str, cpus: int,
            min_jaccard: float = 0.7, dry: bool = False, only: int | None = None) -> None:
    """Three real same-subtype strains per dataset, similarity-filtered."""
    strains, seqs, order = load_db()
    groups, complete = candidate_groups(strains, seqs, min_jaccard)

    print(f"{len(groups)} subtypes supply a mutually-similar group of {N_GENOMES} "
          f"at Jaccard >= {min_jaccard}:")
    for sub, lo, _ in groups:
        print(f"  {sub:6} min pairwise {lo:.3f}")
    if len(groups) < N_DATASETS:
        sys.exit(f"ERROR: need {N_DATASETS} distinct subtypes, only {len(groups)} qualify. "
                 f"Lower --min-jaccard.")
    print()

    out_root = Path(out_root)
    seeds: dict[int, int] = {}
    for i in indices(only):
        seeds[i] = write_subtype_dataset(
            out_root / f"Dataset-{i}", groups[i - 1], complete, order,
            dataset_rng(seed, i), reads, cpus, dry)

    write_seeds(out_root, seed, seeds)
    print("done:", out_root)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="mode", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--reads", default="1M", help="total reads (ISS -n); 1M = 500k pairs")
    common.add_argument("--cpus", type=int, default=_iss_int("cpus"),
                        help="(default: params.toml [iss].cpus)")
    common.add_argument("--only", type=int, default=None, metavar="I",
                        help="write only Dataset-I (identical to its slot in a full run)")

    p = sub.add_parser("syn-iso", parents=[common], help="1 chimeric genome, uniform abundance")
    p.add_argument("out_root", help="split root, e.g. samples/synthetic/isolate")
    p.add_argument("--seed", type=int, default=_iss_int("seed_syn_iso"),
                   help="parent seed (default: params.toml [iss].seed_syn_iso)")
    p.set_defaults(mode="syn-iso")

    p = sub.add_parser("syn-mix", parents=[common], help="3 chimeric genomes, Dirichlet shares")
    p.add_argument("out_root", help="split root, e.g. samples/synthetic/mixed")
    p.add_argument("--seed", type=int, default=_iss_int("seed_syn_mix"),
                   help="parent seed (default: params.toml [iss].seed_syn_mix)")
    p.set_defaults(mode="syn-mix")

    p = sub.add_parser("syn-mix-subtype", parents=[common],
                       help="3 real same-subtype strains per dataset, similarity-filtered")
    p.add_argument("out_root", help="e.g. samples/synthetic/mixed-subtype")
    p.add_argument("--min-jaccard", type=float, default=0.7)
    p.add_argument("--seed", type=int, default=_iss_int("seed_syn_mix_subtype"),
                   help="parent seed (default: params.toml [iss].seed_syn_mix_subtype)")
    p.add_argument("--dry-run", action="store_true",
                   help="write sources and manifests but skip read simulation")
    p.set_defaults(mode="syn-mix-subtype")

    args = ap.parse_args()
    if args.mode == "syn-iso":
        isolate(Path(args.out_root), args.seed, args.reads, args.cpus, args.only)
    elif args.mode == "syn-mix":
        mixed(Path(args.out_root), args.seed, args.reads, args.cpus, args.only)
    else:
        subtype(Path(args.out_root), args.seed, args.reads, args.cpus,
                args.min_jaccard, args.dry_run, args.only)
    return 0


if __name__ == "__main__":
    sys.exit(main())
