#!/usr/bin/env python3
"""Fetch the real samples from NCBI SRA into samples/real/{isolate,mixed}/<SRR>/.

Usage:
    python3 scripts/fetch_real_samples.py                     # both splits, skipping what exists
    python3 scripts/fetch_real_samples.py --split real-iso    # real-iso runs only
    python3 scripts/fetch_real_samples.py --threads 16
    python3 scripts/fetch_real_samples.py --dry-run
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import utils

SPLITS: dict[str, tuple[str, tuple[str, ...]]] = {
    "real-iso": ("real/isolate", ("SRR31013463", "SRR31013465", "SRR31013467", "SRR31013473")),
    "real-mix": ("real/mixed", ("SRR3360139", "SRR3360140", "SRR3360145", "SRR3360146")),
}


def run(cmd: list[str], **kw) -> None:
    print("   $ " + " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, **kw)


def tool_version(tool: str) -> str:
    try:
        out = subprocess.run([tool, "--version"], capture_output=True, text=True, timeout=30)
        return (out.stdout + out.stderr).strip().splitlines()[-1].strip()
    except Exception:
        return "?"


def have_reads(dest: Path, srr: str) -> bool:
    """Returns True if both mates are present and non-empty."""
    mates = (dest / f"{srr}_1.fastq", dest / f"{srr}_2.fastq")
    return all(f.exists() and f.stat().st_size > 0 for f in mates)


def fetch(srr: str, dest: Path, threads: int, prefetch: bool, clean_sra: bool) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    sra = dest / f"{srr}.sra"

    if prefetch and not sra.exists():
        run(["prefetch", "--max-size", "u", "--output-file", str(sra), srr])

    src = str(sra) if sra.exists() else srr
    run(["fasterq-dump", "--split-files", "--threads", str(threads),
         "--outdir", str(dest), "--temp", str(dest), src])

    for mate in (1, 2):
        f = dest / f"{srr}_{mate}.fastq"
        if not f.exists() or f.stat().st_size == 0:
            raise SystemExit(f"{srr}: fasterq-dump produced no {f.name} — "
                             f"the run may not be paired-end")
    if clean_sra and sra.exists():
        sra.unlink()


def fetch_all(split: str = "all", accessions: list[str] | None = None, threads: int = 8,
              prefetch: bool = True, clean_sra: bool = False, force: bool = False,
              dry: bool = False) -> int:
    """Fetch the real SRA runs for `split` under $BENCH_DATA.
    """
    for tool in ("fasterq-dump",) + (() if not prefetch else ("prefetch",)):
        if shutil.which(tool) is None:
            raise ValueError(f"{tool} not on PATH — enter the pinned toolchain: "
                             f"nix develop ..#benchmark")

    if accessions:
        known = {s for _, accs in SPLITS.values() for s in accs}
        for a in accessions:
            if a not in known:
                raise ValueError(f"unknown accession: {a} (not one of {', '.join(sorted(known))})")

    wanted = sorted(SPLITS) if split == "all" else [split]
    jobs: list[tuple[str, Path]] = []
    for sp in wanted:
        subdir, accs = SPLITS[sp]
        for srr in accs:
            if accessions and srr not in accessions:
                continue
            jobs.append((srr, utils.samples() / subdir / srr))

    print(f"data root   : {utils.data_root()}")
    print(f"sra-tools   : {tool_version('fasterq-dump')}")
    print(f"runs        : {len(jobs)}\n")

    todo = [(srr, d) for srr, d in jobs if force or not have_reads(d, srr)]
    pending = {srr for srr, _ in todo}
    for srr, _dest in jobs:
        if srr not in pending:
            print(f"== {srr}: already present — skipping (--force to re-fetch)")
    if dry:
        for srr, dest in todo:
            print(f"== {srr}: would fetch -> {dest}")
        return 0

    for i, (srr, dest) in enumerate(todo, 1):
        print(f"== [{i}/{len(todo)}] {srr} -> {dest}")
        try:
            fetch(srr, dest, threads, prefetch, clean_sra)
        except subprocess.CalledProcessError as e:
            print(f"\nFETCH FAILED for {srr} (exit {e.returncode}). Re-run to resume; already "
                  f"completed runs are skipped.", file=sys.stderr)
            return 1

    print(f"\nfetched {len(todo)} run(s). Next: preprocess them —")
    print("  python3 scripts/prepare_real_samples.py   # trimmed + filtered reads, and their truth")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--split", choices=sorted(SPLITS) + ["all"], default="all")
    ap.add_argument("--accessions", nargs="*", default=None,
                    help="fetch only these accessions (must belong to the selected split)")
    ap.add_argument("--threads", type=int, default=int(os.environ.get("THREADS", "8")))
    ap.add_argument("--no-prefetch", action="store_true", help="stream instead of prefetching .sra")
    ap.add_argument("--clean-sra", action="store_true", help="delete each .sra after conversion")
    ap.add_argument("--force", action="store_true", help="re-fetch runs whose FASTQ already exists")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    try:
        return fetch_all(args.split, args.accessions, args.threads, not args.no_prefetch,
                         args.clean_sra, args.force, args.dry_run)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
