#!/usr/bin/env python3
"""prepare_real_samples.py — turn the raw SRA reads into analysis-ready reads + per-read truth.

Both real splits run the same pipeline.

Usage:
  python3 scripts/prepare_real_samples.py
  python3 scripts/prepare_real_samples.py --split real-iso
  python3 scripts/prepare_real_samples.py --threads 16 --force
  python3 scripts/prepare_real_samples.py --samples SRR31013463 SRR31013465
  python3 scripts/prepare_real_samples.py --dry-run
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import NamedTuple

import read_truth
import utils
from utils import count_pairs, die, nonempty

CODE = utils.code_root()

def root_str() -> str:
    return utils.data_root_str()


def root() -> Path:
    return utils.data_root()


def tmpidx_str() -> str:
    return os.environ.get("TMPIDX") or f"{root_str()}/.tmp-prepidx"

SRC_FASTA = "true_sources.fasta"

TRIM_SFX = ("_1.ca.fastq", "_2.ca.fastq")
OUT_SFX = ("_1-filtered.ca.fastq", "_2-filtered.ca.fastq")
TRUTH_TSV = read_truth.TRUTH_NAME
UNCL_TXT = read_truth.UNCL_NAME

BWA = "bwa-mem2"

ISO_PRIMER_5 = "GACCATCTAGCGACCTCCACNNNNNNNN"
MIX_ADAPT = "AGATCGGAAGAGC"


@dataclass(frozen=True)
class Split:
    """Everything that differs between the two real splits.
    """
    name: str
    sub: str
    cutadapt: tuple[str, ...]


SPLITS: dict[str, Split] = {
    "real-iso": Split(
        name="real-iso",
        sub="samples/real/isolate",
        cutadapt=("-g", ISO_PRIMER_5, "-G", ISO_PRIMER_5,
                  "--overlap", "8", "--minimum-length", "20", "--trim-n"),
    ),
    "real-mix": Split(
        name="real-mix",
        sub="samples/real/mixed",
        cutadapt=("-a", MIX_ADAPT, "-A", MIX_ADAPT, "-m", "30"),
    ),
}

ORDER = ("real-iso", "real-mix")

class Sample(NamedTuple):
    base: str
    dir: Path
    raw1: Path
    raw2: Path
    tr1: Path
    tr2: Path
    out1: Path
    out2: Path
    srcfa: Path


def ts() -> str:
    """The `date '+%F %T'` analogue, spelled out (%F/%T are glibc extensions)."""
    return time.strftime("%Y-%m-%d %H:%M:%S")

def say(msg: str) -> None:
    print(f"[{ts()}] {msg}", flush=True)


def warn(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def rel(p: Path) -> str:
    """Path as the shell would have it: relative to $ROOT, which every child sees as cwd."""
    return str(p.relative_to(root()))


def fasta_headers(path: Path) -> tuple[int, list[str]]:
    """Returns (count of '>' lines, first whitespace field of each, in file order).
    """
    n, names = 0, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                n += 1
                names.append((line[1:].split() or [""])[0])
    return n, names


def check_tools(dry: bool) -> str:
    """Collect every missing tool before reporting, as the shell does, and resolve bwa-mem2."""
    missing = [t for t in ("bwa-mem2", "cutadapt") if shutil.which(t) is None]
    if missing and not dry:
        die(2, f"ERROR: not on PATH: {' '.join(missing)} — enter the pinned toolchain first: "
               f"nix develop ..#benchmark")
    return shutil.which("bwa-mem2") or "bwa-mem2"


def echo(cmd: list[str]) -> None:
    print("   $ " + " ".join(cmd), flush=True)


def ensure_index(src: Path, prefix: str, dry: bool) -> None:
    Path(prefix).parent.mkdir(parents=True, exist_ok=True)
    if Path(f"{prefix}.bwt.2bit.64").is_file():
        return
    say(f"indexing {rel(src)} -> {prefix}")
    cmd = [BWA, "index", "-p", prefix, rel(src)]
    if dry:
        echo(cmd)
        return
    try:
        subprocess.run(cmd, cwd=root_str(), check=True)
    except subprocess.CalledProcessError as e:
        die(e.returncode, f"ERROR: bwa-mem2 index exited {e.returncode} for {rel(src)}")


def run_cutadapt(spec: Split, s: Sample, threads: int, dry: bool) -> None:
    cmd = (["cutadapt", "-j", str(threads), *spec.cutadapt,
            "-o", rel(s.tr1), "-p", rel(s.tr2), rel(s.raw1), rel(s.raw2)])
    if dry:
        echo(cmd)
        return
    with open(s.dir / f"{s.base}.cutadapt.log", "w") as log:
        try:
            subprocess.run(cmd, cwd=root_str(), stdout=log, check=True)
        except subprocess.CalledProcessError as e:
            die(e.returncode, f"ERROR: {s.base}: cutadapt exited {e.returncode} — "
                              f"see {rel(s.dir)}/{s.base}.cutadapt.log")

def align_and_derive_truth(s: Sample, idx: str, threads: int) -> read_truth.TruthResult:
    cmd = [BWA, "mem", "-t", str(threads), idx, rel(s.tr1), rel(s.tr2)]
    with open(s.dir / f"{s.base}.bwa.log", "w") as errlog:
        proc = subprocess.Popen(cmd, cwd=root_str(), stdout=subprocess.PIPE, stderr=errlog)
        try:
            result = read_truth.run(proc.stdout, s.dir, collect_ids=True)
        except BaseException:
            proc.stdout.close()
            proc.kill()
            proc.wait()
            raise
        proc.stdout.close()
        rc = proc.wait()
    if rc != 0:
        die(rc, f"ERROR: {s.base}: bwa-mem2 mem exited {rc} — "
                f"see {rel(s.dir)}/{s.base}.bwa.log")
    return result


def filter_fastq(inp: Path, out: Path, names: set[str], invert: bool = False) -> tuple[int, int]:
    """Copy records from `inp` to `out`, dropping (or with invert=True, keeping) listed IDs.
    """
    kept = dropped = 0
    with open(inp) as fi, open(out, "w") as fo:
        while True:
            header = fi.readline()
            if not header:
                break
            seq, _plus, qual = fi.readline(), fi.readline(), fi.readline()
            if not qual:
                die(1, f"ERROR: {inp} ends mid-record (truncated FASTQ)")
            if not header.startswith("@"):
                die(1, f"ERROR: {inp}: expected '@' header, got {header[:40]!r}")
            if not _plus.startswith("+"):
                die(1, f"ERROR: {inp}: expected '+' separator, got {_plus[:40]!r}")
            read_id = header[1:].split(None, 1)[0]
            if (read_id in names) != invert:
                dropped += 1
            else:
                fo.write(header); fo.write(seq); fo.write("+\n"); fo.write(qual)
                kept += 1
    return kept, dropped


def subset_reads(s: Sample, keep: set[str]) -> tuple[int, int]:
    """Returns (kept, dropped) for R1.
    """
    counts = []
    for inp, out in ((s.tr1, s.out1), (s.tr2, s.out2)):
        kept, dropped = filter_fastq(inp, out, keep, invert=True)
        print(f"fastq_exclude: {inp.name}: kept {kept}, dropped {dropped} "
              f"({len(keep)} names in keep-list)", file=sys.stderr, flush=True)
        counts.append((kept, dropped))
    if counts[0][0] != counts[1][0]:
        die(1, f"ERROR: {s.base}: mates out of step — R1 kept {counts[0][0]}, "
               f"R2 kept {counts[1][0]}")
    return counts[0]


def prepare_sample(spec: Split, s: Sample, threads: int, force: bool, dry: bool) -> None:
    if not force and nonempty(s.out1) and nonempty(s.out2) and nonempty(s.dir / TRUTH_TSV):
        say(f"{s.base}: already prepared (--force to redo)")
        return
    if not (nonempty(s.raw1) and nonempty(s.raw2)):
        say(f"{s.base}: raw reads missing — run scripts/fetch_real_samples.py first, skip")
        return

    if not nonempty(s.srcfa):
        die(2, f"ERROR: {root_str()}/{rel(s.srcfa)} not found.",
               "  It is the reference set this sample's truth is derived from, and truth cannot",
               "  be built without it: it lists exactly the references this sample contains.")

    idx = f"{tmpidx_str()}/idx-{s.base}"

    say(f"{s.base} [1/3] cutadapt trim")
    run_cutadapt(spec, s, threads, dry)

    say(f"{s.base} [2/3] align to {SRC_FASTA} + derive truth")
    ensure_index(s.srcfa, idx, dry)
    if dry:
        echo([BWA, "mem", "-t", str(threads), idx, rel(s.tr1), rel(s.tr2)])
        return
    result = align_and_derive_truth(s, idx, threads)

    declared, names = fasta_headers(s.srcfa)
    seen = len(result.refs)
    if seen != declared:
        warn(f"  WARN: {s.base}: {seen}/{declared} declared references carry reads —")
        for n in sorted(names):
            if n not in result.refs:
                warn(f"    no reads: {n}")

    say(f"{s.base} [3/3] subset reads to the truth set")
    assert result.read_ids is not None
    got, dropped = subset_reads(s, result.read_ids)
    trimmed = got + dropped

    if got != result.kept:
        die(1, f"ERROR: {s.base}: {got} read pairs but {result.kept} truth rows")
    if result.pairs != trimmed:
        die(1, f"ERROR: {s.base}: {result.pairs - result.kept} unclassified but "
               f"{trimmed - got} trimmed pairs lack truth")

    uncl = result.pairs - result.kept
    pct = float("%.6g" % (100.0 * got / trimmed)) if trimmed else 0.0
    say(f"{s.base} DONE  raw={count_pairs(s.raw1)}  trimmed={trimmed}  refs={declared}  "
        f"final={got}  unclassified={uncl} ({pct:.2f}% of trimmed kept)")


def prepare_split(spec: Split, threads: int, force: bool, only: list[str], dry: bool) -> None:
    sdir = root() / spec.sub
    for d in sorted(sdir.glob("SRR*")):
        if not d.is_dir():
            continue
        base = d.name
        if only and base not in only:
            continue
        prepare_sample(spec, Sample(
            base=base, dir=d,
            raw1=d / f"{base}_1.fastq", raw2=d / f"{base}_2.fastq",
            tr1=d / f"{base}{TRIM_SFX[0]}", tr2=d / f"{base}{TRIM_SFX[1]}",
            out1=d / f"{base}{OUT_SFX[0]}", out2=d / f"{base}{OUT_SFX[1]}",
            srcfa=d / SRC_FASTA,
        ), threads, force, dry)


def run(split: str = "all", threads: int = 8, only: list[str] | None = None,
        force: bool = False, dry: bool = False) -> int:
    only = only or []
    if not root().is_dir():
        die(2, f"ERROR: BENCH_DATA does not exist: {root_str()}")

    global BWA
    BWA = check_tools(dry)
    if not dry:
        Path(tmpidx_str()).mkdir(parents=True, exist_ok=True)

    print(f"data root : {root_str()}", flush=True)
    print(f"threads   : {threads}", flush=True)
    if only:
        print(f"samples   : {' '.join(only)}", flush=True)

    for name in ORDER:
        if split in (name, "all"):
            prepare_split(SPLITS[name], threads, force, only, dry)

    say("ALL DONE")
    if split != "real-mix":
        print(f"  real-iso -> samples/real/isolate/<SRR>/{{<SRR>_{{1,2}}-filtered.ca.fastq,"
              f"{TRUTH_TSV},{UNCL_TXT}}}", flush=True)
    if split != "real-iso":
        print(f"  real-mix -> samples/real/mixed/<SRR>/{{<SRR>_{{1,2}}-filtered.ca.fastq,"
              f"{TRUTH_TSV},{UNCL_TXT}}}", flush=True)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--split", choices=["real-iso", "real-mix", "all"], default="all")
    ap.add_argument("--threads", type=int, default=int(os.environ.get("THREADS", "8")))
    ap.add_argument("--samples", nargs="*", default=None,
                    help="only these samples (also honours $SAMPLES, whitespace-separated)")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    raw = args.samples if args.samples is not None else [os.environ.get("SAMPLES", "")]
    only = [s for chunk in raw for s in chunk.split()]
    force = args.force or os.environ.get("FORCE") == "1"
    return run(args.split, args.threads, only, force, args.dry_run)


if __name__ == "__main__":
    sys.exit(main())
