#!/usr/bin/env python3
"""utils.py — shared helpers: where the benchmark's data lives, and the small
functions more than one script needs.

Two different roots, which earlier versions of these scripts conflated:

    code_root()   the checkout: run-analysis.py, params.toml, scripts/
    data_root()   $BENCH_DATA: indexes/, samples/, results/, truth/, ablation/
"""
from __future__ import annotations

import gzip
import os
import sys
from pathlib import Path
from typing import NoReturn

MARKERS = ("run-analysis.py", "run-all.sh", "run_ablation.sh")
PARAMS = "params.toml"

BENCH_DATA_ENV = "BENCH_DATA"


def code_root(start: Path | None = None) -> Path:
    """The benchmark code root
    """
    start = (start or Path(__file__)).resolve()
    for p in [start, *start.parents]:
        if (p / PARAMS).exists() and any((p / m).exists() for m in MARKERS):
            return p
    raise SystemExit(f"benchmark code root (the directory holding {PARAMS} and one of "
                     f"{'/'.join(MARKERS)}) not found above {start}")


def data_root_str() -> str:
    """The RAW $BENCH_DATA string, unnormalised, falling back to the code root.
    """
    return os.environ.get(BENCH_DATA_ENV) or str(code_root())


def data_root() -> Path:
    """The data root: $BENCH_DATA, or the code root when it is unset."""
    return Path(data_root_str())


def indexes() -> Path:
    return data_root() / "indexes"


def samples() -> Path:
    return data_root() / "samples"


def results() -> Path:
    return data_root() / "results"


def die(code: int, *lines: str) -> NoReturn:
    """Print to stderr and exit with `code`.
    """
    for line in lines:
        print(line, file=sys.stderr, flush=True)
    raise SystemExit(code)



def nonempty(p: Path) -> bool:
    """`[ -s "$p" ]` — exists AND non-zero. Path.exists() alone is not the same test."""
    p = Path(p)
    return p.is_file() and p.stat().st_size > 0


def open_maybe_gz(path):
    """Open `path` for text reading, decompressing iff its name ends `.gz`.
    """
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def count_pairs(path: Path) -> int:
    """Read pairs in a FASTQ: integer division, so a truncated file floors silently rather
    than reporting a fractional pair."""
    with open_maybe_gz(path) as f:
        return sum(1 for _ in f) // 4



def strip_version(acc: str | None) -> str | None:
    """`gsub('[.].*', '', acc)`: drop the FIRST '.' and everything after it.
    """
    if acc is None:
        return None
    i = acc.find(".")
    return acc[:i] if i >= 0 else acc


def strain_of(ref: str | None) -> str | None:
    """Collapse a predicted segment/accession to its influenza strain, for real-mix scoring.
    """
    if not ref:
        return None
    if ref.startswith("PR8"):
        return "PR8"
    if ref.startswith("WSN33"):
        return "WSN33"
    return "other"


def parse_timemem(path: Path) -> tuple[float | None, float | None]:
    """(wall_seconds, rss_gb) from a `/usr/bin/time -v` file.
    """
    if not Path(path).exists():
        return None, None
    wall = rss = None
    for line in Path(path).read_text(errors="ignore").splitlines():
        if "Elapsed (wall clock)" in line:
            parts = line.split(": ")[-1].strip().replace("s", "").split(":")
            try:
                wall = float(parts[-1])
                if len(parts) >= 2:
                    wall += float(parts[-2]) * 60
                if len(parts) >= 3:
                    wall += float(parts[-3]) * 3600
            except ValueError:
                pass
        elif "Maximum resident set size" in line:
            try:
                rss = int(line.split(": ")[-1].strip()) / 1048576.0   # kB -> GB
            except ValueError:
                pass
    return wall, rss
