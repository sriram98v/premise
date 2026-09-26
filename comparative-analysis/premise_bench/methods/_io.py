"""Parsing helpers shared by the per-method result loaders."""
from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

from ..utils import is_unclassified

def _warn(msg: str):
    print(f"analyze: WARNING: {msg}", file=sys.stderr)


def _safe_lines(f, path, state):
    """Yield lines from a (possibly gzip) handle, stopping at a truncated stream.
    """
    try:
        for line in f:
            yield line
    except (EOFError, OSError) as e:
        state["truncated"] = True
        _warn(f"{path}: truncated/corrupt stream ({type(e).__name__}: {e})")


def _norm_abund(counts: dict) -> dict:
    total = sum(counts.values())
    return {k: v / total for k, v in counts.items()} if total else {}


def _load_props(path: Path):
    """A seqID<TAB>proportion file -> normalised abundance, or None when absent or empty.
    """
    if not path.exists() or not path.stat().st_size:
        return None
    c = Counter()
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                try:
                    c[p[0]] = float(p[1])
                except ValueError:
                    pass
    return _norm_abund(c)


def _read_2col_tsv(path: Path, strip_mate=False):
    """readID<TAB>seqID with header -> {readID: stripped ref | None}."""
    out = {}
    if not path.exists() or path.stat().st_size == 0:
        return None
    with open(path) as f:
        next(f, None)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2 or not p[0]:
                continue
            rid = p[0].split("/", 1)[0] if strip_mate else p[0]
            ref = p[1].strip()
            out[rid] = None if is_unclassified(ref) else ref
    return out


def load_2col(d: Path, base: str, strip_mate=False):
    """mora/ganon: normalized 2-col .tsv; abundance from assignment counts."""
    assign = _read_2col_tsv(d / f"{base}.tsv", strip_mate=strip_mate)
    if assign is None:
        return None, None
    abund = _norm_abund(Counter(r for r in assign.values() if r is not None))
    return assign, abund
