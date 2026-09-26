"""Centrifuger runner"""
from __future__ import annotations

from collections import Counter
from pathlib import Path

from ..runner import CLEANED, Ctx, Job, fresh, run_cmd, timed
from ..utils import is_unclassified
from ._io import _norm_abund

def build(ctx: Ctx) -> list[int]:
    fresh(ctx, "indexes/centrifuger")
    return [run_cmd(ctx, ["centrifuger-build", "--name-table", "indexes/names.dmp",
                          "--taxonomy-tree", "indexes/nodes.dmp", "-r", CLEANED,
                          "--conversion-table", "indexes/seqid2taxid.map",
                          "-o", "indexes/centrifuger/sequences-custom"])]


def classify(ctx: Ctx, j: Job) -> int:
    return timed(ctx, j.tm, ["centrifuger", "-t", ctx.threads,
                             "-x", "indexes/centrifuger/sequences-custom", "-1", j.r1, "-2", j.r2],
                 out=f"{j.outdir}/{j.ob}.tsv", err=f"{j.outdir}/{j.ob}.log")


def load(d: Path, base: str):
    tsv = d / f"{base}.tsv"
    if not tsv.exists() or not tsv.stat().st_size:
        return None, None
    assign = {}
    counts = Counter()
    with open(tsv) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2:
                continue
            ref = p[1]
            if is_unclassified(ref) or ref in ("no rank", "species"):
                ref = None
            assign[p[0]] = ref
            if ref:
                counts[ref] += 1
    return assign, _norm_abund(counts)
