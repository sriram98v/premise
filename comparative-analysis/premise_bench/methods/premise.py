"""PREMISE runner"""
from __future__ import annotations

from pathlib import Path

from ..runner import CLEANED, Ctx, Job, fresh, run_cmd, timed
from ..utils import is_unclassified
from ._io import _load_props

def build(ctx: Ctx) -> list[int]:
    fresh(ctx, "indexes/premise")
    return [run_cmd(ctx, ["premise", "build", "-s", CLEANED,
                          "-o", "indexes/premise/sequences.fmidx",
                          "--sa_sample_rate", ctx.p("pre", "sa_sample_rate")])]


def classify(ctx: Ctx, j: Job) -> int:
    iters = ctx.p("pre", "iters_real") if j.split.code == "real-iso" else ctx.p("pre", "iters")
    return timed(ctx, j.tm, ["premise", "query", "-s", "indexes/premise/sequences.fmidx",
                             "-t", ctx.p("pre", "threads"), "-1", j.r1, "-2", j.r2,
                             "-o", f"{j.outdir}/{j.ob}", "-m", ctx.p("pre", "mem"), "-i", iters,
                             "--eps_1", ctx.p("pre", "eps_1"), "--eps_2", ctx.p("pre", "eps_2"),
                             "--em_threshold", ctx.p("pre", "em_threshold"),
                             "--rho", ctx.p("pre", "rho"), "--omega", ctx.p("pre", "omega")],
                 out=f"{j.outdir}/{j.ob}.log", err_to_out=True)

def load(d: Path, base: str):
    matches = d / f"{base}.matches"
    assign = None
    if matches.exists() and matches.stat().st_size:
        assign = {}
        with open(matches) as f:
            next(f)
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < 2:
                    continue
                assign[p[0]] = None if is_unclassified(p[1]) else p[1]
    return assign, _load_props(d / f"{base}.props")
