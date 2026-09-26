"""KMCP runner"""
from __future__ import annotations

import gzip
from collections import Counter
from pathlib import Path

from ..runner import CLEANED, Ctx, Job, fresh, run_cmd, timed
from ._io import _norm_abund, _safe_lines, _warn

def build(ctx: Ctx) -> list[int]:
    fresh(ctx, "indexes/kmcp")
    if not ctx.dry_run:
        (ctx.root / "indexes/kmcp/kmers").mkdir(parents=True, exist_ok=True)
    return [
        run_cmd(ctx, ["kmcp", "compute", "--by-seq", "-O", "indexes/kmcp/kmers/",
                      "-k", ctx.p("kmp", "k"), "-j", ctx.threads, "--force", CLEANED]),
        run_cmd(ctx, ["kmcp", "index", "-I", "indexes/kmcp/kmers/",
                      "-O", "indexes/kmcp/sequences.kmcp", "-f", ctx.p("kmp", "f"),
                      "-n", ctx.p("kmp", "n"), "-j", ctx.threads, "--force"]),
    ]


def classify(ctx: Ctx, j: Job) -> int:
    # Three choices:
    #   -w   --load-whole-db: read the index into contiguous memory instead of mmap'ing it.
    #        kmcp recommends it for small databases, and this one is ~3.5 MB (db_build.csv's
    #        size_after_bytes for kmcp), so the extra memory is free. Purely a speed knob.
    #   .gz  the raw match dump is ~18 GB/dataset; writing it uncompressed is what blows past the
    #        timeout. profile + reduce both read .gz transparently.
    dump = f"{j.outdir}/{j.ob}.tsv.gz"
    rc = timed(ctx, j.tm, ["kmcp", "search", "-w", "-j", ctx.kmcp_threads,
                           "-d", "indexes/kmcp/sequences.kmcp", "-1", j.r1, "-2", j.r2,
                           "-o", dump], err=f"{j.outdir}/{j.ob}.log")
    if rc != 0:
        print(f"    kmcp search failed (exit {rc}) — keeping partial dump as "
              f"{j.ob}.tsv.gz.partial, no profile")
        if not ctx.dry_run:
            src = ctx.root / dump
            if src.exists():
                src.replace(ctx.root / f"{dump}.partial")
        return rc
    profile = f"{j.outdir}/{j.ob}.profile"
    run_cmd(ctx, ["kmcp", "profile", "--level", ctx.p("kmp", "profile_level"),
                  "-m", ctx.p("kmp", "profile_m"), dump, "-o", profile],
            err=f"{j.outdir}/{j.ob}.log", err_append=True)
    kmcp_pop_prop(ctx, profile, f"{j.outdir}/{j.ob}.pop-prop")

    return rc


def kmcp_pop_prop(ctx: Ctx, profile: str, dest: str) -> None:
    if ctx.dry_run:
        print(f"   $ {profile} -> {dest}  (fields 1 and 8)")
        return
    src = ctx.root / profile
    with open(ctx.root / dest, "w") as o:
        if not src.exists():
            return
        with open(src) as f:
            for line in f:
                fields = line.split()
                f1 = fields[0] if len(fields) > 0 else ""
                f8 = fields[7] if len(fields) > 7 else ""
                o.write(f"{f1} {f8}\n")

def load(d: Path, base: str):
    prof = d / f"{base}.profile"
    abund = None
    if prof.exists() and prof.stat().st_size:
        c = Counter()
        with open(prof) as f:
            hdr = next(f).rstrip("\n").split("\t")
            try:
                ri, pi = hdr.index("ref"), hdr.index("percentage")
            except ValueError:
                ri = pi = None
            if ri is not None:
                for line in f:
                    p = line.rstrip("\n").split("\t")
                    if len(p) <= max(ri, pi) or p[0].startswith("#"):
                        continue
                    ref = p[ri].replace("sequences.id_", "")
                    try:
                        c[ref] += float(p[pi])
                    except ValueError:
                        pass
        abund = _norm_abund(c)
    gz = d / f"{base}.tsv.gz"
    plain = d / f"{base}.tsv"
    tsv = plain if (plain.exists() and plain.stat().st_size) else gz
    assign = None
    if tsv.exists() and tsv.stat().st_size:
        opener = gzip.open if tsv.suffix == ".gz" else open
        state = {"truncated": False}
        with opener(tsv, "rt") as f:
            best = {}
            try:
                first = f.readline()
            except (EOFError, OSError) as e:
                first = ""
                state["truncated"] = True
                _warn(f"{tsv}: truncated/corrupt stream ({type(e).__name__}: {e})")
            hdr = first.rstrip("\n").split("\t") if first else []
            try:
                qi = hdr.index("qCov")
            except ValueError:
                qi = None
            for line in _safe_lines(f, tsv, state):
                if line.startswith("#"):
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 6:
                    continue
                rid = p[0].split("/", 1)[0]
                ref = p[5].replace("sequences.id_", "") if len(p) > 5 else None
                try:
                    q = float(p[qi]) if qi is not None else 0.0
                except ValueError:
                    q = 0.0
                if rid not in best or q > best[rid][0]:
                    best[rid] = (q, ref)
        if state["truncated"]:
            _warn(f"kmcp result discarded for {tsv.parent.name}: search did not finish")
            assign = None
            abund = None
        else:
            assign = {r: v[1] for r, v in best.items()}
    return assign, abund
