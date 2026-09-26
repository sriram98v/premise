"""Karp runner"""
from __future__ import annotations

from collections import Counter
from pathlib import Path

from ..runner import CLEANED, Ctx, Job, fresh, normalize_output, run_cmd, timed
from ..utils import open_maybe_gz
from ._io import _norm_abund

def normalize(inp: Path, out: Path) -> int:
    """Karp _readinfo.txt[.gz] -> 2-column TSV, each read assigned its highest-probability reference"""
    rows = 0
    with open_maybe_gz(inp) as fin, open(out, 'w') as fout:
        fout.write('readID\tseqID\n')
        for line in fin:
            line = line.rstrip('\n')
            if not line or line.startswith('#') or line[0] in (' ', '\t'):
                continue
            parts = line.split()
            if len(parts) < 2:
                fout.write(f'{parts[0]}\tunclassified\n')
                rows += 1
                continue
            read_id = parts[0]
            best_pct = -1.0
            best_ref = 'unclassified'
            for token in parts[1:]:
                if ',' in token:
                    pct_str, ref_id = token.split(',', 1)
                elif ':' in token:
                    pct_str, ref_id = token.split(':', 1)
                else:
                    continue
                try:
                    pct = float(pct_str)
                except ValueError:
                    continue
                if pct > best_pct:
                    best_pct = pct
                    best_ref = ref_id
            fout.write(f'{read_id}\t{best_ref}\n')
            rows += 1
    return rows

def build(ctx: Ctx) -> list[int]:
    fresh(ctx, "indexes/karp")
    rcs = [run_cmd(ctx, ["samtools", "faidx", CLEANED])]
    if ctx.dry_run:
        print("   $ indexes/seqid2taxid.map -> indexes/sequences.tax  (Root;taxid_<N>; rows)")
    else:
        with open(ctx.root / "indexes/seqid2taxid.map") as f, \
             open(ctx.root / "indexes/sequences.tax", "w") as o:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) == 2:
                    o.write(parts[0] + "\tRoot;taxid_" + parts[1] + ";\n")
    rcs.append(run_cmd(ctx, ["karp", "-c", "index", "-r", CLEANED,
                             "-i", "indexes/karp/sequences.index", "-k", ctx.p("kap", "k")]))
    return rcs


def classify(ctx: Ctx, j: Job) -> int:
    log = f"{j.outdir}/{j.ob}.log"
    rc = timed(ctx, j.tm, ["karp", "-c", "quantify", "-r", CLEANED,
                           "-i", "indexes/karp/sequences.index", "-f", j.r1, "-q", j.r2,
                           "--paired", "-t", "indexes/sequences.tax", "-o", f"{j.outdir}/{j.ob}",
                           "--threads", ctx.threads, "--readinfo", "--no_harp_filter",
                           "--min_freq", ctx.p("kap", "min_freq")], err=log)
    normalize_output(ctx, "karp", normalize, f"{j.outdir}/{j.ob}_readinfo.txt.gz",
               f"{j.outdir}/{j.ob}.tsv", log)
    return rc

def load(d: Path, base: str):
    """Karp EM abundance from <base>.freqs (Label<TAB>ExpectedCounts<TAB>Taxa).
    """
    freqs = d / f"{base}.freqs"
    abund = None
    if freqs.exists() and freqs.stat().st_size:
        c = Counter()
        with open(freqs) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2:
                    try:
                        c[p[0]] = float(p[1])
                    except ValueError:
                        pass
        abund = _norm_abund(c)
    return None, abund
