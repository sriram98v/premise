"""MORA runner"""
from __future__ import annotations

from pathlib import Path

from ..runner import CLEANED, Ctx, Job, fresh, normalize_output, run_cmd, timed
from ..utils import is_unclassified
from ._io import load_2col

def normalize(inp: Path, out: Path) -> int:
    """MORA .txt -> 2-column TSV"""
    rows = 0
    with open(inp) as fin, open(out, 'w') as fout:
        fout.write('readID\tseqID\n')
        header_skipped = False
        for line in fin:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            read_id, ref_id = parts[0], parts[1]
            
            if not header_skipped and read_id.lower() in ('read_id', 'query_id', 'readid'):
                header_skipped = True
                continue
            header_skipped = True

            # "NOT ALIGNED" is an unclassified read
            if is_unclassified(ref_id):
                ref_id = 'unclassified'
            fout.write(f'{read_id}\t{ref_id}\n')
            rows += 1
    return rows

def build(ctx: Ctx) -> list[int]:
    fresh(ctx, "indexes/mora")
    return [run_cmd(ctx, [ctx.bwa, "index", "-p", "indexes/mora/sequences", CLEANED])]


def classify(ctx: Ctx, j: Job) -> int:
    sam = f"{j.outdir}/{j.ob}.sam"
    rc = timed(ctx, j.tm, [ctx.bwa, "mem", "-t", ctx.threads, "indexes/mora/sequences",
                           j.r1, j.r2], out=sam, err=f"{j.outdir}/{j.ob}.bwa.log")
    log = f"{j.outdir}/{j.ob}.mora.log"
    run_cmd(ctx, ["mora", "-s", sam, "-o", f"{j.outdir}/{j.ob}.txt"], err=log)
    normalize_output(ctx, "mora", normalize, f"{j.outdir}/{j.ob}.txt", f"{j.outdir}/{j.ob}.tsv", log)
    if ctx.dry_run:
        print(f"   $ rm -f {sam}")
    else:
        (ctx.root / sam).unlink(missing_ok=True)
    return rc

load = load_2col
