"""Ganon runner"""
from __future__ import annotations

from pathlib import Path

from ..runner import CLEANED, Ctx, Job, fresh, normalize_output, run_cmd, timed
from ._io import load_2col

def strip_mate(read_id: str) -> str:
    """Drop a trailing /1 or /2 mate suffix."""
    if len(read_id) > 2 and read_id[-2] == '/' and read_id[-1] in ('1', '2'):
        return read_id[:-2]
    return read_id


def clean_id(read_id: str) -> str:
    """Bare read ID: first whitespace token, mate suffix removed."""
    return strip_mate(read_id.split()[0] if read_id else read_id)


def normalize(inp: Path, out: Path) -> int:
    """ganon .one -> 2-column TSV. Returns the number of assignment rows written.
    """
    rows = 0
    with open(inp) as fin, open(out, 'w') as fout:
        fout.write('readID\tseqID\n')
        for line in fin:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            read_id, target = clean_id(parts[0]), parts[1]
            fout.write(f'{read_id}\t{target}\n')
            rows += 1
    return rows

def build(ctx: Ctx) -> list[int]:
    fresh(ctx, "indexes/ganon")
    abs_fasta = str(ctx.root / CLEANED)
    if ctx.dry_run:
        print(f"   $ indexes/seqid2taxid.map -> indexes/ganon/input.tsv  ({abs_fasta}\\t$1\\t$2)")
    else:
        with open(ctx.root / "indexes/seqid2taxid.map") as f, \
             open(ctx.root / "indexes/ganon/input.tsv", "w") as o:
            for line in f:
                fields = line.split()
                f1 = fields[0] if len(fields) > 0 else ""
                f2 = fields[1] if len(fields) > 1 else ""
                o.write(f"{abs_fasta}\t{f1}\t{f2}\n")
    return [run_cmd(ctx, ["ganon", "build-custom", "--input-file", "indexes/ganon/input.tsv",
                          "--input-target", "sequence", "--taxonomy", "ncbi",
                          "--taxonomy-files", "indexes/nodes.dmp", "indexes/names.dmp",
                          "--db-prefix", "indexes/ganon/sequences", "--skip-genome-size",
                          "--threads", ctx.threads],
                    err="indexes/ganon/build.log")]


def classify(ctx: Ctx, j: Job) -> int:
    log = f"{j.outdir}/{j.ob}.log"
    rc = timed(ctx, j.tm, ["ganon", "classify", "-d", "indexes/ganon/sequences",
                           "-p", j.r1, j.r2, "-o", f"{j.outdir}/{j.ob}", "--output-one",
                           "--multiple-matches", "em", "--rel-cutoff",
                           ctx.p("gan", "rel_cutoff"), "--fpr-query", ctx.p("gan", "fpr_query"),
                           "-t", ctx.threads], err=log)
    normalize_output(ctx, "ganon", normalize, f"{j.outdir}/{j.ob}.one", f"{j.outdir}/{j.ob}.tsv", log)
    return rc

load = load_2col
