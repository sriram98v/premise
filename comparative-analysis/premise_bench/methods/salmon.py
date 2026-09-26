"""Salmon runner"""
from __future__ import annotations

from pathlib import Path

from ..runner import CLEANED, Ctx, Job, fresh, normalize_output, run_cmd, timed
from ..utils import open_maybe_gz
from ._io import _load_props, _read_2col_tsv

_UNMAPPED = 0x4


_SECONDARY = 0x100


_SUPPLEMENTARY = 0x800


def props_path(out: str | Path) -> str:
    """Where the companion proportions file goes, derived from the output name."""
    out = str(out)
    return out[:-len('.tsv')] + '.props' if out.endswith('.tsv') else out + '.props'


def _alignment_score(fields: list[str]) -> float | None:
    """The AS:i: value of one SAM record, or None when the record carries no AS tag."""
    for tag in fields[11:]:
        if tag.startswith('AS:i:'):
            try:
                return float(tag[5:])
            except ValueError:
                return None
    return None


def normalize(inp: Path, out: Path) -> int:
    """salmon `--writeMappings` SAM -> 2-column TSV, one row per mapped read"""
    best: dict[str, tuple[float, int, str]] = {}
    order: list[str] = []

    with open_maybe_gz(inp) as fin:
        for line in fin:
            if not line or line[0] == '@':
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 3:
                continue
            qname, flag_str, rname = fields[0], fields[1], fields[2]
            try:
                flag = int(flag_str)
            except ValueError:
                continue
            if flag & _UNMAPPED or rname == '*':
                continue
            score = _alignment_score(fields)
            rank = (score if score is not None else float('-inf'),
                    0 if flag & (_SECONDARY | _SUPPLEMENTARY) else 1)
            prev = best.get(qname)
            if prev is None:
                order.append(qname)
            elif rank <= (prev[0], prev[1]):
                continue
            best[qname] = (rank[0], rank[1], rname)

    with open(out, 'w') as fout:
        fout.write('readID\tseqID\n')
        for qname in order:
            fout.write(f'{qname}\t{best[qname][2]}\n')
    return len(order)


def props(inp: Path, out: Path) -> int:
    """salmon quant.sf -> seqID<TAB>proportion, the EM read fraction per target"""
    counts: dict[str, float] = {}
    order: list[str] = []

    with open(inp) as fin:
        header = fin.readline().rstrip('\n').split('\t')
        try:
            name_col = header.index('Name')
            reads_col = header.index('NumReads')
        except ValueError as e:
            raise ValueError(f'unexpected salmon quant.sf header: {e}\nHeader: {header}') from e

        for line in fin:
            parts = line.rstrip('\n').split('\t')
            if len(parts) <= max(name_col, reads_col):
                continue
            seqid = parts[name_col].strip().split()[0] if parts[name_col].strip() else ''
            if not seqid:
                continue
            try:
                n = float(parts[reads_col])
            except ValueError:
                n = 0.0
            if seqid not in counts:
                order.append(seqid)
            counts[seqid] = n

    total = sum(counts.values())
    with open(out, 'w') as fout:
        fout.write('seqID\tproportion\n')
        for seqid in order:
            norm = counts[seqid] / total if total > 0 else 0.0
            fout.write(f'{seqid}\t{norm:.8f}\n')
    return len(order)

def build(ctx: Ctx) -> list[int]:
    fresh(ctx, "indexes/salmon")
    return [run_cmd(ctx, ["salmon", "index", "-t", CLEANED, "-i", "indexes/salmon/sequences",
                          "-k", ctx.p("sal", "k"), "--keepDuplicates", "-p", ctx.threads],
                    err="indexes/salmon/build.log")]


def classify(ctx: Ctx, j: Job) -> int:
    """salmon output"""
    log = f"{j.outdir}/{j.ob}.log"
    sam = f"{j.outdir}/{j.ob}.sam"
    quant = f"{j.outdir}/{j.ob}.salmon"
    rc = timed(ctx, j.tm, ["salmon", "quant", "-i", "indexes/salmon/sequences",
                           "-l", ctx.p("sal", "libtype"), "-1", j.r1, "-2", j.r2,
                           "-p", ctx.threads, f"--writeMappings={sam}", "-o", quant],
               err=log)
    normalize_output(ctx, "salmon", normalize, sam, f"{j.outdir}/{j.ob}.tsv", log)
    normalize_output(ctx, "salmon_props", props, f"{quant}/quant.sf", f"{j.outdir}/{j.ob}.props", log)
    if ctx.dry_run:
        print(f"   $ rm -f {sam}")
    else:
        (ctx.root / sam).unlink(missing_ok=True)
    return rc

def load(d: Path, base: str):
    """Per-read labels from the mappings"""
    assign = _read_2col_tsv(d / f"{base}.tsv")
    return assign, _load_props(d / f"{base}.props")
