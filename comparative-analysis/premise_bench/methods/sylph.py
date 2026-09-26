"""Sylph runner"""
from __future__ import annotations

import os
from pathlib import Path

from ..runner import CLEANED, Ctx, Job, fresh, normalize_output, run_cmd, timed
from ._io import _load_props

def contig_name_to_seqid(contig_name: str) -> str:
    """Extract accession from a Contig_name FASTA header"""
    return contig_name.strip().split()[0].split('|')[0].strip()


def genome_file_to_seqid(genome_file: str) -> str:
    """Extract accession from genome file path."""
    basename = os.path.basename(genome_file)
    for ext in ('.fasta', '.fa', '.fna'):
        if basename.endswith(ext):
            basename = basename[: -len(ext)]
            break
    return basename


def props_path(out: str | Path) -> str:
    """Where the companion proportions file goes, derived from the output name."""
    out = str(out)
    return out.replace('.norm.tsv', '.props') if out.endswith('.norm.tsv') else out + '.props'


def normalize(inp: Path, out: Path) -> int:
    """sylph profile TSV"""
    abundance: dict[str, float] = {}
    order: list[str] = []

    with open(inp) as fin:
        header = fin.readline().rstrip('\n').split('\t')
        try:
            genome_col = header.index('Genome_file')
            seq_abund_col = header.index('Sequence_abundance')
        except ValueError as e:
            raise ValueError(f'unexpected sylph header: {e}\nHeader: {header}') from e
        contig_col = header.index('Contig_name') if 'Contig_name' in header else None

        for line in fin:
            parts = line.rstrip('\n').split('\t')
            if len(parts) <= max(genome_col, seq_abund_col):
                continue
            seqid = ''
            if contig_col is not None and len(parts) > contig_col and parts[contig_col].strip():
                seqid = contig_name_to_seqid(parts[contig_col])
            if not seqid:
                seqid = genome_file_to_seqid(parts[genome_col])
            if not seqid:
                continue
            try:
                proportion = float(parts[seq_abund_col]) / 100.0
            except ValueError:
                proportion = 0.0
            if seqid not in abundance:
                abundance[seqid] = 0.0
                order.append(seqid)
            abundance[seqid] += proportion

    # Write per-read-style assignment file (one row per detected genome, seqID as target)
    with open(out, 'w') as fout:
        fout.write('readID\tseqID\n')
        for seqid in order:
            fout.write(f'{seqid}\t{seqid}\n')

    # Write proportion file for pi/distance analysis
    total = sum(abundance.values())
    with open(props_path(out), 'w') as fout:
        fout.write('seqID\tproportion\n')
        for seqid in order:
            norm = abundance[seqid] / total if total > 0 else 0.0
            fout.write(f'{seqid}\t{norm:.8f}\n')
    return len(order)

def build(ctx: Ctx) -> list[int]:
    fresh(ctx, "indexes/sylph")
    return [run_cmd(ctx, ["sylph", "sketch", CLEANED, "-i", "-c", ctx.p("syl", "c"),
                          "-o", "indexes/sylph/sequences", "-t", ctx.threads],
                    err="indexes/sylph/build.log")]


def classify(ctx: Ctx, j: Job) -> int:
    log = f"{j.outdir}/{j.ob}.log"
    rc = timed(ctx, j.tm, ["sylph", "profile", "indexes/sylph/sequences.syldb", j.r1, j.r2,
                           "-t", ctx.threads, "-c", ctx.p("syl", "c"),
                           "--min-number-kmers", ctx.p("syl", "min_kmers")],
               out=f"{j.outdir}/{j.ob}.tsv", err=log)
    normalize_output(ctx, "sylph", normalize, f"{j.outdir}/{j.ob}.tsv",
               f"{j.outdir}/{j.ob}.norm.tsv", log)
    return rc

def load(d: Path, base: str):
    props = d / f"{base}.norm.tsv.props"
    if not props.exists():
        props = d / f"{base}.props"
    return None, _load_props(props)
