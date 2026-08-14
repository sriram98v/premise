#!/usr/bin/env python3
"""
Normalize sylph profile output to two files:
  1. <output>  — readID<TAB>seqID  (per-read assignment proxy; one row per genome detected)
  2. <output>.props — seqID<TAB>proportion  (Sequence_abundance / 100, the read fraction)
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from utils import strip_version

NAME = "normalize_sylph"


def contig_name_to_seqid(contig_name: str) -> str:
    """Extract accession from a Contig_name FASTA header: 'KT002519.1 |Influenza A ...'."""
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


def normalize(inp: Path, out: Path, strip_versions: bool = True) -> int:
    """sylph profile TSV -> the assignment TSV and its companion .props. Returns the number of
    distinct genomes detected.
    """
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
            if strip_versions:
                seqid = strip_version(seqid)
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


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('-i', '--input', required=True, type=Path, help='Sylph profile TSV output')
    p.add_argument('-o', '--output', required=True, type=Path,
                   help='Output normalised TSV (readID\\tseqID per detected genome)')
    p.add_argument('--no-strip-version', action='store_true',
                   help='Keep accession version suffixes')
    args = p.parse_args()
    try:
        normalize(args.input, args.output, strip_versions=not args.no_strip_version)
    except (ValueError, OSError) as e:
        print(f'{NAME}: {e}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
