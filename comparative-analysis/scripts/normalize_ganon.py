#!/usr/bin/env python3
"""
Normalize ganon classify --output-one output to 2-col TSV: readID<TAB>seqID
Input:  .one file (read_id<TAB>target<TAB>kmers...)
Output: .tsv file with header: readID  seqID
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

NAME = "normalize_ganon"


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


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('-i', '--input', required=True, type=Path, help='Ganon .one file')
    p.add_argument('-o', '--output', required=True, type=Path, help='Output .tsv file')
    args = p.parse_args()
    try:
        normalize(args.input, args.output)
    except (ValueError, OSError) as e:
        print(f'{NAME}: {e}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
