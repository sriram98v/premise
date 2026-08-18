#!/usr/bin/env python3
"""
Normalize MORA output to a 2-col TSV: readID<TAB>seqID
Input:  MORA .txt file (query_id<TAB>reference_id)
Output: .tsv file with header: readID  seqID
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from utils import is_unclassified

NAME = "normalize_mora"


def normalize(inp: Path, out: Path) -> int:
    """MORA .txt -> 2-column TSV. Returns the number of assignment rows written.
    """
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


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('-i', '--input', required=True, type=Path, help='MORA .txt output file')
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
