#!/usr/bin/env python3
"""
Normalize Karp output to a 2-col TSV: readID<TAB>seqID
Input:  Karp .readinfo file: read_id  pct1,ref_id1  pct2,ref_id2 ...
        (space-separated; assign read to highest-probability ref_id)
Output: .tsv file with header: readID  seqID
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from utils import open_maybe_gz

NAME = "normalize_karp"


def normalize(inp: Path, out: Path) -> int:
    """Karp _readinfo.txt[.gz] -> 2-column TSV, each read assigned its highest-probability
    reference. Returns the number of assignment rows written.
    """
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


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('-i', '--input', required=True, type=Path,
                   help='Karp _readinfo.txt[.gz] file')
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
