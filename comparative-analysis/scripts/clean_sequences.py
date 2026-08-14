#!/usr/bin/env python3
"""Produce indexes/sequences-cleaned.fasta

A record is dropped when it is:

  duplicate   byte-identical to a retained record
  substring   a contiguous substring of a longer retained record

Usage:
    python3 scripts/clean_sequences.py                        # indexes/sequences{,-cleaned}.fasta
    python3 scripts/clean_sequences.py --in X.fasta --out Y.fasta
"""
from __future__ import annotations

import argparse
import sys
from bisect import bisect_right
from collections import defaultdict
from pathlib import Path

import fastapy
import utils

SEP = "$"


def read_fasta(path: Path):
    """Returns (accessions, full header lines, uppercased sequences), in file order.
    """
    ids, headers, seqs = [], [], []
    for rec in fastapy.parse(path):
        ids.append(rec.id)
        headers.append(rec.description[1:])   # description is ">id desc"; headers omit the '>'
        seqs.append(rec.seq.upper())
    return ids, headers, seqs


def representative(idxs, ids, seqs):
    """Pick the surviving member of a group of identical records, order-independently."""
    def key(i):
        return (0 if ids[i].startswith("NC_") else 1, -len(seqs[i]), ids[i])
    return min(idxs, key=key)


def clean(src: Path, dst: Path, dropped: Path) -> dict[str, int]:
    """Drop duplicate and substring records; write the cleaned FASTA and the drop manifest.

    Returns {"read": n, "kept": k, "dropped": d}.
    """
    src, dst, dropped = Path(src), Path(dst), Path(dropped)
    ids, headers, seqs = read_fasta(src)
    n = len(ids)
    print(f"read {n} records ({sum(map(len, seqs)):,} bases) from {src}")

    reason: dict[int, tuple[str, int]] = {}

    groups: dict[str, list[int]] = defaultdict(list)
    for i, s in enumerate(seqs):
        groups[s].append(i)
    for g in groups.values():
        if len(g) > 1:
            keep = representative(g, ids, seqs)
            for i in g:
                if i != keep:
                    reason[i] = ("duplicate", keep)
    print(f"  duplicates : {sum(1 for r in reason.values() if r[0] == 'duplicate')} dropped "
          f"({sum(1 for g in groups.values() if len(g) > 1)} groups)")

    big = SEP.join(seqs)
    starts, off = [], 0
    for s in seqs:
        starts.append(off)
        off += len(s) + len(SEP)
    owner = lambda p: bisect_right(starts, p) - 1

    longest = max(map(len, seqs)) if seqs else 0
    for i, s in enumerate(seqs):
        if i in reason or not s or len(s) == longest:
            continue
        start = 0
        while True:
            p = big.find(s, start)
            if p < 0:
                break
            j = owner(p)
            if j != i and len(seqs[j]) > len(s):
                reason[i] = ("substring", j)
                break
            start = p + 1
    print(f"  substrings : {sum(1 for r in reason.values() if r[0] == 'substring')} dropped")

    kept = [i for i in range(n) if i not in reason]
    with open(dst, "w") as w:
        for i in kept:
            w.write(f">{headers[i]}\n")
            s = seqs[i]
            w.writelines(s[k:k + 60] + "\n" for k in range(0, len(s), 60))

    with open(dropped, "w") as w:
        w.write("dropped_accession\treason\tsuperseded_by\tdropped_len\tkept_len\n")
        for i in sorted(reason, key=lambda i: ids[i]):
            why, j = reason[i]
            w.write(f"{ids[i]}\t{why}\t{ids[j]}\t{len(seqs[i])}\t{len(seqs[j])}\n")

    print(f"\nkept    {len(kept)} / {n} records ({100 * len(kept) / n:.1f}%) -> {dst}")
    print(f"dropped {len(reason)} records, manifest -> {dropped}")
    if len(kept) + len(reason) != n:
        raise ValueError(f"kept ({len(kept)}) + dropped ({len(reason)}) != total ({n})")
    return {"read": n, "kept": len(kept), "dropped": len(reason)}


def main() -> int:
    idx = utils.indexes()
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="src", default=str(idx / "sequences.fasta"))
    ap.add_argument("--out", dest="dst", default=str(idx / "sequences-cleaned.fasta"))
    ap.add_argument("--dropped", default=str(idx / "sequences-dropped.tsv"))
    args = ap.parse_args()
    try:
        clean(Path(args.src), Path(args.dst), Path(args.dropped))
    except (ValueError, OSError) as e:
        print(f"clean_sequences: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
