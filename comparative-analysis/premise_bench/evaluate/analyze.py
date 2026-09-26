"""Score every method on every dataset of the chosen splits and write one tidy CSV per split"""
from __future__ import annotations

import argparse
import sys

from ..splits import SPLITS
from .score import compute_dataset, gather_records, write_split_csv


def run(splits) -> str:
    """Score `splits` and write their CSVs; return a line per split written"""
    lines = []
    for split in [s for s in splits if s in SPLITS]:
        per_ds = {}
        for ds in SPLITS[split].datasets:
            try:
                per_ds[ds] = compute_dataset(split, ds)
            except (FileNotFoundError, OSError) as e:
                lines.append(f"  [skip {split}/{ds}: {e}]")
        if not per_ds:
            continue
        path = write_split_csv(split, gather_records(split, list(per_ds), per_ds))
        lines.append(f"[wrote {path}]")
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser(description="Score results -> results/tables/comparative-<split>.csv")
    ap.add_argument("--splits", default=",".join(SPLITS),
                    help=f"comma-separated subset of {', '.join(SPLITS)}")
    args = ap.parse_args()
    try:
        print(run(args.splits.split(",")))
    except (ValueError, OSError) as e:
        print(f"analyze: {e}", file=sys.stderr)
        return 1
    return 0
