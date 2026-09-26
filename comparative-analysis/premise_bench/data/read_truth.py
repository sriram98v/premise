"""Derive per-read ground truth for one real sample"""
from __future__ import annotations

import argparse
import os
import sys
from dataclasses import dataclass
from pathlib import Path

import pysam

from ..utils import die

TRUTH_NAME = "truth_assignments.tsv"
UNCL_NAME = "unclassified_reads.txt"

DROP_REASONS = ("chimeric", "none_mapped", "half_mapped", "cross_ref")


def is_chimeric(aln) -> bool:
    return (aln.is_supplementary
            or aln.has_tag("SA")
            or "H" in (aln.cigarstring or ""))


def stream_groups(alignments):
    cur, buf = None, []
    for aln in alignments:
        if aln.query_name != cur:
            if buf:
                yield buf
            cur, buf = aln.query_name, [aln]
        else:
            buf.append(aln)
    if buf:
        yield buf


@dataclass(frozen=True)
class TruthResult:
    pairs: int
    kept: int
    drops: dict[str, int]
    refs: set[str]
    read_ids: set[str] | None


def resolve_pair(recs: list) -> tuple[str | None, str | None]:
    if any(is_chimeric(r) for r in recs):
        return None, "chimeric"
    mapped = {(1 if r.is_read1 else 2): r for r in recs
              if not r.is_secondary and not r.is_unmapped}
    if len(mapped) < 2:
        return None, ("none_mapped" if not mapped else "half_mapped")
    refs = {r.reference_name for r in mapped.values()}
    if len(refs) != 1:
        return None, "cross_ref"
    return refs.pop(), None


def run(source, sample_dir: Path, collect_ids: bool = False) -> TruthResult:
    """Derive truth from a name-grouped SAM stream.
    """
    truth_out, uncl_out = sample_dir / TRUTH_NAME, sample_dir / UNCL_NAME
    truth_tmp, uncl_tmp = truth_out.with_suffix(truth_out.suffix + ".tmp"), \
        uncl_out.with_suffix(uncl_out.suffix + ".tmp")

    pairs = kept = 0
    drops = dict.fromkeys(DROP_REASONS, 0)
    refs: set[str] = set()
    ids: set[str] | None = set() if collect_ids else None

    try:
        with open(truth_tmp, "w") as fh_truth, open(uncl_tmp, "w") as fh_uncl, \
                pysam.AlignmentFile(source, "r") as af:
            for recs in stream_groups(af):
                pairs += 1
                ref, reason = resolve_pair(recs)
                if ref is None:
                    drops[reason] += 1
                    fh_uncl.write(f"{recs[0].query_name}\n")
                    continue
                fh_truth.write(f"{recs[0].query_name}\t{ref}\n")
                refs.add(ref)
                if ids is not None:
                    ids.add(recs[0].query_name)
                kept += 1
        os.replace(truth_tmp, truth_out)
        os.replace(uncl_tmp, uncl_out)
    except BaseException:
        for tmp in (truth_tmp, uncl_tmp):
            tmp.unlink(missing_ok=True)
        raise

    sys.stderr.write(
        f"pairs={pairs} kept={kept} "
        f"dropped(chimeric={drops['chimeric']} none_mapped={drops['none_mapped']} "
        f"half_mapped={drops['half_mapped']} cross_ref={drops['cross_ref']}) "
        f"retention={100 * kept / pairs:.2f}%\n" if pairs else "no pairs\n")
    return TruthResult(pairs, kept, drops, refs, ids)


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sample_dir", type=Path,
                    help=f"sample directory; {TRUTH_NAME} and {UNCL_NAME} are written here")
    a = ap.parse_args()

    if not a.sample_dir.is_dir():
        die(1, f"read-truth: not a directory: {a.sample_dir}")
    run("-", a.sample_dir)
    return 0
