"""Ground truth per dataset"""
from __future__ import annotations

from collections import Counter

from ..utils import count_pairs, real_truth, synthetic_truth
from ..config import data_root
from ..splits import SPLITS, split_read_path

def load_truth(split: str, ds: str):
    """Return (truth_reads, truth_counts, n_reads)"""
    info = SPLITS[split]
    if not info.real:
        truth = synthetic_truth(split_read_path(split, ds, mate=1))
        return truth, dict(Counter(truth.values())), len(truth)
    tsv = data_root() / "samples" / info.sub / ds / "truth_assignments.tsv"
    try:
        return real_truth(tsv, split)
    except FileNotFoundError:
        raise ValueError(
            f"{tsv} not found — build it with:\n"
            f"  python3 -m premise_bench prepare-real --samples {ds}") from None


def input_read_count(split: str, ds: str):
    """Total read (pair) count in the dataset's R1 input fastq, or None if unavailable.
    """
    try:
        p = split_read_path(split, ds, mate=1)
        if not p.exists() or p.stat().st_size == 0:
            return None
        return count_pairs(p)
    except (OSError, EOFError, KeyError):
        # Coverage degrades to a dash rather than taking the whole table down with it.
        return None
