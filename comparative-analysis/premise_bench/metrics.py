"""Computing comparison metrics"""
from __future__ import annotations

SYNTHETIC_BINS = frozenset({'other', 'unclassified'})


def jaccard_distance(set_a, set_b) -> float:
    inter = len(set_a & set_b)
    union = len(set_a) + len(set_b) - inter
    if union == 0:
        return float('nan')
    return 1.0 - inter / union


def ruzicka_distance(x: dict, y: dict) -> float:
    keys = set(x) | set(y)
    smin = smax = 0.0
    for k in keys:
        a = x.get(k, 0.0)
        b = y.get(k, 0.0)
        smin += min(a, b)
        smax += max(a, b)
    if smax == 0.0:
        return float('nan')
    return 1.0 - smin / smax


def truth_profile(truth_counts: dict) -> dict:
    """Each true source's share of the truth-labelled reads."""
    total = sum(truth_counts.values())
    return {ref: c / total for ref, c in truth_counts.items()} if total > 0 else {}


def method_profile(abund: dict, truth_refs: set, add_other: bool = True) -> dict:
    """A method's abundance rescaled to sum to 1, with every reference outside the truth
    optionally pooled into an extra 'other' bin."""
    total = sum(abund.values())
    prof = {r: (v / total if total else 0.0) for r, v in abund.items()}
    if add_other:
        prof['other'] = sum(v for r, v in prof.items() if r not in truth_refs and r != 'other')
    return prof


def profile_distances(truth_counts: dict, method_abund: dict | None,
                      add_other: bool = True) -> dict:
    """Ruzicka and Jaccard distances of a method's profile from the truth, and the
    reference-level false positives (fp.tx) and false negatives (fn.tx)."""
    if method_abund is None:
        nan = float('nan')
        return {'ruz': nan, 'jac': nan, 'fp.tx': None, 'fn.tx': None}
    truth = truth_profile(truth_counts)
    truth_refs = set(truth_counts)
    prof = method_profile(method_abund, truth_refs, add_other)
    truth_set = {k for k, v in truth.items() if v > 0}
    prof_set = {k for k, v in prof.items() if v > 0}
    detected = prof_set - SYNTHETIC_BINS
    return {
        'ruz': ruzicka_distance(truth, prof),
        'jac': jaccard_distance(truth_set, prof_set),
        'fp.tx': len(detected - truth_refs),
        'fn.tx': len(truth_refs - detected),
    }


def precision(truth_reads: dict, method_reads: dict | None, synthetic: bool) -> float:
    """Per-read precision; NaN when the method has no per-read output.

    synthetic: every read has a source, so a classified read is a true positive when it
               names exactly that source.
    real:      a classified read is a true positive when the read is truth-labelled and it
               names one of the sample's true sources, and a false positive when it names any
               reference outside them.
    """
    if method_reads is None:
        return float('nan')
    tp = fp = 0
    if synthetic:
        for r, t in truth_reads.items():
            m = method_reads.get(r)
            if m is None:
                continue
            if m == t:
                tp += 1
            else:
                fp += 1
    else:
        true_refs = set(truth_reads.values())
        for r, m in method_reads.items():
            if m is None:
                continue
            if m not in true_refs:
                fp += 1
            elif r in truth_reads:
                tp += 1
    return tp / (tp + fp) if tp + fp else float('nan')
