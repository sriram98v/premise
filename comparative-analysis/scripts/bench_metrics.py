"""Computing comparison metrics

  jaccard(A,B) = 1 - |A n B| / |A u B|
  ruzicka(x,y) = 1 - sum(min) / sum(max)
  cosine(x,y)  = <x,y> / (||x|| ||y||)
"""
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


def cosine_distance(x: dict, y: dict) -> float:
    keys = set(x) | set(y)
    dot = sum(x.get(k, 0.0) * y.get(k, 0.0) for k in keys)
    nx = sum(v * v for v in x.values()) ** 0.5
    ny = sum(v * v for v in y.values()) ** 0.5
    if nx == 0.0 or ny == 0.0:
        return float('nan')
    return 1.0 - dot / (nx * ny)


def truth_proportions(truth_counts: dict, n_reads: int):
    """Return (bwa, bwa_uc): the with-unclassified and classified-only truth vectors.

    truth_counts: {ref: count} over classified truth refs (version-stripped). RAW COUNTS, not
                  proportions — the unclassified mass is derived by subtraction below, so passing
                  a normalized profile makes that subtraction meaningless.
    n_reads:      total reads (classified + unclassified). Unclassified truth mass =
                  n_reads - sum(counts). bwa keeps the unclassified proportion; bwa.uc drops it
                  and renormalizes over classified refs.
    """
    classified = sum(truth_counts.values())
    uncl = n_reads - classified
    total = classified + uncl
    bwa = {}
    if total > 0:
        for ref, c in truth_counts.items():
            bwa[ref] = c / total
        bwa['unclassified'] = uncl / total   # 0 for synthetic/mixed (all reads have a source)
    bwa_uc = {}
    if classified > 0:
        for ref, c in truth_counts.items():
            bwa_uc[ref] = c / classified
    return bwa, bwa_uc


def build_method_profile(abund: dict, truth_refs: set, uncl_frac: float, has_uc: bool,
                         add_other: bool = True):
    """Build the comparable profile vectors for one method.
    """
    prof = {r: p * (1.0 - uncl_frac) for r, p in abund.items()}
    prof['unclassified'] = uncl_frac
    cl_sum = sum(v for r, v in prof.items() if r != 'unclassified')
    prof_uc = {r: (v / cl_sum if cl_sum else 0.0) for r, v in prof.items() if r != 'unclassified'}
    prof_uc['unclassified'] = 0.0
    if add_other:
        prof['other'] = sum(v for r, v in prof.items()
                            if r not in truth_refs and r not in ('unclassified', 'other'))
        prof_uc['other'] = sum(v for r, v in prof_uc.items()
                               if r not in truth_refs and r not in ('unclassified', 'other'))
    if not has_uc:
        prof_uc = prof
    return prof, prof_uc


def profile_distances(truth_counts: dict, n_reads: int, method_abund: dict | None,
                      uncl_frac: float = 0.0, has_uc: bool = False, add_other: bool = True) -> dict:
    """Ruzicka/Jaccard/Cosine distances of a method profile vs truth, both framings.

    method_abund: {ref: proportion} (version-stripped, sums to 1 over the method's refs), or None.
    uncl_frac:    method's unclassified fraction (n.uc/n.reads); 0 for pure profilers.
    has_uc:       method exposes a classified-only column (premise, centrifuger, assignment methods).
    add_other:    profile path (True) adds the 'other' aggregate row; assignment path (False) omits it.
    """
    if method_abund is None:
        nan = float('nan')
        out = {k: nan for k in ('cos', 'ruz', 'jac', 'cos.uc', 'ruz.uc', 'jac.uc')}
        out['fp.tx'] = out['fn.tx'] = None
        return out
    bwa, bwa_uc = truth_proportions(truth_counts, n_reads)
    truth_refs = set(truth_counts)
    prof, prof_uc = build_method_profile(method_abund, truth_refs, uncl_frac, has_uc, add_other)
    bwa_set = {k for k, v in bwa.items() if v > 0}
    bwa_uc_set = {k for k, v in bwa_uc.items() if v > 0}
    prof_set = {k for k, v in prof.items() if v > 0}
    prof_uc_set = {k for k, v in prof_uc.items() if v > 0}
    detected = prof_uc_set - SYNTHETIC_BINS
    return {
        'cos': cosine_distance(bwa, prof),
        'ruz': ruzicka_distance(bwa, prof),
        'jac': jaccard_distance(bwa_set, prof_set),
        'cos.uc': cosine_distance(bwa_uc, prof_uc),
        'ruz.uc': ruzicka_distance(bwa_uc, prof_uc),
        'jac.uc': jaccard_distance(bwa_uc_set, prof_uc_set),
        'fp.tx': len(detected - truth_refs),
        'fn.tx': len(truth_refs - detected),
    }


def _safe_div(a, b):
    return a / b if b else float('nan')


def precision_recall(truth_reads: dict, method_reads: dict | None, synthetic: bool) -> dict:
    """Per-read precision/recall counts, mirroring the R do.assignments block.

    truth_reads:  {readID: ref} — synthetic/mixed: every read (source acc from name);
                  real: only bwa-classified reads (ref is a true ref).
    method_reads: {readID: ref|None} where None == unclassified, or None if the method
                  has no per-read output (returns all-NaN).
    """
    if method_reads is None:
        nan = float('nan')
        return dict(tp=0, fp=0, tn=0, fn=0, tp_uc=0, fp_uc=0, tn_uc=0, fn_uc=0,
                    prec=nan, rec=nan, prec_uc=nan, rec_uc=nan)

    if synthetic:
        tp = fp = fn = 0
        for r, t in truth_reads.items():
            m = method_reads.get(r)
            if m is None:
                fn += 1
            elif m == t:
                tp += 1
            else:
                fp += 1
        return dict(tp=tp, fp=fp, tn=0, fn=fn, tp_uc=0, fp_uc=0, tn_uc=0, fn_uc=0,
                    prec=_safe_div(tp, tp + fp), rec=_safe_div(tp, tp + fn),
                    prec_uc=float('nan'), rec_uc=float('nan'))

    true_refs = set(truth_reads.values())
    reads = set(truth_reads) | set(method_reads)
    tp = fp = tn = fn = 0
    tp_uc = fp_uc = tn_uc = fn_uc = 0
    for r in reads:
        bwa_cls = r in truth_reads
        m = method_reads.get(r)
        m_true = (m is not None) and (m in true_refs)
        if bwa_cls and m_true:
            tp += 1
        if (not bwa_cls) and m_true:
            fp += 1
        if (not bwa_cls) and (not m_true):
            tn += 1
        if bwa_cls and (not m_true):
            fn += 1
        if (m is not None) and (m not in true_refs):
            fp_uc += 1
        if (not bwa_cls) and (m is None):
            tn_uc += 1
        if bwa_cls and (m is None):
            fn_uc += 1
    tp_uc = tp
    return dict(tp=tp, fp=fp, tn=tn, fn=fn, tp_uc=tp_uc, fp_uc=fp_uc, tn_uc=tn_uc, fn_uc=fn_uc,
                prec=_safe_div(tp, tp + fp), rec=_safe_div(tp, tp + fn),
                prec_uc=_safe_div(tp_uc, tp_uc + fp_uc), rec_uc=_safe_div(tp_uc, tp_uc + fn_uc))
