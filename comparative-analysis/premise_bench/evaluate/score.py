"""Scoring one dataset"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

from ..utils import parse_timemem, strain_of
from ..config import results
from ..methods import METHOD_BY_CODE, METHODS
from ..metrics import precision, profile_distances
from ..splits import SPLITS
from .truth import input_read_count, load_truth

def _suffix(d: Path, code: str, split: str, ds: str) -> str:
    """'' for the plain <dataset> base the drivers now write, '.ca' for older trees.
    """
    if split != "real-iso" or code not in ("pre", "kmp"):
        return ""
    if any(not p.name.startswith(f"{ds}.ca.") for p in d.glob(f"{ds}.*")):
        return ""
    return ".ca" if any(d.glob(f"{ds}.ca.*")) else ""


def load_method(code: str, split: str, ds: str):
    info = SPLITS[split]
    d = results() / METHOD_BY_CODE[code].tool / info.sub / ds
    base = ds + _suffix(d, code, split, ds)
    return METHOD_BY_CODE[code].load(d, base)


def _check_read_counts(code: str, split: str, ds: str, assign: dict, denom: int | None) -> None:
    """Warn when a method's per-read output and the coverage denominator disagree.
    """
    if denom is None:
        return
    n = len(assign)
    if n > denom:
        print(f"WARNING: {code} {split}/{ds}: {n} reads in output exceeds denominator {denom} "
              f"-- results and input FASTQ disagree; coverage is meaningless", file=sys.stderr)
    elif METHOD_BY_CODE[code].complete_output and n != denom:
        print(f"WARNING: {code} {split}/{ds}: {n} reads in output vs denominator {denom} "
              f"(delta {n - denom:+d}); {code} emits one row per read, so these should match",
              file=sys.stderr)


def compute_dataset(split: str, ds: str):
    info = SPLITS[split]
    truth_reads, truth_counts, _ = load_truth(split, ds)
    rows = {}
    for m in METHODS:
        assign, abund = load_method(m.code, split, ds)
        if split == "real-mix" and assign is not None:
            assign = {r: strain_of(v) for r, v in assign.items()}
        cr = None
        if assign is not None:
            classified = sum(1 for v in assign.values() if v is not None)
            denom = input_read_count(split, ds)
            _check_read_counts(m.code, split, ds, assign, denom)
            if denom:
                cr = 100.0 * classified / denom
        dist = profile_distances(truth_counts, abund, add_other=m.add_other)
        prec = precision(truth_reads, assign if m.perread else None, synthetic=not info.real)
        rows[m.code] = dict(dist=dist, prec=prec, cr=cr,
                            has_data=(abund is not None or assign is not None))
    return rows


def abund_cell(row):
    return row["dist"]["ruz"], row["dist"]["jac"]


def fpfn_cell(row):
    """(false positives, false negatives) at the reference level, or (None, None) to dash.
    """
    dist = row["dist"]
    return (dist.get("fp.tx"), dist.get("fn.tx"))


def pr_cell(row):
    """(precision, coverage)
    """
    return (row["prec"], row["cr"])


def _isnum(x):
    return isinstance(x, (int, float)) and not math.isnan(x)


def _timed_out(tm) -> bool:
    """True when the time-mem file records a `timeout(1)` kill (exit status 124)"""
    tm = Path(tm)
    if not tm.exists():
        return False
    for line in tm.read_text(errors="ignore").splitlines():
        if "Exit status" in line:
            return line.split(":")[-1].strip() == "124"
    return False


def _timemem_path(split, m, ds):
    d = results() / m.tool / SPLITS[split].sub / ds
    tm = d / "time-mem"
    if not tm.exists() and (d / "time-mem.ca").exists():
        tm = d / "time-mem.ca"
    return tm


def _timeout_map(split, datasets):
    """{method code: set of datasets whose classify run was killed at the time limit}."""
    tmap = {}
    for m in METHODS:
        for ds in datasets:
            if _timed_out(_timemem_path(split, m, ds)):
                tmap.setdefault(m.code, set()).add(ds)
    return tmap


CSV_FIELDS = ("method", "method_name", "dataset", "metric", "value", "status")


def split_csv_path(split: str) -> Path:
    """The per-split CSV of every scored value, written by write_split_csv."""
    return results() / "tables" / f"comparative-{split}.csv"


def gather_records(split: str, datasets=None, per_ds=None, tmap=None) -> list[dict]:
    """Tidy (method, method_name, dataset, metric, value, status) records for one split"""
    info = SPLITS[split]
    dss = [d for d in info.datasets if not datasets or d in datasets]
    if per_ds is None:
        per_ds = {ds: compute_dataset(split, ds) for ds in dss}
    dss = [d for d in dss if d in per_ds]
    if tmap is None:
        tmap = _timeout_map(split, dss)
    recs = []

    def rec(m, ds, metric, value, status=None):
        missing = value is None or not _isnum(value)
        recs.append(dict(method=m.code, method_name=m.name, dataset=ds,
                         metric=metric, value=value if not missing else None,
                         status=status or ("missing" if missing else "ok")))

    for ds in dss:
        rows = per_ds[ds]
        for m in METHODS:
            row = rows[m.code]
            if m.prec_rec:
                prec, cov = pr_cell(row)
                rec(m, ds, "precision", prec)
                rec(m, ds, "coverage", cov)
            fp, fn = fpfn_cell(row)
            rec(m, ds, "fp", fp)
            rec(m, ds, "fn", fn)
            ruz, jac = abund_cell(row)
            rec(m, ds, "ruzicka", ruz)
            rec(m, ds, "jaccard", jac)
            if ds in tmap.get(m.code, set()):
                for key in ("wall_s", "rss_gb"):
                    rec(m, ds, key, None, status="timeout")
            else:
                wall, rss = parse_timemem(_timemem_path(split, m, ds))
                rec(m, ds, "wall_s", wall)
                rec(m, ds, "rss_gb", rss)
    return recs


def write_split_csv(split: str, recs: list[dict]) -> Path:
    path = split_csv_path(split)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        w.writeheader()
        w.writerows({**r, "value": "" if r["value"] is None else r["value"]}
                    for r in recs)
    return path
