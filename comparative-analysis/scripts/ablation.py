#!/usr/bin/env python3
"""Ablation study for PREMISE across all four benchmark splits.

Sweeps each tunable parameter in isolation around the benchmark's production configuration and records, per setting.

Parameters swept (see `GRID`):
    mem    (-m)      minimum SMEM seed length -- seeding stage
    eps_2            minimum per-alignment match log-probability -- alignment stage
    eps_1            likelihood cutoff for dropping alignments before EM
    rho              L1 penalty weight -- EM stage
    omega            L1 penalty floor -- EM stage

Usage:
    python3 scripts/ablation.py                     # every split, every parameter
    python3 scripts/ablation.py --split syn-mix    # one split
    python3 scripts/ablation.py --split syn-mix --dataset Dataset-3
    python3 scripts/ablation.py --only mem,rho      # a subset of parameters
    python3 scripts/ablation.py --dry-run           # print the plan only
    python3 scripts/ablation.py --force             # ignore cached points
    python3 scripts/ablation.py --keep-large        # keep .aligns/.posteriors
"""
from __future__ import annotations

import argparse
import contextlib
import csv
import errno
import fcntl
import json
import os
import shutil
import subprocess
import sys
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import load_params
import utils
from bench_metrics import precision_recall, profile_distances
from utils import count_pairs, open_maybe_gz, parse_timemem, strain_of, strip_version


def _pre() -> dict[str, str]:
    return load_params.section("pre")


def _iters(split: str) -> int:
    """-i for a split: iters_real on real-iso, iters everywhere else."""
    pre = _pre()
    return int(pre["iters_real"] if split == "real-iso" else pre["iters"])

PREMISE = shutil.which("premise")


def index() -> Path:
    return utils.indexes() / "premise" / "sequences.fmidx"


def ablation_root() -> Path:
    return utils.results() / "ablation"


SPLITS = {
    "syn-iso":    dict(sub="synthetic/isolate", real=False,
                       datasets=["Dataset-1", "Dataset-2", "Dataset-3", "Dataset-4"],
                       reads=("reads_R1.fastq", "reads_R2.fastq")),
    "syn-mix":    dict(sub="synthetic/mixed", real=False,
                       datasets=["Dataset-1", "Dataset-2", "Dataset-3", "Dataset-4"],
                       reads=("reads_R1.fastq", "reads_R2.fastq")),
    "syn-mix-subtype": dict(sub="synthetic/mixed-subtype", real=False,
                       datasets=["Dataset-1", "Dataset-2", "Dataset-3", "Dataset-4"],
                       reads=("reads_R1.fastq", "reads_R2.fastq")),
    "real-iso":   dict(sub="real/isolate", real=True,
                       datasets=["SRR31013463", "SRR31013465", "SRR31013467", "SRR31013473"],
                       reads=("{ds}_1-filtered.ca.fastq", "{ds}_2-filtered.ca.fastq")),
    "real-mix":   dict(sub="real/mixed", real=True,
                       datasets=["SRR3360139", "SRR3360140", "SRR3360145", "SRR3360146"],
                       reads=("{ds}_1-filtered.ca.fastq", "{ds}_2-filtered.ca.fastq")),
}
SPLIT_ORDER = ["real-iso", "syn-iso", "syn-mix", "syn-mix-subtype", "real-mix"]  # cheapest first

def BASELINE() -> dict[str, str]:
    """The production point, read from params.toml on each call.

    A module-level dict here would read params.toml at import time, so importing this module
    would fail outright when the file is absent and would ignore a later --params.
    """
    pre = _pre()
    return {
        "mem": pre["mem"],
        "iter": pre["iters"],
        "eps_1": pre["eps_1"],
        "eps_2": pre["eps_2"],
        "em_threshold": pre["em_threshold"],
        "rho": pre["rho"],
        "omega": pre["omega"],
    }

GRID = {
    "mem":   [str(v) for v in range(6, 41)],
    "eps_2": ["0", "1e-36", "1e-27", "1e-18", "1e-12", "1e-9", "1e-6", "1e-3"],
    "eps_1": ["0", "1e-256", "1e-200", "1e-160", "1e-128", "1e-96", "1e-64",
              "1e-32", "1e-16", "1e-8", "1e-4"],
    "rho":   ["0", "10", "20", "30", "50", "150", "300", "600", "1200"],
    "omega": ["1e-30", "1e-20", "1e-10", "1e-8", "1e-6", "1e-4"],
}

PARAM_ORDER = ["mem", "eps_2", "eps_1", "rho", "omega"]

EXEC_ORDER = {"mem": lambda vs: sorted(vs, key=int, reverse=True)}

TIMEOUT_S = 7200

CSV_FIELDS = ["split", "dataset", "param", "value", "baseline", "wall_s", "rss_gb",
              "precision", "coverage", "ruzicka", "jaccard", "n_refs"]
METRIC_KEYS = ("wall_s", "rss_gb", "precision", "coverage", "ruzicka", "jaccard", "n_refs")


@dataclass(frozen=True)
class Target:
    """One (split, dataset) pair and every path derived from it."""
    split: str
    dataset: str

    @property
    def cfg(self) -> dict:
        return SPLITS[self.split]

    @property
    def ds_dir(self) -> Path:
        return utils.data_root() / "samples" / self.cfg["sub"] / self.dataset

    @property
    def out_root(self) -> Path:
        return ablation_root() / self.cfg["sub"] / self.dataset

    @property
    def reads(self) -> tuple[Path, Path]:
        r1, r2 = self.cfg["reads"]
        return (self.ds_dir / r1.format(ds=self.dataset),
                self.ds_dir / r2.format(ds=self.dataset))

    def baseline(self) -> dict:
        return dict(BASELINE(), iter=str(_iters(self.split)))

    def __str__(self) -> str:
        return f"{self.cfg['sub']}/{self.dataset}"


def split_csv(split: str) -> Path:
    """One CSV per split, covering that split's four datasets."""
    return ablation_root() / SPLITS[split]["sub"] / "ablation.csv"



def _first_existing(*cands: Path) -> Path | None:
    for c in cands:
        if c.exists():
            return c
    return None



@contextlib.contextmanager
def exclusive_lock():
    """Refuse to start if another sweep is already running, on any split or dataset.
    """
    ablation_root().mkdir(parents=True, exist_ok=True)
    lock = ablation_root() / ".sweep.lock"
    fd = os.open(lock, os.O_CREAT | os.O_RDWR, 0o644)
    try:
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as e:
            if e.errno in (errno.EACCES, errno.EAGAIN):
                holder = os.read(fd, 64).decode(errors="ignore").strip()
                sys.exit(f"another sweep is already running (pid {holder or '?'}); lock: {lock}")
            raise
        os.ftruncate(fd, 0)
        os.write(fd, f"{os.getpid()}\n".encode())
        yield
    finally:
        os.close(fd)



def load_context(t: Target) -> dict:
    """truth_reads, truth_counts, n_reads and n_input for one dataset.
    """
    r1, _ = t.reads
    if not r1.exists():
        raise SystemExit(f"{t}: reads not found at {r1}")

    if t.split in ("syn-iso", "syn-mix", "syn-mix-subtype"):
        truth = {}
        with open_maybe_gz(r1) as f:
            for i, line in enumerate(f):
                if i % 4:
                    continue
                rid = line[1:].rstrip("\n").rsplit("/", 1)[0]
                truth[rid] = strip_version(rid.split("_")[0])
        return dict(truth_reads=truth, truth_counts=dict(Counter(truth.values())),
                    n_reads=len(truth), n_input=len(truth))

    src = _first_existing(t.ds_dir / "truth_assignments.tsv",
                          t.ds_dir / "truth_assignments.tsv.gz")
    if src is None:
        raise SystemExit(f"{t}: {t.ds_dir}/truth_assignments.tsv[.gz] not found — build it with:"
                         f"\n  python3 scripts/prepare_real_samples.py --samples {t.dataset}")
    truth_seg = {}
    with open_maybe_gz(src) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0] and p[1]:
                truth_seg[p[0]] = strip_version(p[1])
    counts = Counter(truth_seg.values())
    if t.split != "real-mix":
        return dict(truth_reads=truth_seg, truth_counts=dict(counts), n_reads=len(truth_seg),
                    n_input=count_pairs(r1))

    truth = {q: strain_of(r) for q, r in truth_seg.items()}
    bad = set(truth.values()) - {"PR8", "WSN33"}
    if bad:
        raise SystemExit(f"{t}: {src}: {sorted(bad)} — every reference must be a PR8_/WSN33_ "
                         "segment on real-mix, or the strain truth collapses to 'other'")
    return dict(truth_reads=truth, truth_counts=dict(counts), n_reads=len(truth),
                n_input=count_pairs(r1))



def run_dir(t: Target, param: str, value: str) -> Path:
    """Directory for a grid point. Baseline settings all collapse to a single shared run."""
    return t.out_root / ("baseline" if value == BASELINE()[param] else f"{param}-{value}")


def build_cmd(t: Target, params: dict, out_base: Path) -> list[str]:
    r1, r2 = t.reads
    return [
        str(PREMISE), "query",
        "-s", str(index()),
        "-t", "0",
        "-1", str(r1),
        "-2", str(r2),
        "-o", str(out_base),
        "-m", params["mem"],
        "-i", params["iter"],
        "--eps_1", params["eps_1"],
        "--eps_2", params["eps_2"],
        "--em_threshold", params["em_threshold"],
        "--rho", params["rho"],
        "--omega", params["omega"],
    ]


def load_output(d: Path, base: str):
    """(assignments, abundance) from one run's .matches/.props, as analyze.load_premise does."""
    assign = None
    matches = d / f"{base}.matches"
    if matches.exists() and matches.stat().st_size:
        assign = {}
        with open(matches) as f:
            next(f)
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2:
                    assign[p[0]] = None if p[1] == "unclassified" else strip_version(p[1])
    abund = None
    props = d / f"{base}.props"
    if props.exists() and props.stat().st_size:
        c = Counter()
        with open(props) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2:
                    try:
                        c[strip_version(p[0])] += float(p[1])
                    except ValueError:
                        pass
        total = sum(c.values())
        abund = {k: v / total for k, v in c.items()} if total > 0 else None
    return assign, abund


def score(t: Target, d: Path, base: str, ctx: dict):
    """Metrics for one completed run, matching analyze.py's `pre` row for this split."""
    real = t.cfg["real"]
    assign, abund = load_output(d, base)
    wall, rss = parse_timemem(d / "time-mem")

    if t.split == "real-mix" and assign is not None:
        assign = {r: strain_of(v) for r, v in assign.items()}

    cov = None
    if assign is not None:
        cov = 100.0 * sum(1 for v in assign.values() if v is not None) / ctx["n_input"]

    dist = profile_distances(ctx["truth_counts"], ctx["n_reads"], abund,
                             uncl_frac=0.0, has_uc=True, add_other=True)
    pr = precision_recall(ctx["truth_reads"], assign, synthetic=not real)
    return {
        "wall_s": wall,
        "rss_gb": rss,
        "precision": pr["prec_uc"] if real else pr["prec"],
        "coverage": cov,
        "ruzicka": dist["ruz.uc"],
        "jaccard": dist["jac.uc"],
        "n_refs": len(abund) if abund else 0,
    }


def execute(t: Target, param: str, value: str, ctx: dict, keep_large: bool):
    """Run one grid point (or reuse an existing one) and return its metrics dict."""
    params = dict(t.baseline(), **{param: value})
    d = run_dir(t, param, value)
    base = "out"
    metrics_path = d / "metrics.json"
    if metrics_path.exists():
        return json.loads(metrics_path.read_text()), True

    d.mkdir(parents=True, exist_ok=True)
    cmd = ["/usr/bin/time", "-v", "-o", str(d / "time-mem")] + build_cmd(t, params, d / base)
    print(f"  $ premise query -m {params['mem']} --eps_1 {params['eps_1']} "
          f"--eps_2 {params['eps_2']} --rho {params['rho']} --omega {params['omega']}",
          flush=True)
    t0 = time.time()
    with open(d / "run.log", "w") as log:
        proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            rc = proc.wait(timeout=TIMEOUT_S)
        except subprocess.TimeoutExpired:
            with contextlib.suppress(ProcessLookupError):
                os.killpg(proc.pid, 15)
            try:
                proc.wait(timeout=30)
            except subprocess.TimeoutExpired:
                with contextlib.suppress(ProcessLookupError):
                    os.killpg(proc.pid, 9)
                proc.wait()
            print(f"    TIMEOUT after {TIMEOUT_S}s -- point abandoned", flush=True)
            return None, False
    if rc != 0:
        print(f"    FAILED rc={rc} after {time.time() - t0:.0f}s -- see {d / 'run.log'}",
              flush=True)
        return None, False

    m = score(t, d, base, ctx)
    m["params"] = params
    metrics_path.write_text(json.dumps(m, indent=2))
    if not keep_large:
        for ext in (".aligns", ".posteriors", ".matches"):
            f = d / f"{base}{ext}"
            if f.exists():
                f.unlink()
    print(f"    {m['wall_s']:.0f}s  {m['rss_gb']:.2f}GB  P={m['precision']:.4f} "
          f"cov={m['coverage']:.2f}%  ruz={m['ruzicka']:.4g}  jac={m['jaccard']:.4g}  "
          f"refs={m['n_refs']}", flush=True)
    return m, False


def write_csv(split: str) -> int:
    """Rebuild one split's CSV from every metrics.json on disk, in (dataset, GRID) order.
    """
    rows = []
    for ds in SPLITS[split]["datasets"]:
        t = Target(split, ds)
        for p in PARAM_ORDER:
            for v in GRID[p]:
                f = run_dir(t, p, v) / "metrics.json"
                if not f.exists():
                    continue
                m = json.loads(f.read_text())
                rows.append(dict(split=split, dataset=ds, param=p, value=v,
                                 baseline=int(v == BASELINE()[p]),
                                 **{k: m[k] for k in METRIC_KEYS}))
    out = split_csv(split)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        w.writeheader()
        w.writerows(rows)
    return len(rows)


def sweep(t: Target, params: list[str], keep_large: bool, failures: list):
    """Run every requested grid point for one dataset."""
    print(f"\n{'=' * 78}\n== {t} ==\n{'=' * 78}", flush=True)
    ctx = load_context(t)
    print(f"  truth: {ctx['n_reads']} reads, {len(ctx['truth_counts'])} references; "
          f"input: {ctx['n_input']} pairs", flush=True)
    for p in params:
        print(f"\n-- sweeping {p} (baseline {BASELINE()[p]}) --", flush=True)
        for v in EXEC_ORDER.get(p, list)(GRID[p]):
            m, cached = execute(t, p, v, ctx, keep_large)
            if m is None:
                failures.append((str(t), p, v))
            elif cached:
                print(f"  {p}={v}: cached ({m['wall_s']:.0f}s, P={m['precision']:.4f})",
                      flush=True)
            write_csv(t.split)



def run(splits_arg: list[str] | None = None, datasets: list[str] | None = None,
        only: list[str] | None = None, dry: bool = False, force: bool = False,
        keep_large: bool = False) -> int:
    """Sweep the parameter grid and write one CSV per split. -> 0, or raises ValueError."""

    splits = splits_arg or SPLIT_ORDER
    unknown = set(splits) - set(SPLITS)
    if unknown:
        raise ValueError(f"unknown split(s): {sorted(unknown)}")
    splits = [s for s in SPLIT_ORDER if s in splits]

    params = [p for p in PARAM_ORDER if p in (only or PARAM_ORDER)]
    unknown = set(only or []) - set(GRID)
    if unknown:
        raise ValueError(f"unknown parameter(s): {sorted(unknown)}")

    targets = []
    for s in splits:
        for ds in SPLITS[s]["datasets"]:
            if datasets and ds not in datasets:
                continue
            targets.append(Target(s, ds))
    if not targets:
        raise ValueError("no (split, dataset) pairs selected")

    n_points = sum(len(GRID[p]) for p in params)
    n_runs = len({run_dir(targets[0], p, v) for p in params for v in GRID[p]})
    print(f"ablation: {len(targets)} dataset(s) x {n_points} grid points over {len(params)} "
          f"parameter(s) = {n_runs * len(targets)} runs", flush=True)
    if dry:
        for t in targets:
            print(f"  {t}")
        for p in params:
            print(f"  {p}: {' '.join(GRID[p])}   (baseline {BASELINE()[p]})")
        return 0

    if PREMISE is None:
        raise ValueError("premise not found on PATH — enter the pinned toolchain first: nix develop ..#benchmark")
    if not index().exists():
        raise ValueError(f"index not found at {index()}")

    failures = []
    with exclusive_lock():
        if force:
            for t in targets:
                for p in params:
                    for v in GRID[p]:
                        (run_dir(t, p, v) / "metrics.json").unlink(missing_ok=True)
        for t in targets:
            sweep(t, params, keep_large, failures)

    print()
    for s in splits:
        print(f"wrote {write_csv(s):>4} rows -> {split_csv(s).relative_to(utils.data_root())}")
    if failures:
        print(f"\nFAILED points ({len(failures)}): {failures}")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--split", help=f"comma-separated subset of {', '.join(SPLIT_ORDER)}")
    ap.add_argument("--dataset", help="comma-separated subset of dataset names")
    ap.add_argument("--only", help="comma-separated subset of parameters to sweep")
    ap.add_argument("--dry-run", action="store_true", help="print the plan, run nothing")
    ap.add_argument("--force", action="store_true", help="re-run points that already have metrics")
    ap.add_argument("--keep-large", action="store_true",
                    help="keep .aligns/.posteriors/.matches (~120 MB/run)")
    args = ap.parse_args()
    try:
        return run(args.split.split(",") if args.split else None,
                   args.dataset.split(",") if args.dataset else None,
                   args.only.split(",") if args.only else None,
                   args.dry_run, args.force, args.keep_large)
    except ValueError as e:
        print(f"ablation: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
