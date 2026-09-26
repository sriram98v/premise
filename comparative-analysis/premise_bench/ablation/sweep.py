"""Ablation study for PREMISE across all four benchmark splits"""
from __future__ import annotations

import argparse
import contextlib
import csv
import errno
import fcntl
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from ..utils import count_pairs, nonempty, parse_timemem, real_truth, strain_of, synthetic_truth
from ..config import data_root, indexes, results, section
from ..methods import premise
from ..metrics import precision, profile_distances
from ..splits import SPLITS, Split



def _pre() -> dict[str, str]:
    return section("pre")


def _iters(split: str) -> int:
    """-i for a split: iters_real on real-iso, iters everywhere else."""
    pre = _pre()
    return int(pre["iters_real"] if split == "real-iso" else pre["iters"])

PREMISE = shutil.which("premise")


def index() -> Path:
    return indexes() / "premise" / "sequences.fmidx"


def reference() -> Path:
    """The FASTA the benchmark builds the production index from (runner.CLEANED)."""
    return indexes() / "sequences-cleaned.fasta"


def index_cache() -> Path:
    """Per-rate indexes, built once and shared by every (split, dataset) point."""
    return ablation_root() / "indexes"


def ablation_root() -> Path:
    return results() / "ablation"


SWEEP_ORDER = ["real-iso", "syn-iso", "syn-mix-sub", "syn-mix-strain", "real-mix"]  # cheapest first

def BASELINE() -> dict[str, str]:
    """The production point, read from params.toml on each call"""
    pre = _pre()
    return {
        "mem": pre["mem"],
        "iter": pre["iters"],
        "eps_1": pre["eps_1"],
        "eps_2": pre["eps_2"],
        "em_threshold": pre["em_threshold"],
        "rho": pre["rho"],
        "omega": pre["omega"],
        "sa_sample_rate": pre["sa_sample_rate"],
    }

GRID = {
    "mem":   ["20", "21", "22", "23", "24"],
    "eps_2": ["1e-36", "1e-27", "1e-18", "1e-12", "1e-9"],
    "eps_1": ["0", "1e-256", "1e-200", "1e-160", "1e-128"],
    "rho":   ["30", "50", "150", "300", "600"],
    "omega": ["1e-30", "1e-20", "1e-10", "1e-8", "1e-6"],
    "sa_sample_rate": ["1", "2", "4", "8", "16", "32", "64"],
}

PARAM_ORDER = ["mem", "eps_2", "eps_1", "rho", "omega", "sa_sample_rate"]

BUILD_PARAMS = {"sa_sample_rate"}

EXEC_ORDER = {"mem": lambda vs: sorted(vs, key=int, reverse=True)}

TIMEOUT_S = 7200

CSV_FIELDS = ["split", "dataset", "param", "value", "baseline", "wall_s", "rss_gb",
              "precision", "coverage", "ruzicka", "jaccard", "n_refs",
              "index_bytes", "build_s"]
METRIC_KEYS = ("wall_s", "rss_gb", "precision", "coverage", "ruzicka", "jaccard", "n_refs",
               "index_bytes", "build_s")


@dataclass(frozen=True)
class Target:
    """One (split, dataset) pair and every path derived from it."""
    split: str
    dataset: str

    @property
    def cfg(self) -> Split:
        return SPLITS[self.split]

    @property
    def ds_dir(self) -> Path:
        return data_root() / "samples" / self.cfg.sub / self.dataset

    @property
    def out_root(self) -> Path:
        return ablation_root() / self.cfg.sub / self.dataset

    @property
    def reads(self) -> tuple[Path, Path]:
        return (self.ds_dir / self.cfg.read_name(self.dataset, 1),
                self.ds_dir / self.cfg.read_name(self.dataset, 2))

    def baseline(self) -> dict:
        return dict(BASELINE(), iter=str(_iters(self.split)))

    def __str__(self) -> str:
        return f"{self.cfg.sub}/{self.dataset}"


def split_csv(split: str) -> Path:
    """One CSV per split, covering that split's four datasets."""
    return ablation_root() / SPLITS[split].sub / "ablation.csv"



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

    if not t.cfg.real:
        truth = synthetic_truth(r1)
        return dict(truth_reads=truth, truth_counts=dict(Counter(truth.values())),
                    n_reads=len(truth), n_input=len(truth))

    tsv = t.ds_dir / "truth_assignments.tsv"
    try:
        truth, counts, n = real_truth(tsv, t.split)
    except FileNotFoundError:
        raise SystemExit(f"{t}: {t.ds_dir}/truth_assignments.tsv[.gz] not found — build it with:"
                         f"\n  python3 -m premise_bench prepare-real --samples {t.dataset}") from None
    except ValueError as e:
        raise SystemExit(f"{t}: {e}") from None
    return dict(truth_reads=truth, truth_counts=counts, n_reads=n, n_input=count_pairs(r1))


def run_dir(t: Target, param: str, value: str) -> Path:
    """Directory for a grid point. Baseline settings all collapse to a single shared run"""
    if param in BUILD_PARAMS:
        return t.out_root / f"{param}-{value}"
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


def align_cmd(t: Target, params: dict, idx: Path, out_base: Path) -> list[str]:
    """The alignment stage alone, against an index built at a specific sampling rate"""
    r1, r2 = t.reads
    return [
        str(PREMISE), "align",
        "-s", str(idx),
        "-t", "0",
        "-1", str(r1),
        "-2", str(r2),
        "-o", f"{out_base}.aligns",
        "-m", params["mem"],
        "--eps_2", params["eps_2"],
    ]


def load_output(d: Path, base: str):
    """(assignments, abundance) from one run's .matches/.props; abundance None when empty."""
    assign, abund = premise.load(d, base)
    return assign, abund or None


def score(t: Target, d: Path, base: str, ctx: dict):
    """Metrics for one completed run, matching the `pre` values `premise_bench analyze` writes for this split."""
    real = t.cfg.real
    assign, abund = load_output(d, base)
    wall, rss = parse_timemem(d / "time-mem")

    if t.split == "real-mix" and assign is not None:
        assign = {r: strain_of(v) for r, v in assign.items()}

    cov = None
    if assign is not None:
        cov = 100.0 * sum(1 for v in assign.values() if v is not None) / ctx["n_input"]

    dist = profile_distances(ctx["truth_counts"], abund, add_other=True)
    return {
        "wall_s": wall,
        "rss_gb": rss,
        "precision": precision(ctx["truth_reads"], assign, synthetic=not real),
        "coverage": cov,
        "ruzicka": dist["ruz"],
        "jaccard": dist["jac"],
        "n_refs": len(abund) if abund else 0,
    }


def spawn(cmd: list[str], log_path: Path) -> int | None:
    """Run `cmd` to completion under TIMEOUT_S. -> exit code, or None if it timed out."""
    with open(log_path, "w") as log:
        proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            return proc.wait(timeout=TIMEOUT_S)
        except subprocess.TimeoutExpired:
            with contextlib.suppress(ProcessLookupError):
                os.killpg(proc.pid, 15)
            try:
                proc.wait(timeout=30)
            except subprocess.TimeoutExpired:
                with contextlib.suppress(ProcessLookupError):
                    os.killpg(proc.pid, 9)
                proc.wait()
            return None


def index_for(rate: str) -> Path | None:
    """The index built at `rate`, building it on first use. -> path, or None if the build failed"""
    cache = index_cache()
    cache.mkdir(parents=True, exist_ok=True)
    idx = cache / f"sequences-sa{rate}.fmidx"
    if nonempty(idx):
        return idx
    if not reference().exists():
        raise ValueError(f"reference FASTA not found at {reference()}")
    print(f"  $ premise build --sa_sample_rate {rate}", flush=True)
    cmd = ["/usr/bin/time", "-v", "-o", str(cache / f"build-sa{rate}.time"),
           str(PREMISE), "build", "-s", str(reference()), "-o", str(idx),
           "--sa_sample_rate", rate]
    rc = spawn(cmd, cache / f"build-sa{rate}.log")
    if rc != 0:
        idx.unlink(missing_ok=True)  # never leave a truncated index to be reused as cached
        print(f"    BUILD FAILED rc={rc} -- see {cache / f'build-sa{rate}.log'}", flush=True)
        return None
    return idx


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def score_build(t: Target, d: Path, base: str, rate: str, idx: Path) -> dict:
    """Resources for one build-param point, plus the invariant that guards the whole sweep"""
    wall, rss = parse_timemem(d / "time-mem")
    build_wall, build_rss = parse_timemem(index_cache() / f"build-sa{rate}.time")
    aligns = d / f"{base}.aligns"
    if not nonempty(aligns):
        print(f"    WARNING: no alignments at {aligns} -- identity gate has nothing to hash",
              flush=True)
    return {
        "wall_s": wall,
        "rss_gb": rss,
        "precision": None,
        "coverage": None,
        "ruzicka": None,
        "jaccard": None,
        "n_refs": None,
        "index_bytes": idx.stat().st_size,
        "build_s": build_wall,
        "build_rss_gb": build_rss,
        "aligns_sha256": sha256(aligns) if nonempty(aligns) else None,
    }


def execute(t: Target, param: str, value: str, ctx: dict, keep_large: bool):
    """Run one grid point (or reuse an existing one) and return its metrics dict."""
    params = dict(t.baseline(), **{param: value})
    d = run_dir(t, param, value)
    base = "out"
    metrics_path = d / "metrics.json"
    if metrics_path.exists():
        return json.loads(metrics_path.read_text()), True

    build_param = param in BUILD_PARAMS
    idx = None
    if build_param:
        idx = index_for(value)
        if idx is None:
            return None, False

    d.mkdir(parents=True, exist_ok=True)
    if build_param:
        cmd = (["/usr/bin/time", "-v", "-o", str(d / "time-mem")]
               + align_cmd(t, params, idx, d / base))
        print(f"  $ premise align -s sequences-sa{value}.fmidx -m {params['mem']} "
              f"--eps_2 {params['eps_2']}", flush=True)
    else:
        cmd = ["/usr/bin/time", "-v", "-o", str(d / "time-mem")] + build_cmd(t, params, d / base)
        print(f"  $ premise query -m {params['mem']} --eps_1 {params['eps_1']} "
              f"--eps_2 {params['eps_2']} --rho {params['rho']} --omega {params['omega']}",
              flush=True)
    t0 = time.time()
    rc = spawn(cmd, d / "run.log")
    if rc is None:
        print(f"    TIMEOUT after {TIMEOUT_S}s -- point abandoned", flush=True)
        return None, False
    if rc != 0:
        print(f"    FAILED rc={rc} after {time.time() - t0:.0f}s -- see {d / 'run.log'}",
              flush=True)
        return None, False

    m = score_build(t, d, base, value, idx) if build_param else score(t, d, base, ctx)
    m["params"] = params
    metrics_path.write_text(json.dumps(m, indent=2))
    if not keep_large:
        for ext in (".aligns", ".posteriors", ".matches"):
            f = d / f"{base}{ext}"
            if f.exists():
                f.unlink()
    if build_param:
        print(f"    {m['wall_s']:.0f}s  {m['rss_gb']:.2f}GB  "
              f"index={m['index_bytes'] / 1e6:.1f}MB  build={m['build_s']:.0f}s", flush=True)
    else:
        print(f"    {m['wall_s']:.0f}s  {m['rss_gb']:.2f}GB  P={m['precision']:.4f} "
              f"cov={m['coverage']:.2f}%  ruz={m['ruzicka']:.4g}  jac={m['jaccard']:.4g}  "
              f"refs={m['n_refs']}", flush=True)
    return m, False


def write_csv(split: str) -> int:
    """Rebuild one split's CSV from every metrics.json on disk, in (dataset, GRID) order.
    """
    rows = []
    for ds in SPLITS[split].datasets:
        t = Target(split, ds)
        for p in PARAM_ORDER:
            for v in GRID[p]:
                f = run_dir(t, p, v) / "metrics.json"
                if not f.exists():
                    continue
                m = json.loads(f.read_text())
                rows.append(dict(split=split, dataset=ds, param=p, value=v,
                                 baseline=int(v == BASELINE()[p]),
                                 **{k: m.get(k) for k in METRIC_KEYS}))
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
    ctx = None
    if any(p not in BUILD_PARAMS for p in params):
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
                extra = (f"index={m['index_bytes'] / 1e6:.1f}MB" if p in BUILD_PARAMS
                         else f"P={m['precision']:.4f}")
                print(f"  {p}={v}: cached ({m['wall_s']:.0f}s, {extra})", flush=True)
            write_csv(t.split)



def check_build_invariants(plan: list, build_params: list[str]) -> None:
    """Report whether every value of a build param produced identical alignments.
    """
    for p in build_params:
        for t, ps in plan:
            if p not in ps:
                continue
            seen: dict[str, list[str]] = {}
            for v in GRID[p]:
                f = run_dir(t, p, v) / "metrics.json"
                if not f.exists():
                    continue
                h = json.loads(f.read_text()).get("aligns_sha256")
                seen.setdefault(h or "MISSING", []).append(v)
            if len(seen) <= 1 and "MISSING" not in seen:
                print(f"{t}: {p} alignments identical across {len(next(iter(seen.values())))} "
                      f"value(s)")
            else:
                groups = "; ".join(f"{k[:12]}: {' '.join(vs)}" for k, vs in seen.items())
                print(f"{t}: !! {p} alignments DIFFER -- {groups}")


def run(splits_arg: list[str] | None = None, datasets: list[str] | None = None,
        only: list[str] | None = None, dry: bool = False, force: bool = False,
        keep_large: bool = False) -> int:
    """Sweep the parameter grid and write one CSV per split. -> 0, or raises ValueError."""

    splits = splits_arg or SWEEP_ORDER
    unknown = set(splits) - set(SPLITS)
    if unknown:
        raise ValueError(f"unknown split(s): {sorted(unknown)}")
    splits = [s for s in SWEEP_ORDER if s in splits]

    params = [p for p in PARAM_ORDER if p in (only or PARAM_ORDER)]
    unknown = set(only or []) - set(GRID)
    if unknown:
        raise ValueError(f"unknown parameter(s): {sorted(unknown)}")

    targets = []
    for s in splits:
        for ds in SPLITS[s].datasets:
            if datasets and ds not in datasets:
                continue
            targets.append(Target(s, ds))
    if not targets:
        raise ValueError("no (split, dataset) pairs selected")

    plan = [(t, params) for t in targets]
    n_runs = sum(len({run_dir(t, p, v) for p in ps for v in GRID[p]}) for t, ps in plan)
    print(f"ablation: {len(plan)} dataset(s) over {len(params)} parameter(s) = "
          f"{n_runs} runs", flush=True)
    if dry:
        for t, ps in plan:
            print(f"  {t}: {' '.join(ps)}")
        for p in params:
            kind = " [build-time]" if p in BUILD_PARAMS else ""
            print(f"  {p}: {' '.join(GRID[p])}   (baseline {BASELINE()[p]}){kind}")
        return 0

    if PREMISE is None:
        raise ValueError("premise not found on PATH — enter the pinned toolchain first: nix develop ..#benchmark")
    if any(p not in BUILD_PARAMS for p in params) and not index().exists():
        raise ValueError(f"index not found at {index()}")

    failures = []
    with exclusive_lock():
        if force:
            for t, ps in plan:
                for p in ps:
                    for v in GRID[p]:
                        (run_dir(t, p, v) / "metrics.json").unlink(missing_ok=True)
                        if p in BUILD_PARAMS:
                            (index_cache() / f"sequences-sa{v}.fmidx").unlink(missing_ok=True)
        for t, ps in plan:
            sweep(t, ps, keep_large, failures)

    print()
    check_build_invariants(plan, [p for p in params if p in BUILD_PARAMS])
    for s in splits:
        print(f"wrote {write_csv(s):>4} rows -> {split_csv(s).relative_to(data_root())}")
    if failures:
        print(f"\nFAILED points ({len(failures)}): {failures}")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--split", help=f"comma-separated subset of {', '.join(SWEEP_ORDER)}")
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
