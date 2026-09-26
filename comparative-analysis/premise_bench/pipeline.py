"""The end-to-end benchmark driver"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
import time
from pathlib import Path

from . import config
from .utils import die, nonempty
from .data import decoys as decoy_refs
from .data import prepare_real as prepare_real_samples
from .evaluate import analyze
from .methods import METHOD_BY_CODE, METHODS, Method
from .runner import Ctx, Job, require_tools, resolve
from .splits import SPLITS as SPLIT_BY_CODE
from .splits import Split

DB_BUILD_CSV = "db_build.csv"


DB_BUILD_HEADER = "method,size_before_bytes,size_after_bytes,build_seconds,exit"


BASE_TOOLS = ("bwa-mem2", "samtools", "cutadapt", "python3", "iss", "fasterq-dump")


def idx_size(ctx: Ctx, m: Method) -> int:
    """On-disk bytes of a method's index"""
    total = 0
    for rel in m.index_paths:
        p = ctx.root / rel
        if p.is_dir():
            total += sum(f.stat().st_size for f in p.rglob("*") if f.is_file())
        elif p.exists():
            total += p.stat().st_size
    return total


def build_index(ctx: Ctx, m: Method) -> int:
    """Build one index, appending its size/time/exit row to db_build.csv. Returns the exit code."""
    before = idx_size(ctx, m)
    t0 = time.monotonic()
    rcs = m.build(ctx)
    secs = time.monotonic() - t0
    after = idx_size(ctx, m)
    rc = next((c for c in rcs if c != 0), 0)
    if not ctx.dry_run:
        with open(ctx.root / DB_BUILD_CSV, "a") as f:
            f.write(f"{m.tool},{before},{after},{secs:.9f},{rc}\n")
    print(f"   built {m.tool} in {secs:.9f}s (index {after / 1e9:.3f} GB, exit {rc})")
    return rc


def prepare_real(ctx: Ctx, splits: list[Split]) -> None:
    """Prepare every selected real split whose filtered reads or truth are missing.
    """
    for code in ("real-iso", "real-mix"):
        sp = next((s for s in splits if s.code == code), None)
        if sp is None:
            continue
        sdir = ctx.root / sp.sample_dir
        if not sdir.is_dir():
            continue
        need = False
        for d in sorted(sdir.glob("SRR*")):
            if not d.is_dir():
                continue
            r1 = d / sp.read_name(d.name, 1)
            if not (nonempty(r1) and nonempty(d / "truth_assignments.tsv")):
                need = True
        if not need:
            continue
        print(f" preparing {code} reads + truth (premise_bench prepare-real)…")
        if ctx.dry_run:
            print(f"   $ prepare_real_samples.run({code})")
            continue
        try:
            prepare_real_samples.run(code, int(ctx.threads), decoys=False)
        except BaseException as e:  # noqa: BLE001 — prepare die()s with SystemExit on a bad tree
            print(f" WARN: prepare_real_samples failed ({e}); {code} may be incomplete")


def add_decoys(ctx: Ctx) -> None:
    """Bring the decoy references in indexes/ in line with [decoy]"""
    print(" decoy references ([decoy] in params.toml)…")
    changed = decoy_refs.update(dry=ctx.dry_run, report=False)
    if changed and ctx.skip_build and not ctx.dry_run:
        die(1, "FATAL: the decoys in indexes/sequences-cleaned.fasta changed, so every index "
            "built from it is stale.",
            "       Rebuild the indexes (drop --skip-build).")


def make_job(ctx: Ctx, m: Method, sp: Split, ds: Path) -> Job | None:
    base = ds.name
    if m.real_only and sp.synthetic and not ctx.kmcp_allow_synthetic:
        print(f"    skip {m.tool} ({sp.code}/{base}): synthetic splits are "
              f"{m.tool}-free by design")
        return None
    outdir = f"results/{m.tool}/{sp.sub}/{base}"
    if not ctx.dry_run:
        (ctx.root / outdir).mkdir(parents=True, exist_ok=True)
    r1, r2 = (f"{sp.sample_dir}/{base}/" + sp.read_name(base, n) for n in (1, 2))
    if not (nonempty(ctx.root / r1) and nonempty(ctx.root / r2)):
        print(f"    skip {m.tool} ({sp.code}/{base}): missing reads")
        return None
    return Job(sp, base, outdir, r1, r2)


def _analyze(ctx: Ctx, splits: list[str]) -> None:
    """Score and write the per-split CSVs"""
    if ctx.dry_run:
        print(f"   $ analyze.run(--splits {','.join(splits)})")
        return
    try:
        print(analyze.run(splits))
    except Exception as e:                                  # noqa: BLE001 — see docstring
        print(f"analyze: {e}", file=sys.stderr)


def classify_split(ctx: Ctx, methods: list[Method], sp: Split) -> None:
    sdir = ctx.root / sp.sample_dir
    if not sdir.is_dir():
        print(f" (skip {sp.code}: {sp.sample_dir} not present)")
        return
    count = 0
    for ds in sorted(p for p in sdir.iterdir() if p.is_dir()):
        print(f" -- {sp.code} / {ds.name} --")
        for m in methods:
            job = make_job(ctx, m, sp, ds)
            if job is None:
                continue
            rc = m.classify(ctx, job)
            if rc == 124:
                print(f"    TIMEOUT ({ctx.timeout}s) {m.tool} ({sp.code}/{ds.name}) — "
                      f"recorded as failed, continuing")
        count += 1
        if ctx.max_ds > 0 and count >= ctx.max_ds:
            print(f" (stopping {sp.code} at {ctx.max_ds} dataset(s))")
            break


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--params", type=Path, default=None,
                    help="params.toml to use (default: the one beside the premise_bench package)")
    ap.add_argument("--bench-data", default=None, help="data root; overrides $BENCH_DATA")
    ap.add_argument("--threads", default=None, help="worker threads for every method")
    ap.add_argument("--methods", default=None, help='whitespace-separated codes, e.g. "pre cen syl"')
    ap.add_argument("--splits", default=None, help="whitespace-separated split codes")
    ap.add_argument("--max-ds", default=None, help="datasets per split; 0 = all")
    ap.add_argument("--skip-build", action="store_const", const="1", default=None,
                    help="reuse existing indexes/ and db_build.csv")
    ap.add_argument("--build-only", action="store_true",
                    help="stop after stage 1 (index construction); no classification or scoring")
    ap.add_argument("--kmcp-threads", default=None)
    ap.add_argument("--kmcp-allow-synthetic", action="store_const", const="1", default=None,
                    help="run kmcp on the synthetic splits too (diagnostic reruns)")
    ap.add_argument("--keep-going", action="store_true",
                    help="carry on after a failed index build instead of aborting")
    ap.add_argument("--dry-run", action="store_true", help="print every command, run nothing")
    return ap.parse_args(argv)


def build_ctx(a: argparse.Namespace) -> tuple[Ctx, list[Method], list[Split]]:
    params_path = a.params if a.params is not None else config.find_file(config.code_root())
    params = config.load(params_path)
    run = params.get(config.RUN_SECTION)
    if run is None:
        die(1, f"FATAL: {params_path} has no [run] section — it holds the driver's knobs "
            f"(bench_data, threads, methods, splits, …)")

    root_str = resolve(a.bench_data, "BENCH_DATA", run, "bench_data") or str(config.code_root())
    root = Path(root_str)
    if not root.is_dir():
        die(1, f"FATAL: BENCH_DATA does not exist: {root_str}")
    os.environ["BENCH_DATA"] = root_str

    max_ds_raw = resolve(a.max_ds, "MAX_DS_PER_SPLIT", run, "max_ds_per_split")
    try:
        max_ds = int(max_ds_raw)
    except ValueError:
        die(1, f"FATAL: max_ds_per_split must be an integer, got {max_ds_raw!r}")

    codes = resolve(a.methods, "RUN_METHODS", run, "methods").split()
    unknown = [c for c in codes if c not in METHOD_BY_CODE]
    if unknown:
        die(1, f"FATAL: unknown method code(s): {' '.join(unknown)} "
            f"(have: {' '.join(m.code for m in METHODS)})")
    scodes = resolve(a.splits, "RUN_SPLITS", run, "splits").split()
    unknown = [c for c in scodes if c not in SPLIT_BY_CODE]
    if unknown:
        die(1, f"FATAL: unknown split code(s): {' '.join(unknown)} "
            f"(have: {' '.join(SPLIT_BY_CODE)})")

    require_tools(BASE_TOOLS)
    bwa = shutil.which("bwa-mem2") or "bwa-mem2"

    ctx = Ctx(
        root=root, params_path=Path(params_path), params=params,
        threads=resolve(a.threads, "THREADS", run, "threads"),
        timeout=resolve(None, "CLASSIFY_TIMEOUT", run, "classify_timeout"),
        max_ds=max_ds,
        skip_build=resolve(a.skip_build, "SKIP_BUILD", run, "skip_build") == "1",
        kmcp_threads=resolve(a.kmcp_threads, "KMCP_THREADS", run, "kmcp_threads"),
        kmcp_allow_synthetic=resolve(a.kmcp_allow_synthetic, "KMCP_ALLOW_SYNTHETIC", run,
                                     "kmcp_allow_synthetic") == "1",
        dry_run=a.dry_run, keep_going=a.keep_going, bwa=bwa,
    )
    return ctx, [METHOD_BY_CODE[c] for c in codes], [SPLIT_BY_CODE[c] for c in scodes]


def main() -> int:
    sys.stdout.reconfigure(line_buffering=True)
    a = parse_args()
    ctx, methods, splits = build_ctx(a)
    if a.build_only and ctx.skip_build:
        die(1, "FATAL: --build-only with skip_build set does nothing")

    stamp = time.strftime("%F %T")
    print(f"=== PREMISE benchmark | {stamp} | threads={ctx.threads} | "
          f"classify timeout={ctx.timeout}s ===")
    print(f"    methods: {' '.join(m.code for m in methods)}")
    scope = f" (first {ctx.max_ds} dataset(s) each)" if ctx.max_ds > 0 else ""
    print(f"    splits:  {' '.join(s.code for s in splits)}{scope}")
    if ctx.dry_run:
        print("    dry run: no command below is executed")

    bins: list[str] = []
    for m in methods:
        bins.append(m.tool)
        bins.extend(m.extra_bins)
    require_tools(bins)

    prepare_real(ctx, splits)
    add_decoys(ctx)

    if not ctx.dry_run:
        (ctx.root / "results").mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ctx.params_path, ctx.root / "results/params.used.toml")

    failed: list[str] = []
    if ctx.skip_build:
        print("== [1/3] Building indexes — SKIPPED (SKIP_BUILD set; reusing existing DBs "
              "+ db_build.csv) ==")
    else:
        print("== [1/3] Building indexes ==")
        if not ctx.dry_run:
            (ctx.root / DB_BUILD_CSV).write_text(DB_BUILD_HEADER + "\n")
        for m in methods:
            print(f" building {m.tool}…")
            if build_index(ctx, m) != 0:
                failed.append(m.tool)
    if failed and not ctx.keep_going:
        die(1, "", f"FATAL: index build failed for: {' '.join(failed)}",
            "       Their index directories are now empty, so every classification against them "
            "would fail.",
            "       See the build logs, or pass --keep-going to run anyway.")

    if a.build_only:
        print(f"=== done (build only) | {DB_BUILD_CSV} ===")
        if failed:
            print(f"    NOTE: index build failed for {' '.join(failed)}")
            return 1
        return 0

    print("== [2/3] Classifying ==")
    for sp in splits:
        classify_split(ctx, methods, sp)

    print("== [3/3] Scoring ==")
    _analyze(ctx, [s.code for s in splits])

    print(f"=== done | {DB_BUILD_CSV} + results/tables/comparative-<split>.csv ===")
    if failed:
        print(f"    NOTE: index build failed for {' '.join(failed)}; "
              f"their results are not trustworthy")
        return 1
    return 0
