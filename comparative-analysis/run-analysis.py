#!/usr/bin/env python3
"""run-analysis.py — PREMISE comparative benchmark, end-to-end driver.

    nix develop ..#benchmark
    python3 run-analysis.py                       # build indexes + all methods + all splits
    python3 run-analysis.py --threads 32
    python3 run-analysis.py --methods "pre cen syl" --splits syn-iso --max-ds 1 --skip-build
    python3 run-analysis.py --dry-run             # print every command, run nothing

Three stages: build one index per method, classify every dataset of every split and evaluate it, then render the LaTeX comparison tables. All configuration lives in params.toml
"""
from __future__ import annotations

import argparse
import contextlib
import os
import shutil
import subprocess
import sys
import time
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

CODE = Path(__file__).resolve().parent
SCRIPTS = CODE / "scripts"

sys.path.insert(0, str(SCRIPTS))

import analyze
import load_params
import normalize_ganon
import normalize_karp
import normalize_mora
import normalize_sylph
import prepare_real_samples
from utils import die, nonempty

TIME = "/usr/bin/time"
KILL_GRACE = "30"
DB_BUILD_CSV = "db_build.csv"
DB_BUILD_HEADER = "method,size_before_bytes,size_after_bytes,build_seconds,exit"

BASE_TOOLS = ("bwa-mem2", "samtools", "cutadapt", "python3", "iss", "fasterq-dump")


@dataclass(frozen=True)
class Ctx:
    """Resolved configuration plus the two roots. `root` is $BENCH_DATA; every tool path below is
    relative to it, because every child runs with cwd=root."""
    root: Path
    params_path: Path
    params: dict[str, dict[str, str]]
    threads: str
    timeout: str
    max_ds: int
    skip_build: bool
    kmcp_threads: str
    kmcp_allow_synthetic: bool
    dry_run: bool
    keep_going: bool
    bwa: str

    def p(self, section: str, key: str) -> str:
        """One method parameter, e.g. p("pre", "mem"). A missing section or key is fatal."""
        try:
            return self.params[section][key]
        except KeyError:
            die(1, f"FATAL: {self.params_path}: no {section}.{key}")
            raise


def resolve(cli: str | None, env: str, params: dict[str, str], key: str) -> str:
    """CLI flag > environment variable > params.toml.
    """
    if cli is not None:
        return cli
    if env in os.environ:
        return os.environ[env]
    if key not in params:
        die(1, f"FATAL: params.toml [run] is missing '{key}' and neither --{key.replace('_', '-')} "
            f"nor ${env} was given")
    return params[key]


def run_cmd(ctx: Ctx, argv: list[str], *, out: str | None = None, err: str | None = None,
            err_append: bool = False, err_to_out: bool = False, err_null: bool = False) -> int:
    """Run one child with cwd=$BENCH_DATA, mirroring the shell's redirections.
    """
    if ctx.dry_run:
        redir = ""
        if out:
            redir += f" > {out}"
        if err_to_out:
            redir += " 2>&1"
        elif err:
            redir += f" 2>{'>' if err_append else ''} {err}"
        elif err_null:
            redir += " 2> /dev/null"
        print("   $ " + " ".join(argv) + redir)
        return 0

    with contextlib.ExitStack() as stack:
        fh_out = stack.enter_context(open(ctx.root / out, "w")) if out else None
        fh_err = None
        if err_to_out:
            fh_err = subprocess.STDOUT
        elif err:
            fh_err = stack.enter_context(open(ctx.root / err, "a" if err_append else "w"))
        elif err_null:
            fh_err = subprocess.DEVNULL
        return subprocess.call(argv, cwd=ctx.root, stdout=fh_out, stderr=fh_err)


def timed(ctx: Ctx, tm: str, argv: list[str], **redir) -> int:
    """The shell's TW(): cap the wall clock and record /usr/bin/time -v to <tm>.
    """
    return run_cmd(ctx, ["timeout", "-k", KILL_GRACE, ctx.timeout,
                         TIME, "-v", "-o", tm] + argv, **redir)


def require_tools(tools: tuple[str, ...] | list[str]) -> None:
    """Fail in the first second rather than three hours in, and name the fix."""
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        die(1, f"FATAL: not on PATH: {' '.join(missing)}",
            "       These come from the pinned toolchain — enter it first:",
            "         nix develop ..#benchmark")


@dataclass(frozen=True)
class Split:
    code: str
    sample_dir: str
    reads: str

    @property
    def res_sub(self) -> str:
        """results/ mirrors samples/ ({real,synthetic}/{isolate,mixed}), derived from sample_dir so
        the mapping lives in one place."""
        return self.sample_dir[len("samples/"):]

    @property
    def synthetic(self) -> bool:
        return self.sample_dir.startswith("samples/synthetic/")


SPLITS = (
    Split("syn-iso", "samples/synthetic/isolate", "reads_R{n}.fastq"),
    Split("syn-mix", "samples/synthetic/mixed", "reads_R{n}.fastq"),
    Split("syn-mix-subtype", "samples/synthetic/mixed-subtype", "reads_R{n}.fastq"),
    Split("real-iso", "samples/real/isolate", "{base}_{n}-filtered.ca.fastq"),
    Split("real-mix", "samples/real/mixed", "{base}_{n}-filtered.ca.fastq"),
)
SPLIT_BY_CODE = {s.code: s for s in SPLITS}


@dataclass(frozen=True)
class Job:
    """One (method, split, dataset) unit of classification work. Paths are relative to the root."""
    split: Split
    base: str
    outdir: str
    r1: str
    r2: str

    @property
    def ob(self) -> str:
        return self.base

    @property
    def tm(self) -> str:
        return f"{self.outdir}/time-mem"


CLEANED = "indexes/sequences-cleaned.fasta"

def fresh(ctx: Ctx, d: str) -> None:
    """rm -rf + mkdir -p. Destructive by design: a build must not merge into a stale index."""
    if ctx.dry_run:
        print(f"   $ rm -rf {d} && mkdir -p {d}")
        return
    shutil.rmtree(ctx.root / d, ignore_errors=True)
    (ctx.root / d).mkdir(parents=True, exist_ok=True)


def build_pre(ctx: Ctx) -> list[int]:
    fresh(ctx, "indexes/premise")
    return [run_cmd(ctx, ["premise", "build", "-s", CLEANED,
                          "-o", "indexes/premise/sequences.fmidx"])]


def build_kmp(ctx: Ctx) -> list[int]:
    fresh(ctx, "indexes/kmcp")
    if not ctx.dry_run:
        (ctx.root / "indexes/kmcp/kmers").mkdir(parents=True, exist_ok=True)
    return [
        run_cmd(ctx, ["kmcp", "compute", "--by-seq", "-O", "indexes/kmcp/kmers/",
                      "-k", ctx.p("kmp", "k"), "-j", ctx.threads, "--force", CLEANED]),
        run_cmd(ctx, ["kmcp", "index", "-I", "indexes/kmcp/kmers/",
                      "-O", "indexes/kmcp/sequences.kmcp", "-f", ctx.p("kmp", "f"),
                      "-n", ctx.p("kmp", "n"), "-j", ctx.threads, "--force"]),
    ]


def build_cen(ctx: Ctx) -> list[int]:
    fresh(ctx, "indexes/centrifuger")
    return [run_cmd(ctx, ["centrifuger-build", "--name-table", "indexes/names.dmp",
                          "--taxonomy-tree", "indexes/nodes.dmp", "-r", CLEANED,
                          "--conversion-table", "indexes/seqid2taxid.map",
                          "-o", "indexes/centrifuger/sequences-custom"])]


def build_mor(ctx: Ctx) -> list[int]:
    fresh(ctx, "indexes/mora")
    return [run_cmd(ctx, [ctx.bwa, "index", "-p", "indexes/mora/sequences", CLEANED])]


def build_kap(ctx: Ctx) -> list[int]:
    fresh(ctx, "indexes/karp")
    rcs = [run_cmd(ctx, ["samtools", "faidx", CLEANED])]
    if ctx.dry_run:
        print("   $ indexes/seqid2taxid.map -> indexes/sequences.tax  (Root;taxid_<N>; rows)")
    else:
        with open(ctx.root / "indexes/seqid2taxid.map") as f, \
             open(ctx.root / "indexes/sequences.tax", "w") as o:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) == 2:
                    o.write(parts[0] + "\tRoot;taxid_" + parts[1] + ";\n")
    rcs.append(run_cmd(ctx, ["karp", "-c", "index", "-r", CLEANED,
                             "-i", "indexes/karp/sequences.index", "-k", ctx.p("kap", "k")]))
    return rcs


def build_syl(ctx: Ctx) -> list[int]:
    fresh(ctx, "indexes/sylph")
    return [run_cmd(ctx, ["sylph", "sketch", CLEANED, "-i", "-c", ctx.p("syl", "c"),
                          "-o", "indexes/sylph/sequences", "-t", ctx.threads],
                    err="indexes/sylph/build.log")]


def build_gan(ctx: Ctx) -> list[int]:
    fresh(ctx, "indexes/ganon")
    abs_fasta = str(ctx.root / CLEANED)
    if ctx.dry_run:
        print(f"   $ indexes/seqid2taxid.map -> indexes/ganon/input.tsv  ({abs_fasta}\\t$1\\t$2)")
    else:
        with open(ctx.root / "indexes/seqid2taxid.map") as f, \
             open(ctx.root / "indexes/ganon/input.tsv", "w") as o:
            for line in f:
                fields = line.split()
                f1 = fields[0] if len(fields) > 0 else ""
                f2 = fields[1] if len(fields) > 1 else ""
                o.write(f"{abs_fasta}\t{f1}\t{f2}\n")
    return [run_cmd(ctx, ["ganon", "build-custom", "--input-file", "indexes/ganon/input.tsv",
                          "--input-target", "sequence", "--taxonomy", "ncbi",
                          "--taxonomy-files", "indexes/nodes.dmp", "indexes/names.dmp",
                          "--db-prefix", "indexes/ganon/sequences", "--skip-genome-size",
                          "--threads", ctx.threads],
                    err="indexes/ganon/build.log")]


def classify_pre(ctx: Ctx, j: Job) -> int:
    iters = ctx.p("pre", "iters_real") if j.split.code == "real-iso" else ctx.p("pre", "iters")
    return timed(ctx, j.tm, ["premise", "query", "-s", "indexes/premise/sequences.fmidx",
                             "-t", ctx.p("pre", "threads"), "-1", j.r1, "-2", j.r2,
                             "-o", f"{j.outdir}/{j.ob}", "-m", ctx.p("pre", "mem"), "-i", iters,
                             "--eps_1", ctx.p("pre", "eps_1"), "--eps_2", ctx.p("pre", "eps_2"),
                             "--em_threshold", ctx.p("pre", "em_threshold"),
                             "--rho", ctx.p("pre", "rho"), "--omega", ctx.p("pre", "omega")],
                 out=f"{j.outdir}/{j.ob}.log", err_to_out=True)


def classify_kmp(ctx: Ctx, j: Job) -> int:
    # Three choices:
    #   -w   --load-whole-db: read the index into contiguous memory instead of mmap'ing it.
    #        kmcp recommends it for small databases, and this one is ~3.5 MB (db_build.csv's
    #        size_after_bytes for kmcp), so the extra memory is free. Purely a speed knob.
    #   .gz  the raw match dump is ~18 GB/dataset; writing it uncompressed is what blows past the
    #        timeout. profile + reduce both read .gz transparently.
    dump = f"{j.outdir}/{j.ob}.tsv.gz"
    rc = timed(ctx, j.tm, ["kmcp", "search", "-w", "-j", ctx.kmcp_threads,
                           "-d", "indexes/kmcp/sequences.kmcp", "-1", j.r1, "-2", j.r2,
                           "-o", dump], err=f"{j.outdir}/{j.ob}.log")
    if rc != 0:
        print(f"    kmcp search failed (exit {rc}) — keeping partial dump as "
              f"{j.ob}.tsv.gz.partial, no profile")
        if not ctx.dry_run:
            src = ctx.root / dump
            if src.exists():
                src.replace(ctx.root / f"{dump}.partial")
        return rc
    profile = f"{j.outdir}/{j.ob}.profile"
    run_cmd(ctx, ["kmcp", "profile", "--level", ctx.p("kmp", "profile_level"),
                  "-m", ctx.p("kmp", "profile_m"), dump, "-o", profile],
            err=f"{j.outdir}/{j.ob}.log", err_append=True)
    kmcp_pop_prop(ctx, profile, f"{j.outdir}/{j.ob}.pop-prop")

    return rc


def kmcp_pop_prop(ctx: Ctx, profile: str, dest: str) -> None:
    if ctx.dry_run:
        print(f"   $ {profile} -> {dest}  (fields 1 and 8)")
        return
    src = ctx.root / profile
    with open(ctx.root / dest, "w") as o:
        if not src.exists():
            return
        with open(src) as f:
            for line in f:
                fields = line.split()
                f1 = fields[0] if len(fields) > 0 else ""
                f8 = fields[7] if len(fields) > 7 else ""
                o.write(f"{f1} {f8}\n")

def classify_cen(ctx: Ctx, j: Job) -> int:
    return timed(ctx, j.tm, ["centrifuger", "-t", ctx.threads,
                             "-x", "indexes/centrifuger/sequences-custom", "-1", j.r1, "-2", j.r2],
                 out=f"{j.outdir}/{j.ob}.tsv", err=f"{j.outdir}/{j.ob}.log")

def classify_mor(ctx: Ctx, j: Job) -> int:
    sam = f"{j.outdir}/{j.ob}.sam"
    rc = timed(ctx, j.tm, [ctx.bwa, "mem", "-t", ctx.threads, "indexes/mora/sequences",
                           j.r1, j.r2], out=sam, err=f"{j.outdir}/{j.ob}.bwa.log")
    log = f"{j.outdir}/{j.ob}.mora.log"
    run_cmd(ctx, ["mora", "-s", sam, "-o", f"{j.outdir}/{j.ob}.txt"], err=log)
    _normalize(ctx, "mora", f"{j.outdir}/{j.ob}.txt", f"{j.outdir}/{j.ob}.tsv", log)
    if ctx.dry_run:
        print(f"   $ rm -f {sam}")
    else:
        (ctx.root / sam).unlink(missing_ok=True)
    return rc


def classify_kap(ctx: Ctx, j: Job) -> int:
    log = f"{j.outdir}/{j.ob}.log"
    rc = timed(ctx, j.tm, ["karp", "-c", "quantify", "-r", CLEANED,
                           "-i", "indexes/karp/sequences.index", "-f", j.r1, "-q", j.r2,
                           "--paired", "-t", "indexes/sequences.tax", "-o", f"{j.outdir}/{j.ob}",
                           "--threads", ctx.threads, "--readinfo", "--no_harp_filter",
                           "--min_freq", ctx.p("kap", "min_freq")], err=log)
    _normalize(ctx, "karp", f"{j.outdir}/{j.ob}_readinfo.txt.gz",
               f"{j.outdir}/{j.ob}.tsv", log)
    return rc


def classify_gan(ctx: Ctx, j: Job) -> int:
    log = f"{j.outdir}/{j.ob}.log"
    rc = timed(ctx, j.tm, ["ganon", "classify", "-d", "indexes/ganon/sequences",
                           "-p", j.r1, j.r2, "-o", f"{j.outdir}/{j.ob}", "--output-one",
                           "--multiple-matches", "em", "--rel-cutoff",
                           ctx.p("gan", "rel_cutoff"), "--fpr-query", ctx.p("gan", "fpr_query"),
                           "-t", ctx.threads], err=log)
    _normalize(ctx, "ganon", f"{j.outdir}/{j.ob}.one", f"{j.outdir}/{j.ob}.tsv", log)
    return rc


def classify_syl(ctx: Ctx, j: Job) -> int:
    log = f"{j.outdir}/{j.ob}.log"
    rc = timed(ctx, j.tm, ["sylph", "profile", "indexes/sylph/sequences.syldb", j.r1, j.r2,
                           "-t", ctx.threads, "-c", ctx.p("syl", "c"),
                           "--min-number-kmers", ctx.p("syl", "min_kmers")],
               out=f"{j.outdir}/{j.ob}.tsv", err=log)
    _normalize(ctx, "sylph", f"{j.outdir}/{j.ob}.tsv",
               f"{j.outdir}/{j.ob}.norm.tsv", log)
    return rc


NORMALIZERS = {
    "mora": normalize_mora.normalize,
    "karp": normalize_karp.normalize,
    "ganon": normalize_ganon.normalize,
    "sylph": normalize_sylph.normalize,
}


def _normalize(ctx: Ctx, which: str, src: str, dest: str, log: str) -> None:
    """Reshape one tool's native output into the TSV analyze.py reads.
    """
    if ctx.dry_run:
        print(f"   $ normalize_{which}: {src} -> {dest}")
        return
    try:
        NORMALIZERS[which](ctx.root / src, ctx.root / dest)
    except Exception as e:
        with open(ctx.root / log, "a") as fh:
            fh.write(f"normalize_{which}: {e}\n")



@dataclass(frozen=True)
class Method:
    """One compared tool, defined once. The name is also the binary required on $PATH."""
    code: str
    name: str
    index_paths: tuple[str, ...]
    build: Callable[[Ctx], list[int]]
    classify: Callable[[Ctx, Job], int]
    extra_bins: tuple[str, ...] = ()
    real_only: bool = False


METHODS = (
    Method("pre", "premise", ("indexes/premise/sequences.fmidx",), build_pre, classify_pre),
    # kmcp is real-splits-only: on the 1M-pair synthetic datasets `kmcp search` runs at ~4.6k
    # queries/min (~3.6 h/dataset), so it cannot finish under the classification cap. No results
    # dir is created, so analyze.py renders '---' for kmcp on synthetic/mixed. Escape hatch:
    # --kmcp-allow-synthetic (KMCP_ALLOW_SYNTHETIC=1) forces it for diagnostic reruns.
    Method("kmp", "kmcp", ("indexes/kmcp/sequences.kmcp",), build_kmp, classify_kmp,
           real_only=True),
    Method("cen", "centrifuger",
           tuple(f"indexes/centrifuger/sequences-custom.{i}.cfr" for i in (1, 2, 3, 4)),
           build_cen, classify_cen),
    Method("mor", "mora",
           ("indexes/mora/sequences.0123", "indexes/mora/sequences.amb",
            "indexes/mora/sequences.ann", "indexes/mora/sequences.bwt.2bit.64",
            "indexes/mora/sequences.pac"), build_mor, classify_mor),
    Method("kap", "karp", ("indexes/karp/sequences.index",), build_kap, classify_kap),
    Method("gan", "ganon",
           ("indexes/ganon/sequences.hibf", "indexes/ganon/sequences.ibf",
            "indexes/ganon/sequences.tax"), build_gan, classify_gan,
           extra_bins=("ganon-build", "ganon-classify", "raptor")),
    Method("syl", "sylph", ("indexes/sylph/sequences.syldb",), build_syl, classify_syl),
)
METHOD_BY_CODE = {m.code: m for m in METHODS}



def idx_size(ctx: Ctx, m: Method) -> int:
    """On-disk bytes of a method's index; 0 when it is absent."""
    total = 0
    for rel in m.index_paths:
        p = ctx.root / rel
        if p.exists():
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
            f.write(f"{m.name},{before},{after},{secs:.9f},{rc}\n")
    print(f"   built {m.name} in {secs:.9f}s (index {after / 1e9:.3f} GB, exit {rc})")
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
            r1 = d / sp.reads.format(base=d.name, n=1)
            if not (nonempty(r1) and nonempty(d / "truth_assignments.tsv")):
                need = True
        if not need:
            continue
        print(f" preparing {code} reads + truth (scripts/prepare_real_samples.py)…")
        if ctx.dry_run:
            print(f"   $ prepare_real_samples.run({code})")
            continue
        try:
            prepare_real_samples.run(code, int(ctx.threads))
        except BaseException as e:  # noqa: BLE001 — prepare die()s with SystemExit on a bad tree
            print(f" WARN: prepare_real_samples failed ({e}); {code} may be incomplete")


def make_job(ctx: Ctx, m: Method, sp: Split, ds: Path) -> Job | None:
    base = ds.name
    if m.real_only and sp.synthetic and not ctx.kmcp_allow_synthetic:
        print(f"    skip {m.name} ({sp.code}/{base}): synthetic splits are "
              f"{m.name}-free by design")
        return None
    outdir = f"results/{m.name}/{sp.res_sub}/{base}"
    if not ctx.dry_run:
        (ctx.root / outdir).mkdir(parents=True, exist_ok=True)
    r1, r2 = (f"{sp.sample_dir}/{base}/" + sp.reads.format(base=base, n=n) for n in (1, 2))
    if not (nonempty(ctx.root / r1) and nonempty(ctx.root / r2)):
        print(f"    skip {m.name} ({sp.code}/{base}): missing reads")
        return None
    return Job(sp, base, outdir, r1, r2)


def _analyze(ctx: Ctx, splits: list[str], datasets: list[str] | None, latex: bool,
             outfile: str | None = None) -> None:
    """Score and print, in-process. Failure is non-fatal, exactly as the subprocess call was:
    one dataset that cannot be scored must not abort a three-hour run."""
    if ctx.dry_run:
        ds = f" --datasets {','.join(datasets)}" if datasets else ""
        print(f"   $ analyze.run(--splits {','.join(splits)}{ds})")
        return
    try:
        print(analyze.run(splits, datasets,
                          ctx.root / outfile if outfile else None, latex=latex))
    except Exception as e:                                  # noqa: BLE001 — see docstring
        print(f"analyze.py: {e}", file=sys.stderr)


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
                print(f"    TIMEOUT ({ctx.timeout}s) {m.name} ({sp.code}/{ds.name}) — "
                      f"recorded as failed, continuing")
        print(f"   evaluating {sp.code}/{ds.name}…")
        _analyze(ctx, [sp.code], [ds.name], latex=False)
        count += 1
        if ctx.max_ds > 0 and count >= ctx.max_ds:
            print(f" (stopping {sp.code} at {ctx.max_ds} dataset(s))")
            break



def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--params", type=Path, default=None,
                    help="params.toml to use (default: the one beside this script)")
    ap.add_argument("--bench-data", default=None, help="data root; overrides $BENCH_DATA")
    ap.add_argument("--threads", default=None, help="worker threads for every method")
    ap.add_argument("--methods", default=None, help='whitespace-separated codes, e.g. "pre cen syl"')
    ap.add_argument("--splits", default=None, help="whitespace-separated split codes")
    ap.add_argument("--max-ds", default=None, help="datasets per split; 0 = all")
    ap.add_argument("--skip-build", action="store_const", const="1", default=None,
                    help="reuse existing indexes/ and db_build.csv")
    ap.add_argument("--kmcp-threads", default=None)
    ap.add_argument("--kmcp-allow-synthetic", action="store_const", const="1", default=None,
                    help="run kmcp on the synthetic splits too (diagnostic reruns)")
    ap.add_argument("--keep-going", action="store_true",
                    help="carry on after a failed index build instead of aborting")
    ap.add_argument("--dry-run", action="store_true", help="print every command, run nothing")
    return ap.parse_args(argv)


def build_ctx(a: argparse.Namespace) -> tuple[Ctx, list[Method], list[Split]]:
    params_path = a.params if a.params is not None else load_params.find_file(CODE)
    params = load_params.load(params_path)
    run = params.get(load_params.RUN_SECTION)
    if run is None:
        die(1, f"FATAL: {params_path} has no [run] section — it holds the driver's knobs "
            f"(bench_data, threads, methods, splits, …)")

    root_str = resolve(a.bench_data, "BENCH_DATA", run, "bench_data") or str(CODE)
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
            f"(have: {' '.join(s.code for s in SPLITS)})")

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
        bins.append(m.name)
        bins.extend(m.extra_bins)
    require_tools(bins)

    prepare_real(ctx, splits)

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
            print(f" building {m.name}…")
            if build_index(ctx, m) != 0:
                failed.append(m.name)
    if failed and not ctx.keep_going:
        die(1, "", f"FATAL: index build failed for: {' '.join(failed)}",
            "       Their index directories are now empty, so every classification against them "
            "would fail.",
            "       See the build logs, or pass --keep-going to run anyway.")

    print("== [2/3] Classifying + per-dataset evaluation ==")
    for sp in splits:
        classify_split(ctx, methods, sp)

    print("== [3/3] Generating LaTeX comparison tables ==")
    _analyze(ctx, [s.code for s in splits], None, latex=True,
             outfile="results/tables.generated.tex")

    print(f"=== done | {DB_BUILD_CSV} + results/tables.generated.tex ===")
    if failed:
        print(f"    NOTE: index build failed for {' '.join(failed)}; "
              f"their results are not trustworthy")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
