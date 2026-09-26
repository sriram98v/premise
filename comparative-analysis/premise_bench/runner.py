"""Running tools for the benchmark"""
from __future__ import annotations

import contextlib
import os
import shutil
import subprocess
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from .config import CLEANED_FASTA
from .utils import die
from .splits import Split

KILL_GRACE = "30"


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
        fh_out, fh_err = _redirect(ctx, stack, out, err, err_append, err_to_out, err_null)
        return subprocess.call(argv, cwd=ctx.root, stdout=fh_out, stderr=fh_err)


def _redirect(ctx: Ctx, stack: contextlib.ExitStack, out, err, err_append, err_to_out, err_null):
    """The shell's redirections, as (stdout, stderr) handles."""
    fh_out = stack.enter_context(open(ctx.root / out, "w")) if out else None
    fh_err = None
    if err_to_out:
        fh_err = subprocess.STDOUT
    elif err:
        fh_err = stack.enter_context(open(ctx.root / err, "a" if err_append else "w"))
    elif err_null:
        fh_err = subprocess.DEVNULL
    return fh_out, fh_err


def timed(ctx: Ctx, tm: str, argv: list[str], **redir) -> int:
    """Cap the wall clock, and record the child's resource usage to <tm>"""
    full = ["timeout", "-k", KILL_GRACE, ctx.timeout] + argv
    if ctx.dry_run:
        return run_cmd(ctx, full, **redir)

    with contextlib.ExitStack() as stack:
        fh_out, fh_err = _redirect(ctx, stack, redir.get("out"), redir.get("err"),
                                   redir.get("err_append", False), redir.get("err_to_out", False),
                                   redir.get("err_null", False))
        t0 = os.times().elapsed
        proc = subprocess.Popen(full, cwd=ctx.root, stdout=fh_out, stderr=fh_err)
        _, status, ru = os.wait4(proc.pid, 0)
        wall = os.times().elapsed - t0

    rc = os.waitstatus_to_exitcode(status)
    proc.returncode = rc          # the child is already reaped; stop Popen waiting on it again
    (ctx.root / tm).write_text(
        f'\tCommand being timed: "{" ".join(argv)}"\n'
        f"\tUser time (seconds): {ru.ru_utime:.2f}\n"
        f"\tSystem time (seconds): {ru.ru_stime:.2f}\n"
        f"\tElapsed (wall clock) time (h:mm:ss or m:ss): {wall:.2f}\n"
        f"\tMaximum resident set size (kbytes): {ru.ru_maxrss}\n"
        f"\tExit status: {rc}\n")
    return rc


def require_tools(tools: tuple[str, ...] | list[str]) -> None:
    """Fail in the first second rather than three hours in, and name the fix."""
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        die(1, f"FATAL: not on PATH: {' '.join(missing)}",
            "       These come from the pinned toolchain — enter it first:",
            "         nix develop ..#benchmark")


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


CLEANED = f"indexes/{CLEANED_FASTA}"


def fresh(ctx: Ctx, d: str) -> None:
    """rm -rf + mkdir -p. Destructive by design: a build must not merge into a stale index."""
    if ctx.dry_run:
        print(f"   $ rm -rf {d} && mkdir -p {d}")
        return
    shutil.rmtree(ctx.root / d, ignore_errors=True)
    (ctx.root / d).mkdir(parents=True, exist_ok=True)


def normalize_output(ctx: Ctx, which: str, fn: Callable[[Path, Path], object], src: str,
                     dest: str, log: str) -> None:
    """Reshape one tool's native output into the TSV the evaluation reads, with `fn`"""
    if ctx.dry_run:
        print(f"   $ normalize_{which}: {src} -> {dest}")
        return
    try:
        fn(ctx.root / src, ctx.root / dest)
    except Exception as e:
        with open(ctx.root / log, "a") as fh:
            fh.write(f"normalize_{which}: {e}\n")
