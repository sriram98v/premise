"""Benchmark configuration: params.toml and the code / data roots"""
from __future__ import annotations

import argparse
import os
import shlex
import sys
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:
    sys.exit("load_params: needs Python 3.11+ for the stdlib tomllib module "
             f"(running {sys.version_info.major}.{sys.version_info.minor})")

FILENAME = "params.toml"
RUN_SECTION = "run"

_ALLOWED = str

MARKERS = ("premise_bench", "run-analysis.py")
PARAMS = FILENAME

BENCH_DATA_ENV = "BENCH_DATA"


def code_root(start: Path | None = None) -> Path:
    """The benchmark code root
    """
    start = (start or Path(__file__)).resolve()
    for p in [start, *start.parents]:
        if (p / PARAMS).exists() and any((p / m).exists() for m in MARKERS):
            return p
    raise SystemExit(f"benchmark code root (the directory holding {PARAMS} and one of "
                     f"{'/'.join(MARKERS)}) not found above {start}")


def data_root_str() -> str:
    """The RAW $BENCH_DATA string, unnormalised, falling back to the code root.
    """
    return os.environ.get(BENCH_DATA_ENV) or str(code_root())


def data_root() -> Path:
    """The data root: $BENCH_DATA, or the code root when it is unset."""
    return Path(data_root_str())


def indexes() -> Path:
    return data_root() / "indexes"

CLEANED_FASTA = "sequences-cleaned.fasta"


def samples() -> Path:
    return data_root() / "samples"


def results() -> Path:
    return data_root() / "results"


def find_file(start: Path | None = None) -> Path:
    """Locate params.toml"""
    here = Path(__file__).resolve()
    candidates = [here.parent.parent / FILENAME, here.parent / FILENAME]
    if start is not None:
        start = Path(start).resolve()
        candidates += [p / FILENAME for p in [start, *start.parents]]
    for c in candidates:
        if c.is_file():
            return c
    raise SystemExit(f"load_params: {FILENAME} not found (looked in "
                     + ", ".join(str(c.parent) for c in candidates) + ")")


def load(path: Path | None = None) -> dict[str, dict[str, str]]:
    """Parse params.toml into {section: {key: str}}, rejecting anything non-scalar."""
    path = Path(path) if path is not None else find_file()
    with open(path, "rb") as fh:
        try:
            raw = tomllib.load(fh)
        except tomllib.TOMLDecodeError as exc:
            raise SystemExit(f"load_params: {path} is not valid TOML: {exc}") from exc

    out: dict[str, dict[str, str]] = {}
    for sect, body in raw.items():
        if not isinstance(body, dict):
            raise SystemExit(f"load_params: {path}: top-level key '{sect}' must be a "
                             "[section] table, not a bare value")
        vals: dict[str, str] = {}
        for key, val in body.items():
            if not isinstance(val, _ALLOWED):
                raise SystemExit(
                    f"load_params: {path}: {sect}.{key} must be a quoted string, got "
                    f"{type(val).__name__} ({val!r}). These values are command-line tokens and "
                    "must reach the method verbatim; TOML typing would rewrite them (1e-6 -> "
                    f'"1e-06"). Write it as {key} = "{val}".')
            vals[key] = val
        out[sect] = vals
    return out


def section(name: str, path: Path | None = None) -> dict[str, str]:
    """One method's parameters, by three-letter code. Missing section is fatal, not empty."""
    params = load(path)
    try:
        return params[name]
    except KeyError:
        raise SystemExit(f"load_params: no [{name}] section in {path or find_file()} "
                         f"(have: {', '.join(sorted(params)) or 'none'})") from None


def shell_assignments(params: dict[str, dict[str, str]]) -> list[str]:
    lines = []
    for sect in sorted(params):
        if sect == RUN_SECTION:
            continue
        for key in sorted(params[sect]):
            lines.append(f"{sect.upper()}_{key.upper()}={shlex.quote(params[sect][key])}")
    return lines


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--shell", action="store_true",
                    help="emit shell assignments for `eval` (default action)")
    ap.add_argument("--get", metavar="SECTION.KEY", default=None,
                    help="print one value, unquoted, for `: \"${VAR:=$(… --get run.threads)}\"`")
    ap.add_argument("--file", type=Path, default=None,
                    help=f"path to {FILENAME} (default: the one beside the premise_bench package)")
    args = ap.parse_args()

    params = load(args.file)
    if not params:
        return 1
    if args.get is not None:
        sect, _, key = args.get.partition(".")
        if not _ or not key:
            print(f"load_params: --get takes SECTION.KEY, got {args.get!r}", file=sys.stderr)
            return 2
        if sect not in params:
            print(f"load_params: no [{sect}] section (have: {', '.join(sorted(params))})",
                  file=sys.stderr)
            return 2
        if key not in params[sect]:
            print(f"load_params: no key '{key}' in [{sect}] "
                  f"(have: {', '.join(sorted(params[sect])) or 'none'})", file=sys.stderr)
            return 2
        print(params[sect][key])
        return 0
    print("\n".join(shell_assignments(params)))
    return 0
