"""The compared tool runners"""
from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from ..runner import Ctx, Job
from . import centrifuger, ganon, karp, kmcp, mora, premise, salmon, sylph


@dataclass(frozen=True)
class Method:
    code: str
    """Three-letter code used in params.toml and on the command line."""
    tool: str
    """The binary required on $PATH, which also names its indexes/ and results/ directories."""
    name: str
    """Display name, written to the CSVs' method_name column."""
    index_paths: tuple[str, ...]
    """Index files (or directories) whose size is reported as the database size."""
    build: Callable[[Ctx], list[int]]
    classify: Callable[[Ctx, Job], int]
    load: Callable[[Path, str], tuple[dict | None, dict | None]]
    extra_bins: tuple[str, ...] = ()
    real_only: bool = False
    """Skipped on the synthetic splits unless kmcp_allow_synthetic is set (see kmcp)."""
    profiler: bool = True
    perread: bool = True
    """Emits per-read assignments, so precision and coverage can be scored."""
    add_other: bool = True
    prec_rec: bool = True
    """Scored on precision and coverage."""
    complete_output: bool = False
    """Emits exactly one row per input read, so its row count must equal the read count."""

METHODS = (
    Method("pre", "premise", "Premise", ("indexes/premise/sequences.fmidx",),
           premise.build, premise.classify, premise.load, complete_output=True),
    Method("cen", "centrifuger", "Centrifuger",
           tuple(f"indexes/centrifuger/sequences-custom.{i}.cfr" for i in (1, 2, 3, 4)),
           centrifuger.build, centrifuger.classify, centrifuger.load, complete_output=True),
    Method("gan", "ganon", "Ganon",
           ("indexes/ganon/sequences.hibf", "indexes/ganon/sequences.ibf",
            "indexes/ganon/sequences.tax"),
           ganon.build, ganon.classify, ganon.load,
           extra_bins=("ganon-build", "ganon-classify", "raptor"),
           profiler=False, add_other=False),
    Method("kap", "karp", "Karp", ("indexes/karp/sequences.index",),
           karp.build, karp.classify, karp.load,
           perread=False, prec_rec=False),
    # kmcp is real-splits-only; on the 1M-pair synthetic datasets `kmcp search` runs at ~4.6k queries/min (~3.6 h/dataset)
    # --kmcp-allow-synthetic (KMCP_ALLOW_SYNTHETIC=1) forces it for diagnostic reruns.
    Method("kmp", "kmcp", "KMCP", ("indexes/kmcp/sequences.kmcp",),
           kmcp.build, kmcp.classify, kmcp.load, real_only=True),
    Method("mor", "mora", "MORA",
           ("indexes/mora/sequences.0123", "indexes/mora/sequences.amb",
            "indexes/mora/sequences.ann", "indexes/mora/sequences.bwt.2bit.64",
            "indexes/mora/sequences.pac"),
           mora.build, mora.classify, mora.load,
           profiler=False, add_other=False, complete_output=True),
    # salmon's index is a directory, not a set of named files; its size is summed recursively.
    Method("sal", "salmon", "Salmon", ("indexes/salmon/sequences",),
           salmon.build, salmon.classify, salmon.load),
    Method("syl", "sylph", "Sylph", ("indexes/sylph/sequences.syldb",),
           sylph.build, sylph.classify, sylph.load,
           perread=False, prec_rec=False),
)
METHOD_BY_CODE = {m.code: m for m in METHODS}

NORMALIZERS = {
    "mora": mora.normalize,
    "karp": karp.normalize,
    "ganon": ganon.normalize,
    "sylph": sylph.normalize,
    "salmon": salmon.normalize,
    "salmon_props": salmon.props,
}


def normalize_main() -> int:
    """`premise_bench normalize TOOL -i IN -o OUT`: reshape one output file by hand."""
    import argparse
    import sys

    ap = argparse.ArgumentParser(description=normalize_main.__doc__)
    ap.add_argument("tool", choices=sorted(NORMALIZERS),
                    help="salmon_props reads quant.sf and writes proportions")
    ap.add_argument("-i", "--input", required=True, type=Path)
    ap.add_argument("-o", "--output", required=True, type=Path)
    a = ap.parse_args()
    try:
        NORMALIZERS[a.tool](a.input, a.output)
    except (ValueError, OSError) as e:
        print(f"normalize_{a.tool}: {e}", file=sys.stderr)
        return 1
    return 0
