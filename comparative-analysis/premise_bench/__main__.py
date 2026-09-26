"""`python3 -m premise_bench <command> [options]` """
from __future__ import annotations

import importlib
import sys

COMMANDS = {
    "run":              ("pipeline", "main", "build every index, classify every split, score them"),
    "analyze":          ("evaluate.analyze", "main", "score results into results/tables/comparative-<split>.csv"),
    "ablation":         ("ablation.sweep", "main", "one-at-a-time sweep of PREMISE's parameters"),
    "fetch-real":       ("data.fetch_real", "main", "download the real samples from NCBI SRA"),
    "prepare-real":     ("data.prepare_real", "main", "trim, filter and derive truth for the real samples"),
    "read-truth":       ("data.read_truth", "main", "per-read truth from a bwa-mem2 stream (pipe stage)"),
    "synth":            ("data.synthetic", "main", "regenerate the synthetic datasets with InSilicoSeq"),
    "clean-db":         ("data.clean_db", "main", "drop duplicate and contained reference records"),
    "normalize":        ("methods", "normalize_main", "reshape one tool's native output by hand"),
    "params":           ("config", "main", "print params.toml as shell assignments"),
}


def usage() -> str:
    width = max(map(len, COMMANDS))
    lines = ["usage: python3 -m premise_bench <command> [options]", "", "commands:"]
    lines += [f"  {name:<{width}}  {summary}" for name, (_, _, summary) in COMMANDS.items()]
    return "\n".join(lines)


def main() -> int:
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        print(usage())
        return 0 if len(sys.argv) >= 2 else 2
    cmd = sys.argv[1]
    if cmd not in COMMANDS:
        print(f"premise_bench: unknown command '{cmd}'\n\n{usage()}", file=sys.stderr)
        return 2
    module, entry, _ = COMMANDS[cmd]
    sys.argv = [f"premise_bench {cmd}", *sys.argv[2:]]
    rc = getattr(importlib.import_module(f"premise_bench.{module}"), entry)()
    return rc or 0


if __name__ == "__main__":
    sys.exit(main())
