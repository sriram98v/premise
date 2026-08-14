#!/usr/bin/env python3
"""Render the PREMISE parameter ablation (ablation.py) as text + LaTeX tables.

Usage:
    python3 scripts/ablation_report.py
    python3 scripts/ablation_report.py --split syn-mix
    python3 scripts/ablation_report.py --outfile results/ablation-tables.tex
    python3 scripts/ablation_report.py --no-latex
"""
from __future__ import annotations

import argparse
import csv
import statistics
import sys
from pathlib import Path

import utils


def default_tex() -> Path:
    return utils.results() / "ablation-tables.generated.tex"


SPLITS = {
    "syn-iso":  ("synthetic/isolate", "synthetic single-strain"),
    "syn-mix":  ("synthetic/mixed", "synthetic mixed-infection"),
    "real-iso": ("real/isolate", "real single-strain"),
    "real-mix": ("real/mixed", "real mixed-infection"),
}
SPLIT_ORDER = ["syn-iso", "syn-mix", "real-iso", "real-mix"]

COLUMNS = [
    ("wall_s",    "Time (s)",   4,    True),
    ("rss_gb",    "Mem (GB)",   3,    True),
    ("precision", "Precision",  4,    False),
    ("coverage",  "Coverage",   None, False),
    ("ruzicka",   "Ruzicka",    4,    True),
    ("jaccard",   "Jaccard",    4,    True),
]

PARAM_META = {
    "mem":   (r"\texttt{-p}",       "minimum SMEM seed length"),
    "eps_2": (r"$\varepsilon_2$",   "minimum per-alignment match probability"),
    "eps_1": (r"$\varepsilon_1$",   "likelihood cutoff for dropping alignments before EM"),
    "rho":   (r"$\rho$",            "$L_1$ penalty weight in the EM"),
    "omega": (r"$\omega$",          "$L_1$ penalty floor in the EM"),
}
PARAM_ORDER = ["mem", "eps_2", "eps_1", "rho", "omega"]
PLAIN_LABEL = {"mem": "-p (seed len)", "eps_2": "eps_2", "eps_1": "eps_1",
               "rho": "rho", "omega": "omega"}


def fmt(v, sig):
    if v is None or v != v:
        return "-"
    return f"{float(v):.{sig}g}" if sig else f"{float(v):.2f}"


def load(split: str):
    """param -> [aggregated row dicts] in grid order, each averaged over the split's samples."""
    csv_path = utils.data_root() / "results" / "ablation" / SPLITS[split][0] / "ablation.csv"
    if not csv_path.exists():
        return None
    groups: dict[tuple, list] = {}
    order: list[tuple] = []
    with open(csv_path) as f:
        for r in csv.DictReader(f):
            k = (r["param"], r["value"])
            if k not in groups:
                groups[k] = []
                order.append(k)
            groups[k].append(r)
    by_param: dict[str, list] = {}
    for param, value in order:
        rs = groups[(param, value)]
        row = {"value": value, "baseline": rs[0]["baseline"] == "1", "n": len(rs)}
        for key, _, _, _ in COLUMNS:
            vals = [float(r[key]) for r in rs if r[key] not in ("", None)]
            row[key] = statistics.fmean(vals) if vals else None
            row[key + "_sd"] = statistics.stdev(vals) if len(vals) > 1 else 0.0
        refs = [int(r["n_refs"]) for r in rs if r["n_refs"]]
        row["n_refs"] = statistics.fmean(refs) if refs else 0
        by_param.setdefault(param, []).append(row)
    return by_param


def best_rows(rows):
    """{column key: index of the best mean} — ties resolved to the first occurrence."""
    best = {}
    for key, _, _, lower in COLUMNS:
        vals = [(i, r[key]) for i, r in enumerate(rows) if r[key] is not None]
        if vals:
            best[key] = (min if lower else max)(vals, key=lambda iv: iv[1])[0]
    return best


def render(split: str, param: str, rows):
    label, desc = PARAM_META[param]
    _, caption = SPLITS[split]
    tag = split          # \label tag is the split code: tab:ablation-real-mix-mem
    best = best_rows(rows)
    n = rows[0]["n"] if rows else 0

    out = [f"\n=== {caption} / {PLAIN_LABEL[param]}: "
           f"{desc.replace('$', '').replace(chr(92), '')}  (mean +/- SD, n={n}) ===",
           f"{'value':>10}" + "".join(f"{h:>21}" for _, h, _, _ in COLUMNS) + f"{'refs':>7}"]
    for r in rows:
        line = f"{r['value'] + ('*' if r['baseline'] else ' '):>10}"
        for key, _, sig, _ in COLUMNS:
            line += f"{fmt(r[key], sig) + ' ± ' + fmt(r[key + '_sd'], sig or 2):>21}"
        out.append(line + f"{r['n_refs']:>7.1f}")
    out.append("  (* = production setting)")

    L = ["\\begin{table}",
         f"\\caption{{PREMISE ablation on the {caption} split: {desc} ({label}). "
         f"Mean $\\pm$ SD over {n} samples. The production setting is marked $\\dagger$."
         f"\\label{{tab:ablation-{tag}-{param.replace('_', '-')}}}}}",
         "\\begin{tabular*}{\\columnwidth}{@{\\extracolsep{\\fill}}l" + "c" * len(COLUMNS) + "@{}}",
         "\\toprule",
         f"\\textbf{{{label}}}" + "".join(f" & \\textbf{{{h}}}" for _, h, _, _ in COLUMNS)
         + " \\\\ \\midrule"]
    for i, r in enumerate(rows):
        cells = f"{r['value']}" + ("$^\\dagger$" if r["baseline"] else "")
        for key, _, sig, _ in COLUMNS:
            s = fmt(r[key], sig)
            if best.get(key) == i and s != "-":
                s = f"\\textbf{{{s}}}"
            cells += f" & {s}\\,\\tiny{{$\\pm$\\,{fmt(r[key + '_sd'], sig or 2)}}}"
        L.append(cells + " \\\\")
    L += ["\\bottomrule", "\\end{tabular*}", "\\end{table}", ""]
    return "\n".join(out), "\n".join(L)


def run(splits_arg: list[str] | None = None, outfile: Path | None = None,
        latex: bool = True) -> int:
    """Render the ablation tables to stdout and, when `latex`, to `outfile`."""
    outfile = Path(outfile) if outfile else default_tex()

    splits = [s for s in SPLIT_ORDER if s in (splits_arg or SPLIT_ORDER)]
    tex = ["% Generated by scripts/ablation_report.py -- do not edit by hand.", ""]
    any_rows = False
    for split in splits:
        by_param = load(split)
        if by_param is None:
            print(f"{split}: no results yet, skipped")
            continue
        for p in PARAM_ORDER:
            if p not in by_param:
                continue
            txt, latex = render(split, p, by_param[p])
            print(txt)
            tex.append(latex)
            any_rows = True

    if not any_rows:
        raise ValueError("no results found -- run the sweep first "
                         "(scripts/ablation.py; see the Parameter ablation section of README.md)")
    if latex:
        outfile.parent.mkdir(parents=True, exist_ok=True)
        outfile.write_text("\n".join(tex))
        print(f"\nwrote {outfile}")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--split", help=f"comma-separated subset of {', '.join(SPLIT_ORDER)}")
    ap.add_argument("--outfile", type=Path, default=None,
                    help="LaTeX output (default: <BENCH_DATA>/results/ablation-tables.generated.tex)")
    ap.add_argument("--no-latex", action="store_true")
    args = ap.parse_args()
    try:
        return run(args.split.split(",") if args.split else None, args.outfile,
                   latex=not args.no_latex)
    except ValueError as e:
        print(f"ablation_report: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
