#!/usr/bin/env python3
"""Emit the LaTeX source listings for the synthetic datasets.
"""
from __future__ import annotations

from pathlib import Path

import argparse
import re
import sys
from collections import OrderedDict

import utils


# Resolved per call, not at import: a constant would bind to whatever $BENCH_DATA held when this
# module was first imported.
def syn() -> "Path":
    return utils.samples() / "synthetic"

STRAIN_RE = re.compile(r"Influenza A virus \((.+?)\) segment")
SEG_RE = re.compile(r"segment (\d)")
SUBTYPE_RE = re.compile(r"\(H(\d+)N(\d+)\)")

HA_SEGMENT = 4
NA_SEGMENT = 6

TEX_ESCAPE = str.maketrans({"&": r"\&", "%": r"\%", "_": r"\_", "#": r"\#"})


def tex(s: str) -> str:
    return s.translate(TEX_ESCAPE)


def parse_src(d: Path):
    """Return [(accession, segment, strain)] in file order."""
    out = []
    for line in (d / "src.fasta").read_text().splitlines():
        if not line.startswith(">"):
            continue
        head = line[1:]
        acc = head.split()[0]
        strain = STRAIN_RE.search(head)
        seg = SEG_RE.search(head)
        out.append((acc, int(seg.group(1)), strain.group(1).replace(" (", "(")))
    return out


def parse_abundance(d: Path) -> dict[str, float]:
    out = {}
    for line in (d / "src-abundance.txt").read_text().splitlines():
        if line.strip():
            acc, prop = line.split()
            out[acc] = float(prop)
    return out


def genomes(d: Path):
    """Group records into genomes by shared proportion, preserving file order.

    Returns [(genome_share, [(seg, acc, strain), ...])] with segments sorted 1..8.
    """
    recs = parse_src(d)
    abund = parse_abundance(d)
    groups: OrderedDict[float, list] = OrderedDict()
    for acc, seg, strain in recs:
        groups.setdefault(abund[acc], []).append((seg, acc, strain))
    return [(p * len(g), sorted(g)) for p, g in groups.items()]


def subtype(segs) -> str:
    """The chimera's own subtype: H from the HA-segment donor, N from the NA-segment donor.

    Generally differs from either donor's subtype, since the two segments come from
    different strains.
    """
    by_seg = {s: strain for s, _, strain in segs}
    h = SUBTYPE_RE.search(by_seg[HA_SEGMENT]).group(1)
    n = SUBTYPE_RE.search(by_seg[NA_SEGMENT]).group(2)
    return f"H{h}N{n}"


def details_rows():
    """The subtype summary strings for the Details column of tab:datasets."""
    iso = [subtype(genomes(d)[0][1]) for d in sorted((syn() / "isolate").glob("Dataset-*"))]
    mixed = []
    for d in sorted((syn() / "mixed").glob("Dataset-*")):
        mixed.append((d.name, [subtype(segs) for _, segs in genomes(d)]))
    return iso, mixed


def emit(split: str, label: str, caption: str, multi: bool):
    """Render one source table. ``multi`` adds the per-genome column, which is
    only meaningful for the mixtures (the single-genome sets would repeat 1.000)."""
    ncol = 5 if multi else 4
    print("\\begin{table*}[t]")
    print(f"\\caption{{{caption}\\label{{{label}}}}}")
    print("\\centering")
    print("\\footnotesize")
    print("\\begin{tabular}{%s}" % ("cclll" if multi else "clll"))
    print("\\toprule")
    head = ["\\textbf{Dataset}"]
    if multi:
        head.append("\\textbf{Genome (share)}")
    head += ["\\textbf{Seg.}", "\\textbf{Accession}", "\\textbf{Donor strain}"]
    print(" & ".join(head) + " \\\\")
    print("\\midrule")

    dirs = sorted((syn() / split).glob("Dataset-*"), key=lambda p: p.name)
    for di, d in enumerate(dirs):
        gs = genomes(d)
        nrow = sum(len(g) for _, g in gs)
        first_of_dataset = True
        for gi, (share, segs) in enumerate(gs):
            gname = chr(ord("A") + gi)
            for si, (seg, acc, strain) in enumerate(segs):
                cells = []
                if first_of_dataset:
                    tag = d.name if multi else f"{d.name} ({subtype(segs)})"
                    cells.append(f"\\multirow{{{nrow}}}{{*}}{{{tag}}}")
                    first_of_dataset = False
                else:
                    cells.append("")
                if multi:
                    cells.append(
                        f"\\multirow{{{len(segs)}}}{{*}}"
                        f"{{{gname}, {subtype(segs)} ({share:.3f})}}"
                        if si == 0 else ""
                    )
                cells += [str(seg), tex(acc), tex(strain)]
                print(" & ".join(cells) + " \\\\")
            if gi < len(gs) - 1:
                print(f"\\cmidrule(l){{2-{ncol}}}")
        if di < len(dirs) - 1:
            print("\\midrule")
    print("\\bottomrule")
    print("\\end{tabular}")
    print("\\end{table*}")


def run() -> None:
    """Print the Details column summary and both source tables to stdout."""
    iso, mixed = details_rows()
    print("%% Details column for tab:datasets, Set 1:")
    print("%%   " + ", ".join(iso))
    print("%% Details column for tab:datasets, Set 2:")
    print("%%   " + "; ".join(f"{n} ({', '.join(s)})" for n, s in mixed))
    print()
    emit(
        "isolate",
        "tab:syn-isolate-sources",
        "Source sequences for the synthetic single-genome datasets (Set 1). Each "
        "dataset is one reassortant-like chimeric genome of eight influenza A "
        "segments, every segment drawn from a different strain, at uniform "
        "segment abundance (0.125 each).",
        multi=False,
    )
    print()
    emit(
        "mixed",
        "tab:syn-mixed-sources",
        "Source sequences for the synthetic mixture datasets (Set 2). Each dataset "
        "contains three reassortant-like chimeric genomes of eight influenza A "
        "segments each; the parenthesised value is the genome's share of the "
        "sample, split uniformly across its eight segments.",
        multi=True,
    )


def main() -> int:
    argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter).parse_args()
    try:
        run()
    except (ValueError, OSError) as e:
        print(f"make_datasets_tex: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
