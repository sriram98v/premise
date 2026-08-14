#!/usr/bin/env python3
"""Generate figures for the PREMISE parameter ablation.

Usage:
    python3 scripts/ablation_figures.py
    python3 scripts/ablation_figures.py --split syn-mix --only mem,rho
    python3 scripts/ablation_figures.py --formats pdf,png,svg
    python3 scripts/ablation_figures.py --no-relative
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import cnsplots as cns
import matplotlib
import pandas as pd
import utils
from matplotlib import font_manager as fm
from matplotlib.text import Text
from matplotlib.ticker import NullFormatter, ScalarFormatter

JOURNAL_FONT = "Arial"
TARGET_WIDTH_PT = 540.0
MIN_PT, MAX_PT = 8.0, 12.0
MATHTEXT_SHRINK = 0.7

_EXPECTED_PS = {
    ("normal", "normal"): "ArialMT",
    ("bold", "normal"): "Arial-BoldMT",
    ("normal", "italic"): "Arial-ItalicMT",
    ("bold", "italic"): "Arial-BoldItalicMT",
}


def assert_arial():
    """Fail loudly if matplotlib resolves anything but real Arial for any of the four faces."""
    for (weight, style), want in _EXPECTED_PS.items():
        fp = fm.FontProperties(family=JOURNAL_FONT, weight=weight, style=style)
        try:
            path = fm.findfont(fp, fallback_to_default=False)
        except ValueError as exc:
            raise SystemExit(
                f"Arial not installed ({weight}/{style}): {exc}\n"
                "Run inside `nix develop .#benchmark`, which supplies corefonts. If Arial was "
                "installed after matplotlib built its font cache, that cache is not invalidated "
                "by new fonts -- clear it with `rm -rf ~/.cache/matplotlib`."
            ) from exc
        got = fm.get_font(path).postscript_name
        if got != want:
            raise SystemExit(
                f"{weight}/{style} resolved to {got} ({path}), not {want} -- that is a metric "
                "clone (Liberation Sans / Arimo / Nimbus Sans), not Arial, and would violate the "
                "journal font policy."
            )


def apply_journal_style():
    """Pin fonts and sizes for journal compliance. Must run before the first cns.multipanel().
    """
    cns.settings.font_family = "sans-serif"
    cns.settings.font_sans_serif = (JOURNAL_FONT,)
    cns.settings.panel_label_fontname = JOURNAL_FONT
    cns.settings.mathtext_fontset = "custom"
    cns.settings.pdf_fonttype = 42
    cns.settings.title_fontsize = 12
    cns.settings.legend_fontsize = 9
    cns.settings.multipanel_title_height_pad = 22

    rc = matplotlib.rcParams
    rc["ps.fonttype"] = 42
    rc["axes.unicode_minus"] = False
    rc["mathtext.rm"] = rc["mathtext.sf"] = rc["mathtext.tt"] = JOURNAL_FONT
    rc["mathtext.it"] = f"{JOURNAL_FONT}:italic"
    rc["mathtext.bf"] = f"{JOURNAL_FONT}:bold"
    rc["mathtext.cal"] = f"{JOURNAL_FONT}:italic"
    if "mathtext.bfit" in rc:
        rc["mathtext.bfit"] = f"{JOURNAL_FONT}:italic:bold"
    rc["mathtext.fallback"] = None

    assert_arial()


def assert_printed_sizes(fig, stem: str):
    """Refuse to emit a figure that would break the width or 8-12 pt contract.

    TARGET_WIDTH_PT is a ceiling, not a target: a figure narrower than the print block is placed
    at its natural size, so its points are printed points. One wider than the block gets scaled
    down by the typesetter, and every size below shrinks with it -- which is the case to catch.
    """
    fig.canvas.draw()
    fig.canvas.draw()
    w_pt = fig.get_tightbbox(fig.canvas.get_renderer()).width * 72
    if w_pt > TARGET_WIDTH_PT:
        raise SystemExit(
            f"{stem}: tight width {w_pt:.1f} pt exceeds the {TARGET_WIDTH_PT:.0f} pt print block, "
            f"so placement would scale it by {TARGET_WIDTH_PT / w_pt:.3f} and drag text below "
            f"{MIN_PT:g} pt. Shrink PANEL_W/MAX_WIDTH or shorten the suptitle.")

    for t in fig.findobj(Text):
        s = t.get_text()
        if not t.get_visible() or not s.strip():
            continue
        eff = t.get_fontsize()
        smallest = eff * MATHTEXT_SHRINK if "$" in s else eff
        if smallest < MIN_PT or eff > MAX_PT:
            raise SystemExit(f"{stem}: {s!r} renders at {smallest:.2f}-{eff:.2f} pt printed, "
                             f"outside {MIN_PT:g}-{MAX_PT:g} pt")
        if t.axes is not None and any(ord(c) > 127 for c in s):
            raise SystemExit(f"{stem}: {s!r} is non-ASCII on an axes; cns.savefig()'s "
                             "apply_unicode_font() would switch it to DejaVu Sans")


ABLATION = Path(__file__).resolve().parents[1]


def figs() -> Path:
    return utils.data_root() / "ablation" / "figs"


SPLITS = {
    "syn-iso":  ("synthetic/isolate", "synthetic single-strain"),
    "syn-mix":  ("synthetic/mixed", "synthetic mixed-infection"),
    "real-iso": ("real/isolate", "real single-strain"),
    "real-mix": ("real/mixed", "real mixed-infection"),
}
SPLIT_ORDER = ["syn-iso", "syn-mix", "real-iso", "real-mix"]

PANELS = [
    ("wall_s",    "Runtime",   "Wall-clock time (s)",  "ratio"),
    ("rss_gb",    "Memory",    "Peak RSS (GB)",        "ratio"),
    ("precision", "Precision", "Precision",            "delta"),
    ("coverage",  "Coverage",  "Reads classified (%)", "delta"),
    ("ruzicka",   "Ruzicka",   "Ruzicka distance",     "ratio"),
    ("jaccard",   "Jaccard",   "Jaccard distance",     "delta"),
]
REL_YLABEL = {"ratio": "x baseline", "delta": "change vs baseline"}
RATIO_KEYS = {k for k, _, _, mode in PANELS if mode == "ratio"}

ZERO_VALUED = {"jaccard"}

PARAM_AXIS = {
    "mem":   ("Minimum seed length (-p)", "numeric"),
    "eps_2": ("$\\varepsilon_2$",         "categorical"),
    "eps_1": ("$\\varepsilon_1$",         "categorical"),
    "rho":   ("$\\rho$",                  "categorical"),
    "omega": ("$\\omega$",                "categorical"),
}
PARAM_TITLE = {
    "mem":   "Minimum seed length",
    "eps_2": "Alignment probability floor",
    "eps_1": "Pre-EM likelihood cutoff",
    "rho":   "$L_1$ penalty weight",
    "omega": "$L_1$ penalty floor",
}
PARAM_ORDER = ["mem", "eps_2", "eps_1", "rho", "omega"]

PANEL_W, PANEL_H, MAX_WIDTH = 115, 90, 540
PANELS_PER_ROW = 3
ROW_GAP = {"numeric": 34, "categorical": 52}


def load(split: str) -> pd.DataFrame | None:
    """One split's CSV as a DataFrame, with an x position per swept value.
    """
    csv_path = utils.data_root() / "results" / "ablation" / SPLITS[split][0] / "ablation.csv"
    if not csv_path.exists():
        return None
    df = pd.read_csv(csv_path, dtype={"value": str})
    if df.empty:
        return None
    df["baseline"] = df["baseline"].astype(bool)
    order = {}
    for param, sub in df.groupby("param", sort=False):
        order[param] = {v: i for i, v in enumerate(dict.fromkeys(sub["value"].tolist()))}
    df["pos"] = [order[p][v] for p, v in zip(df["param"], df["value"])]
    return df


def relativize(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise each metric to that (dataset, param) group's own baseline point.
    """
    out = []
    for _, grp in df.groupby(["dataset", "param"], sort=False):
        base = grp[grp["baseline"]]
        if base.empty:
            continue
        g = grp.copy()
        for key, _, _, mode in PANELS:
            b = float(base[key].iloc[0])
            g[key] = (g[key] / b if b else float("nan")) if mode == "ratio" else g[key] - b
        out.append(g)
    return pd.concat(out, ignore_index=True) if out else df.iloc[0:0]


def _style_axis(ax, df, key, xlabel, xkind, ylabel, title, relative):
    ax.set_title(title, fontsize=10, pad=3)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if xkind == "categorical":
        vals = df.drop_duplicates("pos").sort_values("pos")
        step = 1 if len(vals) <= 8 else 2
        ax.set_xticks(list(vals["pos"])[::step])
        ax.set_xticklabels(list(vals["value"])[::step], rotation=45, ha="right")

    if relative:
        ax.axhline(1.0 if key in RATIO_KEYS else 0.0, ls=":", lw=0.6, color="0.7", zorder=0)
        ax.ticklabel_format(axis="y", style="plain", useOffset=False)
    else:
        v = df[key].dropna()
        pos = v[v > 0]
        if (key not in ZERO_VALUED and len(pos) == len(v) and len(pos) > 1
                and pos.max() / pos.min() > 20):
            ax.set_yscale("log")
            ax.yaxis.set_major_formatter(ScalarFormatter())
            ax.yaxis.get_major_formatter().set_scientific(False)
            ax.yaxis.set_minor_formatter(NullFormatter())
        else:
            ax.ticklabel_format(axis="y", style="plain", useOffset=False)

    base = df[df["baseline"]]
    if len(base):
        x = float(base["pos"].iloc[0]) if xkind == "categorical" else float(base["value"].iloc[0])
        ax.axvline(x, ls=(0, (3, 2)), lw=0.6, color="0.55", zorder=0)


def figure_for(split: str, param: str, df: pd.DataFrame, formats, relative: bool):
    xlabel, xkind = PARAM_AXIS[param]
    xcol = "pos" if xkind == "categorical" else "value"
    if xkind == "numeric":
        df = df.copy()
        df["value"] = df["value"].astype(float)

    n_ds = df["dataset"].nunique()
    kind = "relative to production setting" if relative else "absolute"
    mp = cns.multipanel(
        max_width=MAX_WIDTH,
        title=f"PREMISE ablation ({SPLITS[split][1]}): {PARAM_TITLE[param]}\n"
              f"mean ± SD over {n_ds} samples, {kind}",
        loc="left")

    n_rows = -(-len(PANELS) // PANELS_PER_ROW)
    for i, (key, title, ylabel, mode) in enumerate(PANELS):
        row, col = divmod(i, PANELS_PER_ROW)
        label = "ABCDEF"[i]
        kw = {}
        if col == PANELS_PER_ROW - 1:
            kw["margin_right"] = 0
        if row < n_rows - 1:
            kw["margin_bottom"] = ROW_GAP[xkind]
        mp.panel(label, PANEL_W, PANEL_H, pad_top=4, **kw)
        cns.lineplot(data=df, x=xcol, y=key, marker="o", markersize=2.5, linewidth=1.0,
                     errorbar="sd", err_style="bars",
                     err_kws=dict(elinewidth=0.7, capsize=1.4, capthick=0.7))
        yl = f"{ylabel} ({REL_YLABEL[mode]})" if relative else ylabel
        _style_axis(mp.get_axes(label), df, key, xlabel, xkind, yl, title, relative)

    out_dir = figs() / split
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = f"ablation-{param.replace('_', '')}" + ("-rel" if relative else "")
    assert_printed_sizes(mp.fig, stem)
    written = []
    for ext in formats:
        p = out_dir / f"{stem}.{ext}"
        cns.savefig(p)
        written.append(p)
    return written


def run(splits_arg: list[str] | None = None, only: list[str] | None = None,
        formats: tuple[str, ...] = ("pdf", "png"), relative: bool = True) -> int:
    """Render the ablation figures under <BENCH_DATA>/ablation/figs."""

    apply_journal_style()

    formats = list(formats)
    splits = [s for s in SPLIT_ORDER if s in (splits_arg or SPLIT_ORDER)]
    wanted = only or PARAM_ORDER

    any_written = False
    for split in splits:
        df = load(split)
        if df is None:
            print(f"{split}: no results yet, skipped")
            continue
        rel_df = None if not relative else relativize(df)
        for param in PARAM_ORDER:
            if param not in wanted:
                continue
            sub = df[df["param"] == param]
            if sub.empty:
                continue
            for p in figure_for(split, param, sub, formats, relative=False):
                print(f"wrote {p.relative_to(utils.data_root())}  "
                      f"({sub['pos'].nunique()} points x {sub['dataset'].nunique()} samples)")
                any_written = True
            if rel_df is not None:
                rsub = rel_df[rel_df["param"] == param]
                if not rsub.empty:
                    for p in figure_for(split, param, rsub, formats, relative=True):
                        print(f"wrote {p.relative_to(utils.data_root())}")
    if not any_written:
        raise ValueError("no figures written -- run the sweep first "
                         "(scripts/ablation.py; see the Parameter ablation section of README.md)")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--split", help=f"comma-separated subset of {', '.join(SPLIT_ORDER)}")
    ap.add_argument("--only", help="comma-separated subset of parameters to plot")
    ap.add_argument("--formats", default="pdf,png",
                    help="comma-separated output formats (pdf,png,svg)")
    ap.add_argument("--no-relative", action="store_true",
                    help="skip the baseline-normalised companion figures")
    args = ap.parse_args()
    try:
        return run(args.split.split(",") if args.split else None,
                   args.only.split(",") if args.only else None,
                   tuple(f.strip() for f in args.formats.split(",") if f.strip()),
                   relative=not args.no_relative)
    except ValueError as e:
        print(f"ablation_figures: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
