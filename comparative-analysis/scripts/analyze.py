#!/usr/bin/env python3
"""Comparative analysis for the PREMISE benchmark.

Usage:
    python3 analyze.py [--splits real-iso,syn-iso,syn-mix] [--datasets D1,D2]
                             [--outfile results/tables.generated.tex] [--no-latex]
"""
from __future__ import annotations

import argparse
import gzip
import sys
from collections import Counter
from pathlib import Path

from bench_metrics import precision_recall, profile_distances
from utils import data_root, parse_timemem, results, strain_of, strip_version

SYNTH_DS = ["Dataset-1", "Dataset-2", "Dataset-3", "Dataset-4"]
REAL_DS = ["SRR31013463", "SRR31013465", "SRR31013467", "SRR31013473"]
REALMIXED_DS = ["SRR3360139", "SRR3360140", "SRR3360145", "SRR3360146"]
SPLITS = {
    "syn-iso": {"sub": "synthetic/isolate", "datasets": SYNTH_DS, "real": False},
    "syn-mix": {"sub": "synthetic/mixed", "datasets": SYNTH_DS, "real": False},
    "syn-mix-subtype": {"sub": "synthetic/mixed-subtype", "datasets": SYNTH_DS, "real": False},
    "real-iso": {"sub": "real/isolate", "datasets": REAL_DS, "real": True},
    "real-mix": {"sub": "real/mixed", "datasets": REALMIXED_DS, "real": True},
}

SPLIT_CAPTION = {
    "syn-iso": "synthetic",
    "syn-mix": "mixed-infection synthetic",
    "syn-mix-subtype": "same-subtype mixed-infection synthetic",
    "real-iso": "real",
    "real-mix": "real mixed-infection (influenza)",
}

METHODS = [
    dict(code="pre", name="Premise",     dirn="premise",     profiler=True,  perread=True,  has_uc=True,  add_other=True),
    dict(code="cen", name="Centrifuger", dirn="centrifuger", profiler=True,  perread=True,  has_uc=True,  add_other=True),
    dict(code="kmp", name="KMCP",        dirn="kmcp",        profiler=True,  perread=True,  has_uc=False, add_other=True),
    dict(code="mor", name="MORA",        dirn="mora",        profiler=False, perread=True,  has_uc=True,  add_other=False),
    dict(code="kap", name="Karp",        dirn="karp",        profiler=True,  perread=False, has_uc=False, add_other=True),
    dict(code="gan", name="Ganon",       dirn="ganon",       profiler=False, perread=True,  has_uc=True,  add_other=False),
    dict(code="syl", name="Sylph",       dirn="sylph",       profiler=True,  perread=False, has_uc=False, add_other=True),
]
METHOD_BY_CODE = {m["code"]: m for m in METHODS}

PREC_REC_EXCLUDE = {"syl", "kap"}


def _warn(msg: str):
    print(f"analyze.py: WARNING: {msg}", file=sys.stderr)


def _safe_lines(f, path, state):
    """Yield lines from a (possibly gzip) handle, stopping at a truncated stream.
    """
    try:
        for line in f:
            yield line
    except (EOFError, OSError) as e:
        state["truncated"] = True
        _warn(f"{path}: truncated/corrupt stream ({type(e).__name__}: {e})")


def _open_maybe_gz(path: Path):
    """Open `path`, falling back to `path`.gz — same format, different container.

    Missing is fatal, not empty: an absent truth file used to yield an empty dict and a
    fully rendered table of NaNs, which is indistinguishable from a real result.
    """
    if path.exists():
        return open(path)
    gz = path.with_suffix(path.suffix + ".gz")
    if gz.exists():
        return gzip.open(gz, "rt")
    raise ValueError(
        f"{path} not found — build it with:\n"
        f"  python3 scripts/prepare_real_samples.py --samples {path.parent.name}")


def _suffix(d: Path, code: str, split: str, ds: str) -> str:
    """'' for the plain <dataset> base the drivers now write, '.ca' for older trees.
    """
    if split != "real-iso" or code not in ("pre", "kmp"):
        return ""
    if any(not p.name.startswith(f"{ds}.ca.") for p in d.glob(f"{ds}.*")):
        return ""
    return ".ca" if any(d.glob(f"{ds}.ca.*")) else ""


def _norm_abund(counts: dict) -> dict:
    total = sum(counts.values())
    return {k: v / total for k, v in counts.items()} if total else {}


def load_truth(split: str, ds: str):
    """Return (truth_reads, truth_counts, n_reads). truth_reads maps readID -> stripped ref.
    """
    info = SPLITS[split]
    if info["real"]:
        sd = data_root() / "samples" / info["sub"] / ds / "truth_assignments.tsv"
        truth_seg = {}
        with _open_maybe_gz(sd) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2 and p[0] and p[1]:
                    truth_seg[p[0]] = strip_version(p[1])
        counts = Counter(truth_seg.values())
        if split != "real-mix":
            return truth_seg, dict(counts), len(truth_seg)
        truth = {q: strain_of(r) for q, r in truth_seg.items()}
        bad = set(truth.values()) - {"PR8", "WSN33"}
        if bad:
            raise ValueError(
                f"{sd}: {sorted(bad)} — every reference must be a PR8_/WSN33_ "
                "segment on real-mix, or the strain truth collapses to 'other'")
        return truth, dict(counts), len(truth)
    base = results() / "premise" / info["sub"] / ds / ds
    truth = {}
    with open(f"{base}.matches") as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 1 and p[0]:
                truth[p[0]] = strip_version(p[0].split("_")[0])
    return truth, dict(Counter(truth.values())), len(truth)


_R1_PATH = {
    "syn-iso":         lambda ds: data_root() / "samples" / "synthetic" / "isolate" / ds / "reads_R1.fastq",
    "syn-mix":         lambda ds: data_root() / "samples" / "synthetic" / "mixed"   / ds / "reads_R1.fastq",
    "syn-mix-subtype": lambda ds: data_root() / "samples" / "synthetic" / "mixed-subtype" / ds / "reads_R1.fastq",
    "real-iso":        lambda ds: data_root() / "samples" / "real" / "isolate" / ds / f"{ds}_1-filtered.ca.fastq",
    "real-mix":        lambda ds: data_root() / "samples" / "real" / "mixed"   / ds / f"{ds}_1-filtered.ca.fastq",
}


def input_read_count(split: str, ds: str):
    """Total read (pair) count in the dataset's R1 input fastq, or None if unavailable."""
    p = _R1_PATH[split](ds)
    opener = (lambda: gzip.open(p, "rt")) if p.suffix == ".gz" else (lambda: open(p))
    try:
        if not p.exists() or p.stat().st_size == 0:
            return None
        with opener() as f:
            return sum(1 for _ in f) // 4
    except (OSError, EOFError):
        return None


def _read_2col_tsv(path: Path, strip_mate=False):
    """readID<TAB>seqID with header -> {readID: stripped ref | None}."""
    out = {}
    if not path.exists() or path.stat().st_size == 0:
        return None
    with open(path) as f:
        next(f, None)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2 or not p[0]:
                continue
            rid = p[0].split("/", 1)[0] if strip_mate else p[0]
            ref = p[1].strip()
            out[rid] = strip_version(ref) if ref and ref != "unclassified" else None
    return out


def load_premise(d: Path, base: str):
    matches = d / f"{base}.matches"
    assign = None
    if matches.exists() and matches.stat().st_size:
        assign = {}
        with open(matches) as f:
            next(f)
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < 2:
                    continue
                assign[p[0]] = None if p[1] == "unclassified" else strip_version(p[1])
    props = d / f"{base}.props"
    abund = None
    if props.exists() and props.stat().st_size:
        c = Counter()
        with open(props) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2:
                    try:
                        c[strip_version(p[0])] += float(p[1])
                    except ValueError:
                        pass
        abund = _norm_abund(c)
    return assign, abund


def load_centrifuger(d: Path, base: str):
    tsv = d / f"{base}.tsv"
    if not tsv.exists() or not tsv.stat().st_size:
        return None, None
    assign = {}
    counts = Counter()
    with open(tsv) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2:
                continue
            ref = strip_version(p[1])
            if ref in ("no rank", "species"):
                ref = None
            assign[p[0]] = ref
            if ref:
                counts[ref] += 1
    return assign, _norm_abund(counts)


def load_kmcp(d: Path, base: str):
    prof = d / f"{base}.profile"
    abund = None
    if prof.exists() and prof.stat().st_size:
        c = Counter()
        with open(prof) as f:
            hdr = next(f).rstrip("\n").split("\t")
            try:
                ri, pi = hdr.index("ref"), hdr.index("percentage")
            except ValueError:
                ri = pi = None
            if ri is not None:
                for line in f:
                    p = line.rstrip("\n").split("\t")
                    if len(p) <= max(ri, pi) or p[0].startswith("#"):
                        continue
                    ref = strip_version(p[ri].replace("sequences.id_", ""))
                    try:
                        c[ref] += float(p[pi])
                    except ValueError:
                        pass
        abund = _norm_abund(c)
    gz = d / f"{base}.tsv.gz"
    plain = d / f"{base}.tsv"
    tsv = plain if (plain.exists() and plain.stat().st_size) else gz
    assign = None
    if tsv.exists() and tsv.stat().st_size:
        opener = gzip.open if tsv.suffix == ".gz" else open
        state = {"truncated": False}
        with opener(tsv, "rt") as f:
            best = {}
            try:
                first = f.readline()
            except (EOFError, OSError) as e:
                first = ""
                state["truncated"] = True
                _warn(f"{tsv}: truncated/corrupt stream ({type(e).__name__}: {e})")
            hdr = first.rstrip("\n").split("\t") if first else []
            try:
                qi = hdr.index("qCov")
            except ValueError:
                qi = None
            for line in _safe_lines(f, tsv, state):
                if line.startswith("#"):
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 6:
                    continue
                rid = p[0].split("/", 1)[0]
                ref = strip_version(p[5].replace("sequences.id_", "")) if len(p) > 5 else None
                try:
                    q = float(p[qi]) if qi is not None else 0.0
                except ValueError:
                    q = 0.0
                if rid not in best or q > best[rid][0]:
                    best[rid] = (q, ref)
        if state["truncated"]:
            _warn(f"kmcp result discarded for {tsv.parent.name}: search did not finish")
            assign = None
            abund = None
        else:
            assign = {r: v[1] for r, v in best.items()}
    return assign, abund


def load_sylph(d: Path, base: str):
    props = d / f"{base}.norm.tsv.props"
    if not props.exists():
        props = d / f"{base}.props"
    abund = None
    if props.exists() and props.stat().st_size:
        c = Counter()
        with open(props) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2:
                    try:
                        c[strip_version(p[0])] += float(p[1])
                    except ValueError:
                        pass
        abund = _norm_abund(c)
    return None, abund


def load_karp(d: Path, base: str):
    """Karp EM abundance from <base>.freqs (Label<TAB>ExpectedCounts<TAB>Taxa).
    """
    freqs = d / f"{base}.freqs"
    abund = None
    if freqs.exists() and freqs.stat().st_size:
        c = Counter()
        with open(freqs) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2:
                    try:
                        c[strip_version(p[0])] += float(p[1])
                    except ValueError:
                        pass
        abund = _norm_abund(c)
    return None, abund


def load_assignment_method(d: Path, base: str, strip_mate=False):
    """mora/ganon: normalized 2-col .tsv; abundance from assignment counts."""
    assign = _read_2col_tsv(d / f"{base}.tsv", strip_mate=strip_mate)
    if assign is None:
        return None, None
    abund = _norm_abund(Counter(r for r in assign.values() if r is not None))
    return assign, abund


def load_method(code: str, split: str, ds: str):
    info = SPLITS[split]
    d = results() / METHOD_BY_CODE[code]["dirn"] / info["sub"] / ds
    base = ds + _suffix(d, code, split, ds)
    if code == "pre":
        return load_premise(d, base)
    if code == "cen":
        return load_centrifuger(d, base)
    if code == "kmp":
        return load_kmcp(d, base)
    if code == "syl":
        return load_sylph(d, base)
    if code == "kap":
        return load_karp(d, base)
    if code in ("mor", "gan"):
        return load_assignment_method(d, base)
    return None, None


def compute_dataset(split: str, ds: str):
    info = SPLITS[split]
    truth_reads, truth_counts, n_reads = load_truth(split, ds)
    rows = {}
    for m in METHODS:
        assign, abund = load_method(m["code"], split, ds)
        if split == "real-mix" and assign is not None:
            assign = {r: strain_of(v) for r, v in assign.items()}
        cr = None
        if assign is not None:
            classified = sum(1 for v in assign.values() if v is not None)
            denom = input_read_count(split, ds)
            if denom:
                cr = 100.0 * classified / denom
        dist = profile_distances(truth_counts, n_reads, abund,
                                 uncl_frac=0.0, has_uc=m["has_uc"], add_other=m["add_other"])
        pr = precision_recall(truth_reads, assign if m["perread"] else None, synthetic=not info["real"])
        rows[m["code"]] = dict(dist=dist, pr=pr, cr=cr,
                               has_data=(abund is not None or assign is not None))
    return rows


def abund_cell(row):
    return row["dist"]["ruz.uc"], row["dist"]["jac.uc"]


def pr_cell(row, real):
    """(precision, coverage). Coverage is the percentage of input reads assigned a label —
    the metric previously shown as 'Tagged' in the runtime-memory tables. Recall is no longer
    reported: on the real splits it is largely the complement of coverage (an unassigned read
    is a miss), so the two columns restated the same quantity."""
    pr = row["pr"]
    prec = pr["prec_uc"] if real else pr["prec"]
    return (prec, row["cr"])


def _isnum(x):
    return x is not None


def _numstr(x, sig=None):
    """Format a number: `sig` significant digits (%g) when given, else 2 decimal places."""
    return f"{x:.{sig}g}" if sig else f"{x:.2f}"


def fmt(x, dash="-", sig=None):
    return _numstr(x, sig) if _isnum(x) else dash


def _percol(spec, k):
    """Broadcast a scalar to a length-k list, or pass a list/tuple through unchanged."""
    return list(spec) if isinstance(spec, (list, tuple)) else [spec] * k


def _best_per_column(matrix, lower_better, k, sig=None):
    """matrix: list of rows, each a flat list of numeric-or-None cells (k sub-columns per dataset,
    repeating). lower_better is a bool or a length-k list giving the 'better' direction per
    sub-column. Return the set of (row,col) holding the best finite value in their column.
    `sig` mirrors fmt's significant-digit mode so tie-bolding matches the displayed precision;
    it may also be a length-k list when sub-columns are formatted differently."""
    best = set()
    if not matrix:
        return best
    lb = _percol(lower_better, k)
    sg = _percol(sig, k)
    ncol = len(matrix[0])
    for c in range(ncol):
        vals = [(r, matrix[r][c]) for r in range(len(matrix)) if _isnum(matrix[r][c])]
        if not vals:
            continue
        lower = lb[c % k]
        target = min(v for _, v in vals) if lower else max(v for _, v in vals)
        tstr = _numstr(target, sg[c % k])
        for r, v in vals:
            if _numstr(v, sg[c % k]) == tstr:
                best.add((r, c))
    return best


def render(split_title, datasets, method_rows, col_labels, lower_better, dash,
           caption=None, label=None, sig=None):
    """method_rows: list of (tex_name, plain_name, [cell per dataset]) where each cell is a
    k-tuple aligned to col_labels (k sub-columns). lower_better is a bool or length-k list.
    `sig`: significant digits (%g); None keeps the default 2-decimal formatting. May be a
    length-k list to format sub-columns differently (e.g. a 4-sig-digit precision beside a
    2-decimal percentage). Returns (stdout, latex)."""
    k = len(col_labels)
    sg = _percol(sig, k)
    matrix = [[v for cell in cells for v in cell] for _, _, cells in method_rows]
    best = _best_per_column(matrix, lower_better, k, sig)
    nds = len(datasets)

    out = [f"\n=== {split_title} ({'/'.join(col_labels)}) ==="]
    head = f"{'Method':13}" + "".join(f"{ds[:9 * k]:>{9 * k}}" for ds in datasets)
    out.append(head)
    sub = f"{'':13}" + "".join("".join(f"{cl[:8]:>9}" for cl in col_labels) for _ in datasets)
    out.append(sub)
    for _, name, cells in method_rows:
        line = f"{name:13}"
        for cell in cells:
            for j, v in enumerate(cell):
                line += f"{fmt(v, dash, sg[j]):>9}"
        out.append(line)

    colspec = "l" + "c" * (k * nds)
    cap = caption if caption else split_title
    lab = label if label else _slug(split_title)
    L = ["\\begin{table*}", f"\\caption{{{cap}.\\label{{tab:{lab}}}}}",
         "\\tabcolsep=1pt%%",
         f"\\begin{{tabular*}}{{\\textwidth}}{{@{{\\extracolsep{{\\fill}}}}{colspec}@{{\\extracolsep{{\\fill}}}}}}",
         "\\toprule%"]
    hdr = "".join(f" & \\multicolumn{{{k}}}{{@{{}}c@{{}}}}{{\\textbf{{{ds}}}}}" for ds in datasets)
    L.append(hdr + "\\\\")
    L.append("".join(f"\\cline{{{2 + k * i}-{1 + k + k * i}}}" for i in range(nds)) + "%")
    L.append("\\textbf{Method}" + "".join("".join(f" & \\textbf{{{cl}}}" for cl in col_labels)
                                          for _ in datasets) + " \\\\ \\midrule")
    for ri, (tex_name, _, cells) in enumerate(method_rows):
        vals = [v for cell in cells for v in cell]
        cellstr = ""
        for ci, v in enumerate(vals):
            s = fmt(v, dash, sg[ci % k])
            if (ri, ci) in best and _isnum(v):
                s = f"\\textbf{{{s}}}"
            cellstr += f" & {s}"
        L.append(f"{tex_name}{cellstr} \\\\")
    L += ["\\bottomrule", "\\end{tabular*}", "\\end{table*}", ""]
    return "\n".join(out), "\n".join(L)


def _slug(title):
    return title.lower().replace(" ", "-").replace("(", "").replace(")", "").replace("/", "-")


def _display_datasets(split, dss):
    """Column labels for a split. Synthetic splits use their physical Dataset-1..N names directly
    (no dropping/renumbering); real splits keep their accession names."""
    return list(dss)


def resource_rows(split, datasets):
    """Per-method (wall_s, rss_gb) per dataset, from the time-mem files. The classified rate
    moved to the precision-coverage tables as 'Coverage'."""
    info = SPLITS[split]
    rows = []
    for m in METHODS:
        cells = []
        for ds in datasets:
            d = results() / m["dirn"] / info["sub"] / ds
            tm = d / "time-mem"
            if not tm.exists() and (d / "time-mem.ca").exists():
                tm = d / "time-mem.ca"
            cells.append(parse_timemem(tm))
        tex = "\\premise" if m["code"] == "pre" else m["name"]
        rows.append((tex, m["name"], cells))
    return rows


def database_table():
    """DB size (GB) + build time (s) from db_build.csv. Returns (stdout, latex) or (None, None)."""
    csv = data_root() / "db_build.csv"
    if not csv.exists():
        return None, None
    by = {}
    for line in csv.read_text().splitlines()[1:]:
        p = line.split(",")
        if len(p) >= 4:
            by[p[0]] = (p[2], p[3])
    out = ["\n=== Database construction ===", f"{'Method':13}{'Size(GB)':>10}{'Time(s)':>10}"]
    lat = ["\\begin{table}",
           "\\caption{Comparative analysis of custom database construction. Each database is "
           "built from the same reference collection; size is the final index required for "
           "classification (excluding regenerable intermediates). Size is reported in GB "
           "($10^9$ bytes) and build time in seconds."
           "\\label{tab:database-comp}}",
           "\\centering", "\\begin{tabular}{lcc}", "\\toprule",
           "\\textbf{Method} & \\textbf{Size} & \\textbf{Time} \\\\ \\midrule"]
    for m in METHODS:
        if m["dirn"] not in by:
            continue
        size_after, secs = by[m["dirn"]]
        gb = f"{int(size_after) / 1e9:.4g}" if size_after.isdigit() else "-"
        try:
            tsec = f"{float(secs):.4g}"
        except ValueError:
            tsec = "-"
        out.append(f"{m['name']:13}{gb:>10}{tsec:>10}")
        tex = "\\premise" if m["code"] == "pre" else m["name"]
        lat.append(f"{tex} & {gb} & {tsec} \\\\")
    lat += ["\\bottomrule", "\\end{tabular}", "\\end{table}", ""]
    return "\n".join(out), "\n".join(lat)


def build_tables(splits, datasets_filter):
    stdout_blocks, latex_blocks = [], []
    dbo, dbl = database_table()
    if dbo:
        stdout_blocks.append(dbo); latex_blocks.append(dbl)
    for split in splits:
        info = SPLITS[split]
        dss = [d for d in info["datasets"] if not datasets_filter or d in datasets_filter]
        if not dss:
            continue
        per_ds = {}
        for ds in dss:
            try:
                per_ds[ds] = compute_dataset(split, ds)
            except (FileNotFoundError, OSError) as e:
                print(f"  [skip {split}/{ds}: {e}]")
        dss = [d for d in dss if d in per_ds]
        if not dss:
            continue
        # The \label tag IS the split code — tab:class-comp-real-mix-time, etc. No separate
        # tag table: a second vocabulary for the same five things is what this rename removed.
        cap, tag = SPLIT_CAPTION[split], split
        disp = _display_datasets(split, dss)

        rr = resource_rows(split, dss)
        so, la = render(f"{split.capitalize()} runtime-memory", disp, rr,
                        ("Time", "Mem"), [True, True], "-",
                        caption=f"Comparative analysis of runtime and peak memory usage on {cap} "
                                f"datasets. Time is reported in seconds and peak memory in GB",
                        label=f"class-comp-{tag}-time")
        stdout_blocks.append(so); latex_blocks.append(la)

        pr_rows = []
        for m in METHODS:
            if m["code"] in PREC_REC_EXCLUDE:
                continue
            tex = "\\premise" if m["code"] == "pre" else m["name"]
            if m["perread"]:
                pr_rows.append((tex, m["name"], [pr_cell(per_ds[ds][m["code"]], info["real"]) for ds in dss]))
            else:
                pr_rows.append((tex, m["name"], [(None, None) for _ in dss]))
        so, la = render(f"{split.capitalize()} precision-coverage", disp, pr_rows,
                        ("Precision", "Coverage"), False, "---",
                        caption=f"Comparative analysis of precision and coverage on {cap} "
                                f"datasets. Coverage is the percentage of reads assigned a label",
                        label=f"class-comp-{tag}-prec-rec", sig=[4, None])
        stdout_blocks.append(so); latex_blocks.append(la)

        ab_rows = []
        for m in METHODS:
            tex = "\\premise" if m["code"] == "pre" else m["name"]
            ab_rows.append((tex, m["name"], [abund_cell(per_ds[ds][m["code"]]) for ds in dss]))
        so, la = render(f"{split.capitalize()} abundance", disp, ab_rows,
                        ("Ruzicka", "Jaccard"), True, "-",
                        caption=f"Comparative analysis of source prediction and abundance estimation on {cap} datasets",
                        label=f"class-comp-{tag}-perf-abun", sig=4)
        stdout_blocks.append(so); latex_blocks.append(la)
    return "\n".join(stdout_blocks), "\n".join(latex_blocks)


def run(splits, datasets=None, outfile=None, latex: bool = True) -> str:
    """Score `splits` and return the report text; write the LaTeX tables when `latex`.

    Returns the text rather than printing it, so an importing caller controls where it goes and
    gets exactly what the CLI would have printed. Unknown split codes are dropped, matching what
    --splits has always done.
    """
    splits = [s for s in splits if s in SPLITS]
    dsf = set(datasets) or None if datasets else None
    stdout_text, latex_text = build_tables(splits, dsf)
    if latex:
        outfile = outfile or results() / "tables.generated.tex"
        with open(outfile, "w") as f:
            f.write("% Generated by analyze.py — do not edit by hand.\n")
            f.write(latex_text)
        stdout_text += f"\n\n[wrote LaTeX tables to {outfile}]"
    return stdout_text


def main() -> int:
    ap = argparse.ArgumentParser(description="Comparative analysis -> stdout + LaTeX tables.")
    ap.add_argument("--splits", default="syn-iso,syn-mix,syn-mix-subtype,real-iso,real-mix",
                    help=f"comma-separated subset of {', '.join(SPLITS)}")
    ap.add_argument("--datasets", default="")
    ap.add_argument("--outfile", default=None,
                    help="LaTeX output path (default: <BENCH_DATA>/results/tables.generated.tex)")
    ap.add_argument("--no-latex", action="store_true")
    args = ap.parse_args()
    try:
        print(run(args.splits.split(","),
                  [d for d in args.datasets.split(",") if d],
                  args.outfile, latex=not args.no_latex))
    except (ValueError, OSError) as e:
        print(f"analyze.py: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
