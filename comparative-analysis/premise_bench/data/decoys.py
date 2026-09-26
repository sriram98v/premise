"""Near-twin decoy references for the real splits, added to the shared reference set.
"""
from __future__ import annotations

import os
import random
import re
import statistics
from dataclasses import dataclass
from pathlib import Path

from ..utils import SEG_RE, die, iter_fasta
from ..config import CLEANED_FASTA as CLEANED, data_root, indexes, section
from ..splits import SPLITS

PREFIX = "DECOY_"
GENE_SEGMENT = (("(PB2)", 1), ("(PB1)", 2), ("(PA)", 3), ("(HA)", 4), ("(NP)", 5),
                ("(NA)", 6), ("(M1)", 7), ("(M2)", 7), ("(NS1)", 8), ("(NS2)", 8), ("(NEP)", 8))
STRAIN_SPLIT = "syn-mix-strain"
REAL_SPLITS = ("real-iso", "real-mix")
TAXMAP = "seqid2taxid.map"
TOLERANCE = 0.01
MAX_TRIES = 25
TRANSITION = {"A": "G", "G": "A", "C": "T", "T": "C"}
BASES = "ACGT"


@dataclass(frozen=True)
class Params:
    seed: int
    per_source: int
    ti_frac: float
    k: int


@dataclass(frozen=True)
class Source:
    acc: str
    header: str
    seq: str
    segment: int
    split: str
    samples: tuple[str, ...]


@dataclass(frozen=True)
class Decoy:
    acc: str
    source: Source
    seq: str
    target_c: float
    achieved_c: float
    jaccard21: float
    m: int
    shared: int
    sibling_c: float
    tries: int

def read_fasta(path: Path) -> dict[str, tuple[str, str]]:
    """{accession: (full header without '>', uppercased sequence)} in file order."""
    return {acc: (header, seq) for acc, header, seq in iter_fasta(path)}


def segment_of(header: str) -> int | None:
    m = SEG_RE.search(header)
    if m:
        return int(m.group(1))
    return next((seg for gene, seg in GENE_SEGMENT if gene in header), None)


def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def kmers(seq: str, k: int) -> set[str]:
    """Canonical k-mers (min of the k-mer and its reverse complement)."""
    rc = revcomp(seq)
    n = len(seq)
    return {min(seq[i:i + k], rc[n - i - k:n - i]) for i in range(n - k + 1)}


def containment(a: set[str], b: set[str]) -> float:
    return len(a & b) / min(len(a), len(b))


def jaccard(a: set[str], b: set[str]) -> float:
    return len(a & b) / len(a | b)


def params_from_toml() -> Params:
    """[decoy] of params.toml; every value is a quoted string by design."""
    sec = section("decoy")
    try:
        p = Params(seed=int(sec["seed"]), per_source=int(sec["per_source"]),
                   ti_frac=float(sec["ti_frac"]), k=int(sec["k"]))
    except KeyError as e:
        raise SystemExit(f"params.toml: [decoy] has no {e} "
                         f"(have: {', '.join(sorted(sec)) or 'none'})") from None
    except ValueError as e:
        raise SystemExit(f"params.toml: [decoy] value not numeric: {e}") from None
    if p.per_source < 0 or not 0.0 <= p.ti_frac <= 1.0 or p.k < 11:
        die(1, f"params.toml: bad [decoy] parameters: {p}")
    return p

def target_pool(root: Path, db: dict[str, tuple[str, str]], k: int) -> dict[int, list[tuple]]:
    """Per segment: [(containment, dataset, acc_a, acc_b)] over every syn-mix-strain triple."""
    pool: dict[int, list[tuple]] = {s: [] for s in range(1, 9)}
    for ds in sorted((root / SPLITS[STRAIN_SPLIT].sample_dir).glob("Dataset-*")):
        tsv = ds / "genomes.tsv"
        if not tsv.exists():
            die(1, f"{tsv} missing — regenerate {STRAIN_SPLIT} first")
        by_seg: dict[int, list[str]] = {}
        with open(tsv) as f:
            head = f.readline().rstrip("\n").split("\t")
            si, ai = head.index("segment"), head.index("accession")
            for line in f:
                p = line.rstrip("\n").split("\t")
                by_seg.setdefault(int(p[si]), []).append(p[ai])
        for seg, accs in sorted(by_seg.items()):
            ks = {}
            for a in accs:
                if a not in db:
                    die(1, f"{tsv}: {a} not in {CLEANED}")
                ks[a] = kmers(db[a][1], k)
            for i in range(len(accs)):
                for j in range(i + 1, len(accs)):
                    a, b = accs[i], accs[j]
                    pool[seg].append((containment(ks[a], ks[b]), ds.name, a, b))
    empty = [s for s, v in pool.items() if not v]
    if empty:
        die(1, f"no {STRAIN_SPLIT} pairs for segment(s) {empty}")
    return pool


def real_sources(root: Path, db: dict[str, tuple[str, str]]) -> list[Source]:
    """Every distinct true source over both real splits, with the samples it belongs to."""
    seen: dict[str, tuple[str, list[str]]] = {}
    for split in REAL_SPLITS:
        for d in sorted((root / SPLITS[split].sample_dir).glob("SRR*")):
            fa = d / "true_sources.fasta"
            if not fa.exists():
                die(1, f"{fa} missing — every real sample needs its true_sources.fasta")
            for acc in read_fasta(fa):
                if acc not in db:
                    die(1, f"{fa}: {acc} is not a record of {CLEANED}; a decoy "
                           "must be the twin of a database record")
                seen.setdefault(acc, (split, []))[1].append(d.name)
    out = []
    for acc, (split, samples) in seen.items():
        header, seq = db[acc]
        seg = segment_of(header)
        if seg is None:
            die(1, f"{acc}: header has no 'segment N' field: {header}")
        out.append(Source(acc, header, seq, seg, split, tuple(samples)))
    return sorted(out, key=lambda s: (REAL_SPLITS.index(s.split), s.samples, s.segment, s.acc))


def read_taxmap(path: Path) -> dict[str, str]:
    out = {}
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) >= 2:
                out[p[0]] = p[1]
    return out

def substitute(base: str, rng: random.Random, ti_frac: float) -> str:
    if rng.random() < ti_frac:
        return TRANSITION[base]
    return rng.choice([b for b in BASES if b != base and b != TRANSITION[base]])


def apply(seq: str, subs: dict[int, str]) -> str:
    s = list(seq)
    for pos, b in subs.items():
        s[pos] = b
    return "".join(s)


def draw_subs(seq: str, positions: list[int], rng: random.Random, ti_frac: float) -> dict[int, str]:
    return {p: substitute(seq[p], rng, ti_frac) for p in positions}


def mutation_count(length: int, target_c: float, k: int) -> int:
    return max(1, round(length * (1.0 - target_c ** (1.0 / k))))


def make_family(src: Source, target_c: float, n: int, rng: random.Random,
                p: Params) -> list[Decoy]:
    """`n` decoys of one source; decoy 1 defines M1, each later decoy shares half of it.
    """
    mutable = [i for i, b in enumerate(src.seq) if b in BASES]
    m = mutation_count(len(src.seq), target_c, p.k)
    if m > len(mutable) // 2:
        die(1, f"{src.acc}: {m} substitutions on {len(mutable)} mutable bases")
    src_k = kmers(src.seq, p.k)
    src_k21 = kmers(src.seq, 21)

    best: tuple[float, dict[int, str], str, int] | None = None
    for attempt in range(1, MAX_TRIES + 1):
        pos = sorted(rng.sample(mutable, m))
        subs = draw_subs(src.seq, pos, rng, p.ti_frac)
        seq = apply(src.seq, subs)
        c = containment(src_k, kmers(seq, p.k))
        if best is None or abs(c - target_c) < abs(best[0] - target_c):
            best = (c, subs, seq, attempt)
        if abs(c - target_c) <= TOLERANCE:
            break
    c1, subs1, seq1, tries = best
    fam = [Decoy(f"{PREFIX}{src.acc}_1", src, seq1, target_c, c1,
                 jaccard(src_k21, kmers(seq1, 21)), m, 0, 1.0, tries)]
    k1 = kmers(seq1, p.k)
    m1_positions = sorted(subs1)
    half = m // 2
    fresh_pool = [q for q in mutable if q not in subs1]
    for i in range(2, n + 1):
        best = None
        for attempt in range(1, MAX_TRIES + 1):
            keep = sorted(rng.sample(m1_positions, half))
            fresh = sorted(rng.sample(fresh_pool, m - half))
            subs = {q: subs1[q] for q in keep}
            subs.update(draw_subs(src.seq, fresh, rng, p.ti_frac))
            seq = apply(src.seq, subs)
            ki = kmers(seq, p.k)
            c = containment(src_k, ki)
            if best is None or abs(c - target_c) < abs(best[0] - target_c):
                best = (c, seq, ki, attempt)
            if abs(c - target_c) <= TOLERANCE:
                break
        c, seq, ki, attempt = best
        fam.append(Decoy(f"{PREFIX}{src.acc}_{i}", src, seq, target_c, c,
                         jaccard(src_k21, kmers(seq, 21)), m, half, containment(k1, ki), attempt))
    return fam


def plan_decoys(sources: list[Source], pool: dict[int, list[tuple]], p: Params) -> list[Decoy]:
    out = []
    for src in sources:
        rng = random.Random(f"{p.seed}:{src.acc}")
        target_c = rng.choice(pool[src.segment])[0]
        out.extend(make_family(src, target_c, p.per_source, rng, p))
    return out

def decoy_header(d: Decoy, p: Params) -> str:
    return (f"{d.acc} |decoy of {d.source.acc} target_c={d.target_c:.4f} "
            f"achieved_c={d.achieved_c:.4f} m={d.m} k={p.k} seed={p.seed} "
            f"| {d.source.header.split('|', 1)[-1].strip()}")


def wrap(seq: str, width: int = 70) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def split_records(text: str) -> tuple[str, str]:
    """(non-decoy records, decoy records) of a FASTA text, each verbatim and in file order."""
    base, decoy = [], []
    for rec in re.split(r"(?m)^(?=>)", text):
        if rec:
            (decoy if rec.startswith(">" + PREFIX) else base).append(rec)
    return "".join(base), "".join(decoy)


def split_taxmap(text: str) -> tuple[str, str]:
    base, decoy = [], []
    for line in text.splitlines(keepends=True):
        (decoy if line.startswith(PREFIX) else base).append(line)
    return "".join(base), "".join(decoy)


def _write(path: Path, text: str) -> None:
    """Replace `path` atomically, so an interrupted run never leaves half a reference set."""
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text)
    os.replace(tmp, path)


def write_decoys(idx: Path, base_fa: str, base_tax: str, decoys: list[Decoy],
                 pool: dict[int, list[tuple]], p: Params) -> None:
    taxmap = read_taxmap(idx / TAXMAP)
    recs = "".join(f">{decoy_header(d, p)}\n{wrap(d.seq)}\n" for d in decoys)
    rows = []
    for d in decoys:
        taxid = taxmap.get(d.source.acc)
        if taxid is None:
            die(1, f"{d.source.acc} has no row in {TAXMAP}")
        rows.append(f"{d.acc}\t{taxid}\n")
    _write(idx / CLEANED, base_fa + recs)
    _write(idx / TAXMAP, base_tax + "".join(rows))
    (idx / "decoys.fasta").write_text(recs)

    with open(idx / "decoys.tsv", "w") as f:
        f.write("decoy\tsource\tsplit\tsamples\tsegment\tlength\ttarget_c\tachieved_c\t"
                "achieved_jaccard21\tm\tshared_with_decoy1\tc_to_decoy1\ttries\tseed\n")
        for d in decoys:
            f.write(f"{d.acc}\t{d.source.acc}\t{d.source.split}\t{','.join(d.source.samples)}\t"
                    f"{d.source.segment}\t{len(d.seq)}\t{d.target_c:.4f}\t{d.achieved_c:.4f}\t"
                    f"{d.jaccard21:.4f}\t{d.m}\t{d.shared}\t{d.sibling_c:.4f}\t{d.tries}\t{p.seed}\n")

    with open(idx / "decoy-targets.tsv", "w") as f:
        f.write("segment\tcontainment\tdataset\tacc_a\tacc_b\n")
        for seg in sorted(pool):
            for c, ds, a, b in pool[seg]:
                f.write(f"{seg}\t{c:.4f}\t{ds}\t{a}\t{b}\n")


def remove_decoys(idx: Path, base_fa: str, base_tax: str) -> None:
    _write(idx / CLEANED, base_fa)
    _write(idx / TAXMAP, base_tax)
    for name in ("decoys.fasta", "decoys.tsv", "decoy-targets.tsv"):
        (idx / name).unlink(missing_ok=True)

def report_pool(pool: dict[int, list[tuple]]) -> None:
    print(f"\n{STRAIN_SPLIT} co-occurring-strain containment, per segment (target pool):")
    for seg in sorted(pool):
        v = sorted(c for c, *_ in pool[seg])
        print(f"  segment {seg}: n={len(v)} min={v[0]:.3f} median={statistics.median(v):.3f} "
              f"max={v[-1]:.3f}")
    allv = [c for v in pool.values() for c, *_ in v]
    print(f"  all: n={len(allv)} median={statistics.median(allv):.3f}")


def report_decoys(decoys: list[Decoy]) -> None:
    print(f"\n{len(decoys)} decoys:")
    print("  decoy                         split     seg  len   target  achieved  m   c(D1)  tries")
    for d in decoys:
        print(f"  {d.acc:<29} {d.source.split:<9} {d.source.segment}    {len(d.seq):<5} "
              f"{d.target_c:.3f}   {d.achieved_c:.3f}     {d.m:<3} {d.sibling_c:.3f}  {d.tries}")
    off = [d for d in decoys if abs(d.achieved_c - d.target_c) > TOLERANCE]
    print(f"\n  within +/-{TOLERANCE} of target: {len(decoys) - len(off)}/{len(decoys)}")
    for d in off:
        print(f"    off-target: {d.acc} target={d.target_c:.3f} achieved={d.achieved_c:.3f}")


def verify(db: dict[str, tuple[str, str]], sources: list[Source], decoys: list[Decoy],
           k: int) -> None:
    """Nearest-neighbour containment of every real source against the decoyed database."""
    print("\nverifying: nearest absent neighbours of each real source in the decoyed DB "
          f"(k-merising {len(db)} + {len(decoys)} records)…")
    ks = {acc: kmers(seq, k) for acc, (_, seq) in db.items()}
    ks.update({d.acc: kmers(d.seq, k) for d in decoys})
    print("  source         split     NN1              NN2              NN3")
    nn1 = []
    for s in sources:
        a = ks[s.acc]
        top = sorted(((containment(a, b), acc) for acc, b in ks.items() if acc != s.acc),
                     reverse=True)[:3]
        nn1.append(top[0][0])
        cells = "  ".join(f"{c:.3f} {acc[:10]:<10}" for c, acc in top)
        print(f"  {s.acc:<14} {s.split:<9} {cells}")
    print(f"  NN1 median over real sources: {statistics.median(nn1):.3f}")

def update(dry: bool = False, report: bool = True, check: bool = False) -> bool:
    """Bring the decoys in indexes/ in line with [decoy]. Returns True when the reference
    set changed (so every index built from it is stale)."""
    p = params_from_toml()
    idx = indexes()
    if not (idx / CLEANED).exists():
        die(1, f"{idx / CLEANED} not found; run `premise_bench clean-db` first")
    base_fa, old_fa = split_records((idx / CLEANED).read_text())
    base_tax, old_tax = split_taxmap((idx / TAXMAP).read_text())

    if p.per_source == 0:
        if not old_fa and not old_tax:
            print("decoys: none configured ([decoy] per_source = 0), none present")
            return False
        print(f"decoys: removing {old_fa.count('>')} decoy record(s) ([decoy] per_source = 0)")
        if not dry:
            remove_decoys(idx, base_fa, base_tax)
        return True

    print(f"decoys: {p}")
    db = read_fasta(idx / CLEANED)
    db = {acc: v for acc, v in db.items() if not acc.startswith(PREFIX)}
    pool = target_pool(data_root(), db, p.k)
    sources = real_sources(data_root(), db)
    decoys = plan_decoys(sources, pool, p)
    if report:
        report_pool(pool)
        print(f"\nreal true sources: {len(sources)} "
              f"({', '.join(f'{s}={sum(x.split == s for x in sources)}' for s in REAL_SPLITS)})")
        report_decoys(decoys)
    if check:
        verify(db, sources, decoys, p.k)

    new_fa = "".join(f">{decoy_header(d, p)}\n{wrap(d.seq)}\n" for d in decoys)
    if new_fa == old_fa and old_tax.count("\n") == len(decoys):
        print(f"decoys: {len(decoys)} decoys in {CLEANED} are up to date")
        return False
    print(f"decoys: {'would write' if dry else 'writing'} {len(decoys)} decoys to "
          f"indexes/{CLEANED} ({len(db)} reference records) and indexes/{TAXMAP}")
    if not dry:
        write_decoys(idx, base_fa, base_tax, decoys, pool, p)
    return True
