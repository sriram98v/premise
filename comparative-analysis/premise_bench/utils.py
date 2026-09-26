"""Shared helpers"""
from __future__ import annotations

import gzip
import re
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import NoReturn

def die(code: int, *lines: str) -> NoReturn:
    """Print to stderr and exit with `code`.
    """
    for line in lines:
        print(line, file=sys.stderr, flush=True)
    raise SystemExit(code)

def nonempty(p: Path) -> bool:
    """`[ -s "$p" ]` — exists AND non-zero. Path.exists() alone is not the same test."""
    p = Path(p)
    return p.is_file() and p.stat().st_size > 0


def open_maybe_gz(path):
    """Open `path` for text reading, decompressing iff its name ends `.gz`.
    """
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def open_or_gz(path: Path):
    """Open `path`, falling back to `path`.gz — same format, different container.
    Raises FileNotFoundError when neither exists."""
    path = Path(path)
    if path.exists():
        return open(path)
    gz = path.with_suffix(path.suffix + ".gz")
    if gz.exists():
        return gzip.open(gz, "rt")
    raise FileNotFoundError(f"{path}[.gz] not found")


def warn(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def count_pairs(path: Path) -> int:
    """Read pairs in a FASTQ: integer division, so a truncated file floors silently rather
    than reporting a fractional pair."""
    with open_maybe_gz(path) as f:
        return sum(1 for _ in f) // 4

# Labels of unclassified reads.
UNCLASSIFIED_LABELS = frozenset({
    "", "*", "-", "unclassified", "not aligned", "not_aligned", "mapfail",
})

def is_unclassified(ref: str | None) -> bool:
    """True when a method's reference field means "this read got no assignment"."""
    return ref is None or ref.strip().lower() in UNCLASSIFIED_LABELS

# Regex for syn read true source from header
_ISS_READ_SUFFIX = re.compile(r"_\d+(?:_\d+)?$")

def iss_source(read_id: str) -> str | None:
    """The source accession encoded in an InSilicoSeq read name
    """
    return _ISS_READ_SUFFIX.sub("", read_id)


def strain_of(ref: str | None) -> str | None:
    """Collapse a predicted segment/accession to its influenza strain, for real-mix scoring.
    """
    if not ref:
        return None
    if ref.startswith("PR8"):
        return "PR8"
    if ref.startswith("WSN33"):
        return "WSN33"
    return "other"

REAL_MIX_STRAINS = frozenset({"PR8", "WSN33"})

def real_truth(path: Path, split: str) -> tuple[dict, dict, int]:
    """(truth_reads, truth_counts, n_reads) from a real sample's truth_assignments.tsv[.gz].
    """
    truth_seg = {}
    with open_or_gz(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0] and p[1]:
                truth_seg[p[0]] = p[1]
    counts = dict(Counter(truth_seg.values()))
    if split != "real-mix":
        return truth_seg, counts, len(truth_seg)
    truth = {q: strain_of(r) for q, r in truth_seg.items()}
    bad = set(truth.values()) - REAL_MIX_STRAINS
    if bad:
        raise ValueError(f"{path}: {sorted(bad)} — every reference must be a PR8_/WSN33_ "
                         "segment on real-mix, or the strain truth collapses to 'other'")
    return truth, counts, len(truth)


def synthetic_truth(r1: Path) -> dict[str, str | None]:
    """{read_id: source accession} for an InSilicoSeq sample, from its R1 FASTQ headers."""
    truth = {}
    with open_maybe_gz(r1) as f:
        for i, line in enumerate(f):
            if i % 4:
                continue
            rid = line[1:].rstrip("\n").rsplit("/", 1)[0]
            truth[rid] = iss_source(rid)
    return truth


def parse_timemem(path: Path) -> tuple[float | None, float | None]:
    """(wall_seconds, rss_gb) from a `/usr/bin/time -v` file.
    """
    if not Path(path).exists():
        return None, None
    wall = rss = None
    for line in Path(path).read_text(errors="ignore").splitlines():
        if "Elapsed (wall clock)" in line:
            parts = line.split(": ")[-1].strip().replace("s", "").split(":")
            try:
                wall = float(parts[-1])
                if len(parts) >= 2:
                    wall += float(parts[-2]) * 60
                if len(parts) >= 3:
                    wall += float(parts[-3]) * 3600
            except ValueError:
                pass
        elif "Maximum resident set size" in line:
            try:
                # kB to GB
                rss = int(line.split(": ")[-1].strip()) / 1048576.0
            except ValueError:
                pass
    return wall, rss

@dataclass(frozen=True)
class FastaRecord:
    """One fasta record: `header` without the leading '>', `seq` as given."""

    header: str
    seq: str

    @property
    def id(self) -> str:
        return self.header.split(maxsplit=1)[0] if self.header else ""


SEG_RE = re.compile(r"segment (\d)")


def iter_fasta(path: Path):
    """Yield (accession, header without '>', uppercased sequence) per record, in file order."""
    import fastapy  # only the reference-building commands need it

    for rec in fastapy.parse(path):
        yield rec.id, rec.description[1:], rec.seq.upper()


def iter_fasta_headers(path: Path):
    """Yield each header line without the leading '>' or the newline, in file order."""
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                yield line[1:].rstrip("\n")


_GAP_CHARS = "-._"


def kmer_set(seq: str, k: int) -> set[str]:
    """Forward-strand unique k-mers of `seq` after uppercasing and stripping gap chars (-._).
    """
    if k < 1:
        raise ValueError(f"k must be >= 1, got {k}")
    s = seq.upper().translate(str.maketrans("", "", _GAP_CHARS))
    return {s[i : i + k] for i in range(len(s) - k + 1)}


def kmer_containment(rec_a: FastaRecord, rec_b: FastaRecord, k: int) -> float:
    r"""Compute k-mer containment as $\gamma_k(P,Q):=\frac{|S_k(P)\cap S_k(Q)|}{\min(|S_k(P)|,|S_k(Q)|)}$, where $S_k$ is the set of unique k-mers in a string
    """
    a = kmer_set(rec_a.seq, k)
    b = kmer_set(rec_b.seq, k)
    denom = min(len(a), len(b))
    if denom == 0:
        return 0.0
    return len(a & b) / denom
