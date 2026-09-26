"""The five benchmark splits"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .config import data_root

SYNTH_DS = ("Dataset-1", "Dataset-2", "Dataset-3", "Dataset-4")
REAL_DS = ("SRR31013463", "SRR31013465", "SRR31013467", "SRR31013473")
REALMIXED_DS = ("SRR3360139", "SRR3360140", "SRR3360145", "SRR3360146")


@dataclass(frozen=True)
class Split:
    code: str
    sub: str
    """Sub-path shared by samples/<sub>/ and results/<method>/<sub>/, e.g. "real/isolate"."""
    reads: str
    """Read-file name inside a dataset directory; `{base}` is the dataset, `{n}` the mate."""
    datasets: tuple[str, ...]

    @property
    def real(self) -> bool:
        return self.sub.startswith("real/")

    @property
    def synthetic(self) -> bool:
        return not self.real

    @property
    def sample_dir(self) -> str:
        """samples/<sub>, relative to the data root."""
        return f"samples/{self.sub}"

    def read_name(self, ds: str, mate: int) -> str:
        return self.reads.format(base=ds, n=mate)

    def read_path(self, ds: str, mate: int = 1) -> Path:
        """Absolute path to one mate of a dataset's input FASTQ."""
        return data_root() / self.sample_dir / ds / self.read_name(ds, mate)


_SYN_READS = "reads_R{n}.fastq"
_REAL_READS = "{base}_{n}-filtered.ca.fastq"

SPLITS: dict[str, Split] = {s.code: s for s in (
    Split("syn-iso", "synthetic/isolate", _SYN_READS, SYNTH_DS),
    Split("syn-mix-sub", "synthetic/mixed", _SYN_READS, SYNTH_DS),
    Split("syn-mix-strain", "synthetic/mixed-subtype", _SYN_READS, SYNTH_DS),
    Split("real-iso", "real/isolate", _REAL_READS, REAL_DS),
    Split("real-mix", "real/mixed", _REAL_READS, REALMIXED_DS),
)}
SPLIT_ORDER = tuple(SPLITS)


def split_read_path(split: str, ds: str, mate: int = 1) -> Path:
    """Absolute path to one mate of a dataset's input FASTQ."""
    return SPLITS[split].read_path(ds, mate)
