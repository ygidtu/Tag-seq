from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

@dataclass
class PipelineConfig:
    prefix: str
    outdir: Path
    lib_r1: Path
    lib_r2: Path
    forward_tag: str
    reverse_tag: str
    index: Path
    ref: Path
    chromsize: Path
    genome: str = "hg38"
    minlen: int = 50
    readlen: int = 150
    maxins: int = 1000
    threads: int = 16
    umi_len: int = 8
    umi_offset: int = 0
    umi_prefix: str = ""
    min_support_readcount: int = 1
    min_cutting_event_count: int = 2
    max_mismatch: int = 6
    max_gap: int = 2
    max_gap_mismatch: int = 4
    grna_file: Optional[Path] = None
    ctrl: Optional[Path] = None
    adapter: Optional[Path] = None
    blacklist: Optional[Path] = None
    bin_dir: Optional[Path] = None

    def sample_name(self, strand: str) -> str:
        return f"{self.prefix}_{strand}"
    def data_dir(self, sample: str) -> Path:
        return self.outdir / sample / "00datafilter"
    def align_dir(self, sample: str) -> Path:
        return self.outdir / sample / "01alignment"
    def target_dir(self, sample: str) -> Path:
        return self.outdir / sample / "02potentialTargets"
    def find_target_dir(self, grna_id: str) -> Path:
        return self.outdir / f"{grna_id}.find.target"
    @property
    def samples(self) -> list[tuple[str, str]]:
        return [(self.sample_name("plus"), self.forward_tag),
                (self.sample_name("minus"), self.reverse_tag)]

@dataclass
class GrnaEntry:
    id: str; sequence: str; pam: str

@dataclass
class AlignmentStats:
    input_reads: int = 0; unique_mapped: int = 0; multi_mapped: int = 0; too_many_loci: int = 0
    @property
    def total_aligned(self) -> int:
        return self.unique_mapped + self.multi_mapped + self.too_many_loci

@dataclass
class SampleReport:
    name: str = ""; total_reads: int = 0; odn_passed: int = 0
    trim_input: int = 0; trim_passed: int = 0
    alignment: AlignmentStats = field(default_factory=AlignmentStats)
    dedup_unique: int = 0; usable_reads: int = 0
    forward_targets: int = 0; reverse_targets: int = 0
