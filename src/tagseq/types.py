from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

@dataclass
class LibrarySpec:
    """文库 5' 结构定义 (见 library.py 与 TAG-seq信息.xlsx)。

    改造后文库 R1: primer + UMI + barcode + junction + genome + Tag
    原版文库   R1: P5(primer) + UMI + barcode + junction + genome + Tag + P7
    """
    lib_type: str = "modified"
    primer: str = ""            # 5' 固定前缀 (universal primer / P5)
    umi_len: int = 8
    umi_pattern: str = ""       # IUPAC 简并模式, 如 NNWNNWNN
    barcode: str = ""           # 期望 barcode; 空=不校验
    barcode_len: int = 8
    junction: str = ""          # 5' 侧 adapter junction
    tag_len: int = 25
    strip_5p: bool = True       # 比对前剥离 primer+umi+barcode+junction
    strict: bool = False        # 校验失败: True=报错, False=告警计数


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
    lib_type: str = "legacy"
    odn_mode: str = "rmODN"  # 比对入口过滤口径: rmODN=tag 锚定(既有口径) / mixed=parse_mixed(实验)
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
    lib_spec: Optional["LibrarySpec"] = None

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
