"""Pipeline statistics — pure Python (pysam instead of samtools flagstat)."""

from __future__ import annotations

import re
from pathlib import Path

from .alignment import parse_star_log
from .types import SampleReport  # noqa: F401


def _read_odn(path: Path) -> tuple[int, int, int, int]:
    """Read rmODN stat, return (total, passed, r2_count, r1_count)."""
    t = p = r2 = r1 = 0
    if not path.exists():
        return t, p, r2, r1
    with open(path) as f:
        for line in f:
            m = re.search(r"Raw flagment count:\s+(\d+)", line)
            if m: t = int(m.group(1))
            m = re.search(r"Read count with ODN:\s+(\d+)", line)
            if m: p = int(m.group(1))
            m = re.search(r"Tag in R2:\s+(\d+)", line)
            if m: r2 = int(m.group(1))
            m = re.search(r"Tag in R1:\s+(\d+)", line)
            if m: r1 = int(m.group(1))
    return t, p, r2, r1


def _read_adapter_setting(path: Path) -> int:
    """Read 'Total read pairs processed' from a cutadapt report file."""
    if not path.exists():
        return 0
    with open(path) as f:
        for line in f:
            m = re.search(r"Total read pairs processed:\s+([\d,]+)", line)
            if m: return int(m.group(1).replace(",", ""))
    return 0


def _bam_stat(bam_path: Path) -> dict:
    """Get BAM statistics with pysam."""
    import pysam
    s = {"total": 0, "mapped": 0, "paired": 0}
    if not bam_path.exists():
        return s
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam:
            s["total"] += 1
            if not read.is_unmapped:
                s["mapped"] += 1
            if read.is_paired:
                s["paired"] += 1
    return s


def _count_min(path: Path, col: int = 5, minv: int = 5) -> int:
    if not path.exists():
        return 0
    c = 0
    with open(path) as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= col and p[col - 1].isdigit() and int(p[col - 1]) >= minv:
                c += 1
    return c


def generate_report(prefix: str, outdir: Path) -> list[SampleReport]:
    reports = []
    for sample in [f"{prefix}_plus", f"{prefix}_minus"]:
        sr = SampleReport(name=sample)
        dd = outdir / sample / "00datafilter"
        ad = outdir / sample / "01alignment"
        td = outdir / sample / "02potentialTargets"

        sr.total_reads, sr.odn_passed, r2c, r1c = _read_odn(dd / f"{sample}.rmODN.stat")
        sr.trim_input = _read_adapter_setting(dd / f"{sample}.trim_report.txt")
        sr.alignment = parse_star_log(ad / f"{sample}.Log.final.out")
        sr.trim_passed = sr.alignment.input_reads

        dedup_bam = ad / f"{sample}.Aligned.sortedByCoord.out.dedup.bam"
        dedup_stat = _bam_stat(dedup_bam)
        sr.dedup_unique = dedup_stat["paired"]

        filter_bam = td / f"{sample}.Aligned.sortedByCoord.out.dedup.single.filtered.bam"
        filter_stat = _bam_stat(filter_bam)
        sr.usable_reads = filter_stat["mapped"]

        sr.forward_targets = _count_min(td / f"{sample}.plus.proximal")
        sr.reverse_targets = _count_min(td / f"{sample}.minus.proximal")
        reports.append(sr)
    return reports


def print_report(reports: list[SampleReport]) -> None:
    import sys
    for sr in reports:
        ta = sr.alignment.total_aligned
        print(f"\n{'='*50}\n  {sr.name}\n{'='*50}", file=sys.stderr)
        print(f"  ODN pass:    {sr.odn_passed}/{sr.total_reads} ({sr.odn_passed/max(sr.total_reads,1)*100:.2f}%)", file=sys.stderr)
        print(f"  Trim:        {sr.trim_passed}/{sr.trim_input} ({sr.trim_passed/max(sr.trim_input,1)*100:.2f}%)", file=sys.stderr)
        print(f"  Aligned:     {ta}/{sr.trim_passed} ({ta/max(sr.trim_passed,1)*100:.2f}%)", file=sys.stderr)
        print(f"  Unique:      {sr.dedup_unique}", file=sys.stderr)
        print(f"  Dup rate:    {(ta-sr.dedup_unique)/max(ta,1)*100:.1f}%", file=sys.stderr)
        print(f"  Usable:      {sr.usable_reads}", file=sys.stderr)
        print(f"  Targets+:    {sr.forward_targets}", file=sys.stderr)
        print(f"  Targets-:    {sr.reverse_targets}", file=sys.stderr)
