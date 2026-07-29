"""STAR alignment wrapper + log parser — pure Python (no samtools subprocess)."""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path

from .exceptions import ExternalToolError
from .types import AlignmentStats

logger = logging.getLogger(__name__)


def run_star(
    r1: Path,
    r2: Path,
    index: Path,
    outdir: Path,
    sample_name: str,
    threads: int,
    out_prefix: str = "",
) -> Path:
    """Run STAR alignment; returns path to sorted BAM.

    STAR is the only external tool still required (no Python equivalent).
    Indexing uses pysam.
    """
    import pysam

    align_dir = outdir / sample_name / "01alignment"
    align_dir.mkdir(parents=True, exist_ok=True)
    prefix = out_prefix or sample_name
    bam = align_dir / f"{prefix}.Aligned.sortedByCoord.out.bam"

    cmd = [
        "STAR",
        "--genomeDir", str(index),
        "--runThreadN", str(threads),
        "--readFilesIn", str(r1), str(r2),
        "--outFileNamePrefix", str(align_dir / f"{prefix}."),
        "--outSAMtype", "BAM", "SortedByCoordinate",
        "--outReadsUnmapped", "Fastx",
        "--alignIntronMax", "50",
        "--outFilterScoreMinOverLread", "0.5",
    ]
    log = align_dir / f"{prefix}.star.log"
    err = align_dir / f"{prefix}.star.err"

    logger.info("STAR: %s", sample_name)
    with open(log, "w") as lf, open(err, "w") as ef:
        r = subprocess.run(cmd, stdout=lf, stderr=ef)
    if r.returncode != 0:
        raise ExternalToolError(f"STAR failed (exit={r.returncode}); check {err}")

    # Index with pysam instead of samtools
    pysam.index(str(bam))
    return bam


def parse_star_log(path: Path) -> AlignmentStats:
    s = AlignmentStats()
    if not path.exists():
        return s
    with open(path) as f:
        for line in f:
            l = line.strip()
            if "Number of input reads" in l:
                s.input_reads = _int_after_pipe(l)
            elif "Uniquely mapped reads number" in l:
                s.unique_mapped = _int_after_pipe(l)
            elif "Number of reads mapped to multiple loci" in l:
                s.multi_mapped = _int_after_pipe(l)
            elif "Number of reads mapped to too many loci" in l:
                s.too_many_loci = _int_after_pipe(l)
    return s


def compute_flagstat(bam_path: Path) -> dict[str, int]:
    """Compute BAM statistics with pysam instead of samtools flagstat."""
    import pysam
    stats = {"total": 0, "mapped": 0, "paired": 0, "proper_pair": 0}
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam:
            stats["total"] += 1
            if not read.is_unmapped:
                stats["mapped"] += 1
            if read.is_paired:
                stats["paired"] += 1
                if read.is_proper_pair:
                    stats["proper_pair"] += 1
    return stats


def _int_after_pipe(line: str) -> int:
    try:
        return int(line.split("|")[-1].strip())
    except (ValueError, IndexError):
        return 0
