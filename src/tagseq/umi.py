"""UMI handling — extract & dedup, pure Python (pysam instead of samtools subprocess)."""

from __future__ import annotations

import logging
from pathlib import Path

from .fastq import open_fastq

logger = logging.getLogger(__name__)


def extract_umi(
    r1_path: Path,
    r2_path: Path,
    outdir: Path,
    prefix: str,
    umi_len: int = 8,
) -> tuple[Path, Path]:
    """Extract first N bases of R1 as UMI, append to read header."""
    r1_out = outdir / f"{prefix}.Trim.R1.umis.fq"
    r2_out = outdir / f"{prefix}.Trim.R2.umis.fq"

    total = 0
    with open_fastq(r1_path) as f1, open_fastq(r2_path) as f2, \
         open(r1_out, "w") as o1, open(r2_out, "w") as o2:
        while True:
            h1 = f1.readline()
            if not h1:
                break
            s1 = f1.readline()
            p1 = f1.readline()
            q1 = f1.readline()
            h2 = f2.readline()
            s2 = f2.readline()
            p2 = f2.readline()
            q2 = f2.readline()
            total += 1
            umi = s1[:umi_len]
            o1.write(f"{h1.split()[0]}:UMI_{umi}\n{s1[umi_len:]}{p1}{q1[umi_len:]}")
            o2.write(f"{h2}{s2}{p2}{q2}")

    logger.info("UMI extracted from %d reads", total)
    return r1_out, r2_out


def dedup_umi(bam_path: Path) -> Path:
    """Remove PCR duplicates by UMI tag using pysam."""
    import pysam

    dedup_bam = bam_path.parent / bam_path.name.replace(".bam", ".dedup.bam")
    flagstat = dedup_bam.with_suffix(".bam.flagstat")

    seen: set[tuple[str, int, str, bool]] = set()
    kept = total = 0

    with pysam.AlignmentFile(str(bam_path), "rb") as ib, \
         pysam.AlignmentFile(str(dedup_bam), "wb", header=ib.header) as ob:
        for read in ib:
            total += 1
            if read.is_unmapped or read.mate_is_unmapped:
                continue
            name = read.query_name or ""
            umi = name.split(":UMI_")[-1] if ":UMI_" in name else ""
            key = (read.reference_name or "", read.reference_start, umi, read.is_read1)
            if key not in seen:
                seen.add(key)
                ob.write(read)
                kept += 1

    # Compute flagstat with pysam
    fs = {"total": 0, "mapped": 0, "paired": 0}
    with pysam.AlignmentFile(str(dedup_bam), "rb") as bam:
        for read in bam:
            fs["total"] += 1
            if not read.is_unmapped:
                fs["mapped"] += 1
            if read.is_paired:
                fs["paired"] += 1
    with open(flagstat, "w") as f:
        f.write(f"{fs['total']} + 0 in total\n")
        f.write(f"{fs['mapped']} + 0 mapped\n")
        f.write(f"{fs['paired']} + 0 paired in sequencing\n")

    logger.info("Dedup: %d / %d kept (%.1f%%)", kept, total, (kept / max(total, 1)) * 100)
    return dedup_bam
