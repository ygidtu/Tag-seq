from __future__ import annotations
from loguru import logger
from pathlib import Path
from .fastq import open_fastq, open_fastq_write


def extract_umi(r1_path: Path, r2_path: Path, outdir: Path, prefix: str,
                umi_len: int = 8, umi_offset: int = 0, umi_prefix: str = "") -> tuple[Path, Path]:
    """Extract UMI from R1.

    When umi_prefix is set (e.g. the Universal primer of a modified library),
    only reads whose R1 starts with that prefix get their UMI removed at
    [umi_offset, umi_offset+umi_len); other reads are kept intact (no UMI).
    """
    r1_out = outdir / f"{prefix}.umis.R1.fq.gz"
    r2_out = outdir / f"{prefix}.umis.R2.fq.gz"
    total = 0
    n_umi = 0
    with open_fastq(r1_path) as f1, open_fastq(r2_path) as f2, \
         open_fastq_write(r1_out) as o1, open_fastq_write(r2_out) as o2:
        while True:
            h1 = f1.readline()
            if not h1: break
            s1 = f1.readline(); p1 = f1.readline(); q1 = f1.readline()
            h2 = f2.readline(); s2 = f2.readline(); p2 = f2.readline(); q2 = f2.readline()
            total += 1
            if umi_prefix:
                if s1.startswith(umi_prefix) and len(s1) >= umi_offset + umi_len:
                    umi = s1[umi_offset:umi_offset + umi_len]
                    seq = s1[:umi_offset] + s1[umi_offset + umi_len:]
                    qual = q1[:umi_offset] + q1[umi_offset + umi_len:]
                    n_umi += 1
                else:
                    umi = ""
                    seq, qual = s1, q1
            else:
                if len(s1) >= umi_len:
                    umi = s1[:umi_len]
                    seq, qual = s1[umi_len:], q1[umi_len:]
                    n_umi += 1
                else:
                    umi = ""
                    seq, qual = s1, q1
            o1.write(f"{h1.split()[0]}:UMI_{umi}\n{seq}{p1}{qual}")
            o2.write(f"{h2}{s2}{p2}{q2}")
    logger.info("UMI extracted from {}/{} reads", n_umi, total)
    return r1_out, r2_out


def dedup_umi(bam_path: Path) -> Path:
    import pysam
    dedup_bam = bam_path.parent / bam_path.name.replace(".bam", ".dedup.bam")
    flagstat = dedup_bam.with_suffix(".bam.flagstat")
    seen: set[tuple[str, int, str, bool]] = set()
    kept = total = 0
    with pysam.AlignmentFile(str(bam_path), "rb") as ib, pysam.AlignmentFile(str(dedup_bam), "wb", header=ib.header) as ob:
        for read in ib:
            total += 1
            if read.is_unmapped or read.mate_is_unmapped: continue
            name = read.query_name or ""
            umi = name.split(":UMI_")[-1] if ":UMI_" in name else ""
            key = (read.reference_name or "", read.reference_start, umi, read.is_read1)
            if key not in seen:
                seen.add(key); ob.write(read); kept += 1
    fs = {"total": 0, "mapped": 0, "paired": 0}
    with pysam.AlignmentFile(str(dedup_bam), "rb") as bam:
        for read in bam:
            fs["total"] += 1
            if not read.is_unmapped: fs["mapped"] += 1
            if read.is_paired: fs["paired"] += 1
    with open(flagstat, "w") as f:
        f.write(f"{fs['total']} + 0 in total\n{fs['mapped']} + 0 mapped\n{fs['paired']} + 0 paired in sequencing\n")
    logger.info("Dedup: {} / {} kept ({:.1f}%)", kept, total, (kept/max(total, 1))*100)
    return dedup_bam
