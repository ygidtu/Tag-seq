from __future__ import annotations
from loguru import logger
from pathlib import Path
from .fastq import open_fastq, open_fastq_write


def extract_umi(r1_path: Path, r2_path: Path, outdir: Path, prefix: str,
                umi_len: int = 8, umi_offset: int = 0, umi_prefix: str = "") -> tuple[Path, Path]:
    """Extract UMI from R1.  (生产口径, 与 TAG-1/3/5 既有结果逐位一致)

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


def extract_umi_spec(r1_path: Path, r2_path: Path, outdir: Path, prefix: str, spec) -> tuple[Path, Path]:
    """按 LibrarySpec 从 R1 提取 UMI/barcode 并剥离 5' 人工序列 (实验路径, 非既有口径)。

    - R1 以 primer 开头: 提取 UMI(+barcode), 剥离 primer+UMI+barcode+junction,
      保留 insert(genome+tag)用于比对; UMI 写入 read 名供 dedup 使用。
    - R1 不以 primer 开头: 视为反向文库/背景 read, 不剥离, 计入 primer_mismatch。
    对 UMI 简并模式 / barcode / junction 做校验并计数; spec.strict 时对异常报错。
    """
    from .library import parse_5p
    r1_out = outdir / f"{prefix}.umis.R1.fq.gz"
    r2_out = outdir / f"{prefix}.umis.R2.fq.gz"
    total = parsed = primer_mm = umi_bad = bc_bad = junc_bad = 0
    with open_fastq(r1_path) as f1, open_fastq(r2_path) as f2, \
         open_fastq_write(r1_out) as o1, open_fastq_write(r2_out) as o2:
        while True:
            h1 = f1.readline()
            if not h1: break
            s1 = f1.readline().rstrip("\n"); p1 = f1.readline(); q1 = f1.readline().rstrip("\n")
            h2 = f2.readline(); s2 = f2.readline(); p2 = f2.readline(); q2 = f2.readline()
            total += 1
            rp = parse_5p(s1, spec)
            if rp.primer_ok:
                parsed += 1
                if rp.umi_ok is False: umi_bad += 1
                if rp.barcode_ok is False: bc_bad += 1
                if rp.junction_ok is False: junc_bad += 1
                if spec.strip_5p:
                    seq, qual = s1[rp.fp_off:], q1[rp.fp_off:]
                else:
                    seq, qual = s1[:rp.umi_start] + s1[rp.umi_end:], q1[:rp.umi_start] + q1[rp.umi_end:]
                name = h1.rstrip("\n").split()[0]
                o1.write(f"{name}:UMI_{rp.umi}\n{seq}\n{p1}{qual}\n")
            else:
                primer_mm += 1
                o1.write(f"{h1}{s1}\n{p1}{q1}\n")
            o2.write(f"{h2}{s2}{p2}{q2}")
    logger.info("UMI: primer_ok {}/{}; primer_mismatch={}; umi_bad={}; barcode_bad={}; junction_bad={}",
                parsed, total, primer_mm, umi_bad, bc_bad, junc_bad)
    if spec.strict and (umi_bad or bc_bad or junc_bad):
        raise ValueError(f"{prefix}: 严格模式下文库结构校验失败 "
                         f"(umi_bad={umi_bad}, barcode_bad={bc_bad}, junction_bad={junc_bad})")
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
