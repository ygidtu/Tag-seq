from __future__ import annotations
from loguru import logger
import subprocess
from pathlib import Path
from .exceptions import ExternalToolError
from .types import AlignmentStats

def run_star(r1: Path, r2: Path, index: Path, outdir: Path, sample_name: str, threads: int, out_prefix: str = "") -> Path:
    import pysam
    align_dir = outdir / sample_name / "01alignment"
    align_dir.mkdir(parents=True, exist_ok=True)
    prefix = out_prefix or sample_name
    bam = align_dir / f"{prefix}.Aligned.sortedByCoord.out.bam"
    cmd = ["STAR", "--genomeDir", str(index), "--runThreadN", str(threads),
           "--readFilesIn", str(r1), str(r2)]
    if str(r1).endswith(".gz") or str(r2).endswith(".gz"):
        cmd += ["--readFilesCommand", "zcat"]
    cmd += ["--outFileNamePrefix", str(align_dir / f"{prefix}."),
            "--outSAMtype", "BAM", "SortedByCoordinate", "--outReadsUnmapped", "Fastx",
            "--alignIntronMax", "50", "--outFilterScoreMinOverLread", "0.5"]
    log = align_dir / f"{prefix}.star.log"; err = align_dir / f"{prefix}.star.err"
    logger.info("STAR: {}", sample_name)
    with open(log, "w") as lf, open(err, "w") as ef:
        r = subprocess.run(cmd, stdout=lf, stderr=ef, timeout=3600)
    if r.returncode != 0:
        raise ExternalToolError(f"STAR failed (exit={r.returncode}); check {err}")
    pysam.index(str(bam))
    return bam

def parse_star_log(path: Path) -> AlignmentStats:
    s = AlignmentStats()
    if not path.exists(): return s
    with open(path) as f:
        for line in f:
            l = line.strip()
            if "Number of input reads" in l: s.input_reads = _int_after_pipe(l)
            elif "Uniquely mapped reads number" in l: s.unique_mapped = _int_after_pipe(l)
            elif "Number of reads mapped to multiple loci" in l: s.multi_mapped = _int_after_pipe(l)
            elif "Number of reads mapped to too many loci" in l: s.too_many_loci = _int_after_pipe(l)
    return s

def _int_after_pipe(line: str) -> int:
    try: return int(line.split("|")[-1].strip())
    except: return 0
