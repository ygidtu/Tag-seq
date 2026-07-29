from __future__ import annotations
import re
from pathlib import Path
from loguru import logger
from tqdm import tqdm
from .fastq import open_fastq, open_fastq_write, count_reads


def _find_odn(seq: str, odn: str) -> int:
    pos = 0; start = 0
    while True:
        m = re.search(re.escape(odn), seq[start:])
        if not m: break
        start += m.end(); pos = start
    return pos


def remove_odn(r1_path: Path, r2_path: Path, odn_seq: str, outdir: Path, prefix: str) -> tuple[Path, Path, Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    r1_out = outdir / f"{prefix}.rmODN.R1.fq.gz"
    r2_out = outdir / f"{prefix}.rmODN.R2.fq.gz"
    stat_file = outdir / f"{prefix}.rmODN.stat"

    r2_total = count_reads(r2_path)

    kept: set[str] = set()
    with open_fastq(r2_path) as fin, open_fastq_write(r2_out) as fout:
        pbar = tqdm(total=r2_total, unit="reads", desc="  R2 ODN scan", leave=False)
        while True:
            header = fin.readline()
            if not header: break
            seq = fin.readline(); plus = fin.readline(); qual = fin.readline()
            base = header.split()[0].lstrip("@")
            base = re.sub(r"/\d+$", "", base)
            pos = _find_odn(seq, odn_seq)
            if pos > 0:
                new_seq = seq[pos:]; new_qual = qual[pos:]
                if new_seq.strip():
                    kept.add(base)
                    fout.write(f"{header.split()[0]}\n{new_seq}{plus}{new_qual}")
            pbar.update(1)
        pbar.close()

    total = 0; passed = 0
    with open_fastq(r1_path) as fin, open_fastq_write(r1_out) as fout:
        pbar = tqdm(total=r2_total, unit="reads", desc="  R1 filter", leave=False)
        while True:
            header = fin.readline()
            if not header: break
            seq = fin.readline(); plus = fin.readline(); qual = fin.readline()
            total += 1
            base = header.split()[0].lstrip("@")
            base = re.sub(r"/\d+$", "", base)
            if base in kept:
                passed += 1
                fout.write(f"{header.split()[0]}\n{seq}{plus}{qual}")
            pbar.update(1)
        pbar.close()

    pct = (passed / max(total, 1)) * 100
    with open(stat_file, "w") as f:
        f.write(f"Raw flagment count: {total}\nRead count with ODN: {passed}\n")
    logger.info("ODN: {}/{} passed ({:.2f}%)", passed, total, pct)
    return r1_out, r2_out, stat_file
