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


def _process_read(fin):
    h = fin.readline()
    if not h: return None, None, None, None
    return h, fin.readline(), fin.readline(), fin.readline()


def remove_odn(r1_path: Path, r2_path: Path, odn_seq: str, outdir: Path, prefix: str) -> tuple:
    """Remove ODN tag. Checks both R1/R2, both forward and rev-comp."""
    outdir.mkdir(parents=True, exist_ok=True)
    r1_out = outdir / f"{prefix}.rmODN.R1.fq.gz"
    r2_out = outdir / f"{prefix}.rmODN.R2.fq.gz"
    stat_file = outdir / f"{prefix}.rmODN.stat"

    n_total = count_reads(r2_path)
    rev_comp = odn_seq.translate(str.maketrans("ATCGatcg", "TAGCtagc"))[::-1]

    def _scan(path, label):
        kept, seqs, quals = set(), {}, {}
        with open_fastq(path) as fin:
            pbar = tqdm(total=n_total, unit="reads", desc=f"  {label}", leave=False)
            while True:
                h, s, p, q = _process_read(fin)
                if h is None: break
                bid = re.sub(r"/\d+$", "", h.split()[0].lstrip("@"))
                pos = _find_odn(s, odn_seq)
                if pos == 0:
                    pos = _find_odn(s, rev_comp)
                if pos > 0 and s[pos:].strip():
                    kept.add(bid)
                    seqs[bid] = s[pos:]
                    quals[bid] = q[pos:]
                pbar.update(1)
            pbar.close()
        return kept, seqs, quals

    r2_kept, r2_seq, r2_qual = _scan(r2_path, "R2 ODN scan")
    r1_kept, r1_seq, r1_qual = _scan(r1_path, "R1 ODN scan")

    n_r2 = 0; n_r1 = 0
    with open_fastq(r1_path) as f1, open_fastq(r2_path) as f2, \
         open_fastq_write(r1_out) as o1, open_fastq_write(r2_out) as o2:
        pbar = tqdm(total=n_total, unit="reads", desc="  Merge", leave=False)
        while True:
            h1, s1, p1, q1 = _process_read(f1)
            if h1 is None: break
            h2, s2, p2, q2 = _process_read(f2)
            bid = re.sub(r"/\d+$", "", h1.split()[0].lstrip("@"))
            if bid in r2_kept:
                o1.write(f"{h1.split()[0]}\n{s1}{p1}{q1}")
                o2.write(f"{h2.split()[0]}\n{r2_seq[bid]}{p2}{r2_qual[bid]}")
                n_r2 += 1
            elif bid in r1_kept:
                o1.write(f"{h1.split()[0]}\n{r1_seq[bid]}{p1}{r1_qual[bid]}")
                o2.write(f"{h2.split()[0]}\n{s2}{p2}{q2}")
                n_r1 += 1
            pbar.update(1)
        pbar.close()

    total_kept = n_r2 + n_r1
    tag_in_r2 = n_r2 >= n_r1
    pct = (total_kept / max(n_total, 1)) * 100

    with open(stat_file, "w") as f:
        f.write(f"Raw flagment count: {n_total}\nRead count with ODN: {total_kept}\n")
        f.write(f"  Tag in R2: {n_r2}\n  Tag in R1: {n_r1}\n")
    logger.info("ODN: {}/{} passed ({:.2f}%) — R2={}, R1={}",
                total_kept, n_total, pct, n_r2, n_r1)
    return r1_out, r2_out, stat_file, tag_in_r2
