from __future__ import annotations
import re
from pathlib import Path
import pysam
from loguru import logger
from tqdm import tqdm
from .fastq import open_fastq_write


def _find_odn(seq: str, odn: str) -> int:
    pos = -1; start = 0
    while True:
        m = re.search(re.escape(odn), seq[start:])
        if not m: break
        start += m.end(); pos = start
    return pos


def _check_tags(seq, tags):
    """Check all tags, return (tag_name, position) for first match or (None, None)."""
    for tag_name, tag_seq in tags.items():
        pos = _find_odn(seq, tag_seq)
        if pos >= 0 and seq[pos + len(tag_seq):].strip():
            return tag_name, pos
    return None, None


def remove_odn(r1_path: Path, r2_path: Path, fwd_tag: str, rev_tag: str, outdir: Path, prefix: str) -> tuple:
    """Remove ODN tag in a single pass through R1/R2. Checks all 4 tags (fwd/rc + rev/rc).

    Reads FASTQ with pysam (supports gzip transparently).
    """
    outdir.mkdir(parents=True, exist_ok=True)
    r1_out = outdir / f"{prefix}.rmODN.R1.fq.gz"
    r2_out = outdir / f"{prefix}.rmODN.R2.fq.gz"
    stat_file = outdir / f"{prefix}.rmODN.stat"

    fwd_rc = fwd_tag.translate(str.maketrans("ATCGatcg", "TAGCtagc"))[::-1]
    rev_rc = rev_tag.translate(str.maketrans("ATCGatcg", "TAGCtagc"))[::-1]

    all_tags = {
        "fwd": fwd_tag,
        "fwd_rc": fwd_rc,
        "rev": rev_tag,
        "rev_rc": rev_rc,
    }

    n_total = 0; n_r2 = 0; n_r1 = 0
    pbar = tqdm(unit="reads", desc=f"  ODN scan ({prefix})", leave=False)
    with pysam.FastxFile(str(r1_path)) as f1, pysam.FastxFile(str(r2_path)) as f2, \
         open_fastq_write(r1_out) as o1, open_fastq_write(r2_out) as o2:
        for e1, e2 in zip(f1, f2):
            n_total += 1
            bid = re.sub(r"/\d+$", "", e1.name)
            h1 = f"@{e1.name}"
            h2 = f"@{e2.name}"
            s1, q1 = e1.sequence, e1.quality
            s2, q2 = e2.sequence, e2.quality

            _, pos = _check_tags(s2, all_tags)
            if pos is not None:
                o1.write(f"{h1}\n{s1}\n+\n{q1}\n")
                o2.write(f"{h2}\n{s2[pos:]}\n+\n{q2[pos:]}\n")
                n_r2 += 1
            else:
                _, pos = _check_tags(s1, all_tags)
                if pos is not None:
                    o1.write(f"{h1}\n{s1[pos:]}\n+\n{q1[pos:]}\n")
                    o2.write(f"{h2}\n{s2}\n+\n{q2}\n")
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
