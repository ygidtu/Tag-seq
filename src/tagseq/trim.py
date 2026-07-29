"""ODN removal — pure Python, handles .gz transparently."""

from __future__ import annotations

import logging
import re
from pathlib import Path

from .fastq import open_fastq, open_for_write

logger = logging.getLogger(__name__)


def _find_odn(seq: str, odn: str) -> int:
    """Return position after the *last* ODN occurrence, or 0."""
    pos = 0
    start = 0
    while True:
        m = re.search(re.escape(odn), seq[start:])
        if not m:
            break
        start += m.end()
        pos = start
    return pos


def remove_odn(
    r1_path: Path,
    r2_path: Path,
    odn_seq: str,
    outdir: Path,
    prefix: str,
) -> tuple[Path, Path, Path]:
    """Remove ODN tag, keep only R1/R2 pairs where R2 contained the tag.

    Returns (r1_out, r2_out, stat_file).
    """
    outdir.mkdir(parents=True, exist_ok=True)
    r1_out = outdir / f"{prefix}.rmODN.R1.fq"
    r2_out = outdir / f"{prefix}.rmODN.R2.fq"
    stat_file = outdir / f"{prefix}.rmODN.stat"

    # ── pass 1: scan R2 ──
    kept: set[str] = set()
    with open_fastq(r2_path) as fin, open(r2_out, "w") as fout:
        while True:
            header = fin.readline()
            if not header:
                break
            seq = fin.readline()
            plus = fin.readline()
            qual = fin.readline()

            # base ID (strip /1 /2 suffix if present)
            base = header.split()[0].lstrip("@")
            base = re.sub(r"/\d+$", "", base)

            pos = _find_odn(seq, odn_seq)
            if pos == 0:
                continue
            new_seq = seq[pos:]
            new_qual = qual[pos:]
            if not new_seq.strip():
                continue

            kept.add(base)
            fout.write(f"{header.split()[0]}\n{new_seq}{plus}{new_qual}")

    # ── pass 2: filter R1 ──
    total = 0
    passed = 0
    with open_fastq(r1_path) as fin, open(r1_out, "w") as fout:
        while True:
            header = fin.readline()
            if not header:
                break
            seq = fin.readline()
            plus = fin.readline()
            qual = fin.readline()
            total += 1
            base = header.split()[0].lstrip("@")
            base = re.sub(r"/\d+$", "", base)
            if base in kept:
                passed += 1
                fout.write(f"{header.split()[0]}\n{seq}{plus}{qual}")

    pct = (passed / max(total, 1)) * 100
    with open(stat_file, "w") as f:
        f.write(f"Raw flagment count: {total}\nRead count with ODN: {passed}\n")

    logger.info("ODN: %d/%d passed (%.2f%%)", passed, total, pct)
    return r1_out, r2_out, stat_file
