from __future__ import annotations
import gzip
from pathlib import Path

def open_fastq(path: str | Path, mode: str = "rt"):
    p = str(path)
    return gzip.open(p, mode) if p.endswith(".gz") else open(p, mode)

def open_fastq_write(path: str | Path, mode: str = "wt"):
    p = str(path)
    return gzip.open(p, mode) if p.endswith(".gz") else open(p, mode)

def count_reads(path: str | Path) -> int:
    with open_fastq(path) as f:
        for n, _ in enumerate(f): pass
    return (n + 1) // 4

def detect_phred(path: str | Path, n_check: int = 10000) -> int:
    min_q = 200
    with open_fastq(path) as f:
        for i in range(n_check):
            h = f.readline()
            if not h: break
            s = f.readline(); p = f.readline(); q = f.readline().strip()
            for c in q:
                v = ord(c)
                if v < min_q: min_q = v
    return 64 if min_q >= 64 else 33
