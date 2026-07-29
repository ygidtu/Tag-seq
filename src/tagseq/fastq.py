"""FASTQ I/O — handles .fq, .fq.gz, .fastq, .fastq.gz transparently."""

from __future__ import annotations

import gzip
from collections.abc import Iterator
from pathlib import Path
from typing import Optional


def open_fastq(path: str | Path, mode: str = "rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def count_reads(path: str | Path) -> int:
    with open_fastq(path) as f:
        for n, _ in enumerate(f):
            pass
    return (n + 1) // 4


def detect_phred(path: str | Path, n_check: int = 10000) -> int:
    min_q = 200
    with open_fastq(path) as f:
        for i in range(n_check):
            h = f.readline()
            if not h:
                break
            s = f.readline()
            p = f.readline()
            q = f.readline().strip()
            for c in q:
                v = ord(c)
                if v < min_q:
                    min_q = v
    return 64 if min_q >= 64 else 33


def iter_fastq(path: str | Path) -> Iterator[tuple[str, str, str, str]]:
    """Yield (header, sequence, plus, quality) for each record.

    Using raw tuples instead of dataclass avoids allocation overhead.
    """
    with open_fastq(path) as f:
        while True:
            header = f.readline()
            if not header:
                return
            seq = f.readline()
            plus = f.readline()
            qual = f.readline()
            yield header.rstrip("\n"), seq.rstrip("\n"), plus.rstrip("\n"), qual.rstrip("\n")


def write_fastq(records: Iterator[tuple[str, str, str, str]], path: str | Path) -> int:
    count = 0
    with open(path, "w") as f:
        for h, s, p, q in records:
            f.write(f"{h}\n{s}\n{p}\n{q}\n")
            count += 1
    return count


def open_for_write(path: str | Path, mode: str = "w"):
    path = str(path)
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)
