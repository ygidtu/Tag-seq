"""Off-target identification — pure Python (pysam instead of bedtools/samtools)."""

from __future__ import annotations

import logging
from pathlib import Path

logger = logging.getLogger(__name__)


# ── interval helpers ──

def _read_bed(path: Path):
    """Yield (chrom, start, end, *rest) tuples from a BED file."""
    if not path.exists() or path.stat().st_size == 0:
        return
    with open(path) as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 3:
                yield p[0], int(p[1]), int(p[2]), *p[3:]


def _sort_bed(intervals):
    """Sort intervals by (chrom, start)."""
    return sorted(intervals, key=lambda x: (x[0], x[1]))


def _merge_bed(intervals, dist=0):
    """Merge intervals within dist bp."""
    intervals = _sort_bed(intervals)
    if not intervals:
        return
    cur = list(intervals[0])
    for iv in intervals[1:]:
        if iv[0] == cur[0] and iv[1] <= cur[2] + dist:
            cur[2] = max(cur[2], iv[2])
            if len(iv) > 3 and len(cur) > 3:
                cur[3] = f"{cur[3]},{iv[3]}"
        else:
            yield tuple(cur)
            cur = list(iv)
    yield tuple(cur)


def _slop(iv, chromsize_map, bp=40):
    """Extend interval by bp on both sides, clipped to chromosome."""
    chrom, start, end = iv[0], iv[1], iv[2]
    chrom_len = chromsize_map.get(chrom, 10**9)
    return (chrom, max(0, start - bp), min(chrom_len, end + bp), *iv[3:])


def _subtract(intervals, exclude):
    """Remove intervals that overlap with exclude list."""
    for iv in intervals:
        overlap = False
        for ex in exclude:
            if iv[0] == ex[0] and iv[1] < ex[2] and iv[2] > ex[1]:
                overlap = True
                break
        if not overlap:
            yield iv


# ── public API ──

def find_offtargets(
    plus_plus: Path, plus_minus: Path,
    minus_plus: Path, minus_minus: Path,
    control: Path | None,
    blacklist: Path | None,
    grna_seq: str,
    ref_fa: Path,
    chromsize: Path,
    prefix: str,
    outdir: Path,
    min_support: int = 1,
    cut_events: int = 2,
    max_mismatch: int = 6,
) -> Path | None:
    """Full off-target detection — pure Python (pysam + built-ins)."""
    import pysam

    outdir.mkdir(parents=True, exist_ok=True)

    # ── 0. Parse chrom.sizes ──
    chrom_len: dict[str, int] = {}
    with open(chromsize) as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 2:
                chrom_len[p[0]] = int(p[1])

    # ── 1. Read counts from strand-specific BEDs ──
    tcount: dict[str, int] = {}
    for f in [plus_plus, plus_minus, minus_plus, minus_minus]:
        for iv in _read_bed(f):
            if len(iv) >= 5:
                tcount[iv[3]] = int(iv[4])
    if not tcount:
        logger.warning("%s: no sites", prefix)
        return None

    # ── 2. Merge + filter ──
    all_ivs = list(_read_bed(plus_plus)) + list(_read_bed(plus_minus)) + \
              list(_read_bed(minus_plus)) + list(_read_bed(minus_minus))
    all_ivs = _sort_bed(all_ivs)

    # Write all_sites.bed
    (outdir / f"{prefix}.all.sites.bed").write_text(
        "\n".join("\t".join(map(str, iv[:3]) + [iv[3]] if len(iv) > 3 else map(str, iv[:3])) for iv in all_ivs) + "\n"
    )

    merged = list(_merge_bed(all_ivs, dist=10))
    merged_out = outdir / f"{prefix}.all.sites.merged"
    merged_out.write_text(
        "\n".join("\t".join(map(str, iv[:3]) + [iv[3]] if len(iv) > 3 else map(str, iv[:3])) for iv in merged) + "\n"
    )

    # ── 3. Filter by support ──
    confirmed_raw = []
    for iv in merged:
        ids = iv[3].split(",") if len(iv) > 3 else []
        mid = (iv[1] + iv[2]) // 2
        c = {"plus_plus": 0, "plus_minus": 0, "minus_plus": 0, "minus_minus": 0}
        for sid in ids:
            tt = sid.split("_")
            if len(tt) >= 3:
                idx = f"{tt[1]}_{tt[2]}"
                if idx in c:
                    c[idx] += tcount.get(sid, 0)
        if sum(1 for v in c.values() if v >= min_support) >= cut_events:
            confirmed_raw.append(
                (iv[0], mid, mid + 1, iv[3],
                 c["plus_plus"], c["plus_minus"], c["minus_plus"], c["minus_minus"])
            )

    if not confirmed_raw:
        return None

    # ── 4. Subtract blacklist + control ──
    blacklist_ivs = list(_read_bed(blacklist)) if blacklist and blacklist.exists() else []
    control_ivs = list(_read_bed(control)) if control and control.exists() else []

    confirmed = list(_subtract(confirmed_raw, blacklist_ivs))
    if not confirmed:
        return None

    # Extend
    extended = [_slop(iv, chrom_len, 40) for iv in confirmed]

    # Subtract control
    filtered = list(_subtract(extended, control_ivs))
    if not filtered:
        return None

    # ── 5. Extract FASTA ──
    ext_bed = outdir / f"{prefix}.ext.bed"
    with open(ext_bed, "w") as f:
        for iv in filtered:
            f.write(f"{iv[0]}\t{iv[1]}\t{iv[2]}\t{iv[3]}\n")

    ext_fa = outdir / f"{prefix}.ext.fa"
    with pysam.FastaFile(str(ref_fa)) as fa, open(ext_fa, "w") as f:
        for iv in filtered:
            chrom, start, end = iv[0], iv[1], iv[2]
            seq = fa.fetch(chrom, start, end).upper()
            if seq:
                f.write(f">{chrom}:{start}-{end}|{iv[3]}\n{seq}\n")

    if ext_fa.stat().st_size == 0:
        return None

    # ── 6. Smith-Waterman ──
    from Bio import SeqIO
    from Bio.Align import PairwiseAligner

    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 5
    aligner.mismatch_score = -4
    aligner.gap_open_score = -10
    aligner.gap_extend_score = -0.5

    results = []
    for record in SeqIO.parse(ext_fa, "fasta"):
        for strand, tseq in [
            ("+", str(record.seq).upper()),
            ("-", str(record.seq.reverse_complement()).upper()),
        ]:
            alns = aligner.align(grna_seq.upper(), tseq)
            if not alns:
                continue
            best = alns[0]
            mm = sum(1 for a, b in zip(*best.aligned) if a != b)
            if mm <= max_mismatch:
                results.append((record.id, best.score, mm, strand))

    if not results:
        return None

    bed = outdir / f"{prefix}.parsing_water_for_visualization.offtarget.bed"
    with open(bed, "w") as f:
        f.write("#chrom\tstart\tend\tid\tscore\tmm\tstrand\n")
        for rid, sc, mm, strand in sorted(results, key=lambda x: -x[1]):
            chrom = rid.split(":")[0] if ":" in rid else rid
            f.write(f"{chrom}\t0\t0\t{rid}\t{sc}\t{mm}\t{strand}\n")

    logger.info("%s: %d off-targets", prefix, len(results))
    return bed
