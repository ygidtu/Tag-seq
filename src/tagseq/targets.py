from __future__ import annotations
from loguru import logger
from pathlib import Path

def _bed_merge(intervals, dist=0):
    if not intervals: return
    cur_chrom, cur_start, cur_end = intervals[0][:3]
    cur_rest = list(intervals[0][3:]) if len(intervals[0]) > 3 else []
    for iv in intervals[1:]:
        c, s, e = iv[0], int(iv[1]), int(iv[2])
        if c == cur_chrom and s <= cur_end + dist:
            cur_end = max(cur_end, e)
        else:
            yield [cur_chrom, cur_start, cur_end] + cur_rest
            cur_chrom, cur_start, cur_end = c, s, e
            cur_rest = list(iv[3:]) if len(iv) > 3 else []
    yield [cur_chrom, cur_start, cur_end] + cur_rest

def _overlap_count(clusters, bed_intervals):
    for cl in clusters:
        cc, cs, ce = cl[:3]; cnt = 0
        for bi in bed_intervals:
            bc, bs, be = bi[:3]
            if bc == cc and bs < ce and be > cs: cnt += 1
        yield cl, cnt

def detect_targets(bam_path: Path, sample_dir: Path, sample_name: str, chromsize: Path, min_qual: int = 30) -> None:
    import pysam
    td = sample_dir / "02potentialTargets"; td.mkdir(parents=True, exist_ok=True)
    p = td / sample_name
    loci_sorted = f"{p}.loci.sorted.bed"; clustered = f"{p}.clustered.region.bed"
    records = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam:
            if read.is_unmapped or read.mate_is_unmapped: continue
            if read.mapping_quality < min_qual: continue
            if not read.is_read2: continue
            strand = "-" if read.is_reverse else "+"
            records.append((read.reference_name, read.reference_start, read.reference_end, read.query_name, read.mapping_quality, strand))
    mids = []
    for chrom, start, end, name, qual, strand in records:
        mid = (start + end) // 2
        mids.append((chrom, mid, mid + 1, name, qual, strand))
    mids.sort(key=lambda x: (x[0], x[1], x[2]))
    with open(loci_sorted, "w") as f:
        for iv in mids: f.write(f"{iv[0]}\t{iv[1]}\t{iv[2]}\t{iv[3]}\t{iv[4]}\t{iv[5]}\n")
    clusters = list(_bed_merge(mids, dist=10))
    with open(clustered, "w") as f:
        for cl in clusters: f.write(f"{cl[0]}\t{cl[1]}\t{cl[2]}\t{len(clusters)}\n")
    plus = [r for r in mids if r[5] == "+"]; minus = [r for r in mids if r[5] == "-"]
    with open(f"{p}.plus.bed", "w") as f:
        for iv in plus: f.write("\t".join(map(str, iv)) + "\n")
    with open(f"{p}.minus.bed", "w") as f:
        for iv in minus: f.write("\t".join(map(str, iv)) + "\n")
    for strand_label, ivs in [("plus", plus), ("minus", minus)]:
        prox_raw = []
        for cl, cnt in _overlap_count(clusters, ivs):
            prox_raw.append((cl[0], cl[1], cl[2], cl[3], cnt, strand_label))
        prox_raw.sort(key=lambda x: (x[0], x[1]))
        with open(f"{p}.{strand_label}.proximal", "w") as f:
            for iv in prox_raw: f.write("\t".join(map(str, iv)) + "\n")
    logger.info("targets {}: {} clusters, {} plus, {} minus", sample_name, len(clusters), len(plus), len(minus))
