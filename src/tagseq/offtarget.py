"""Off-target identification."""

from __future__ import annotations
from pathlib import Path
from loguru import logger
from tqdm import tqdm


def _read_bed(path: Path, tag: str = ""):
    """Yield (chrom, start, end, id, count, tag) from proximal BED."""
    if not path.exists() or path.stat().st_size == 0:
        return
    with open(path) as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 5:
                yield p[0], int(p[1]), int(p[2]), p[3], int(p[4]), tag


def _merge(intervals, max_dist=10):
    """Merge overlapping/nearby intervals, sum counts."""
    if not intervals:
        return []
    s = sorted(intervals, key=lambda x: (x[0], x[1]))
    m = [list(s[0])]
    for iv in s[1:]:
        if iv[0] == m[-1][0] and iv[1] <= m[-1][2] + max_dist:
            m[-1][2] = max(m[-1][2], iv[2])
            m[-1][3] += f",{iv[3]}"
            m[-1][4] += iv[4]
        else:
            m.append(list(iv))
    return [tuple(x) for x in m]


def _slop(iv, chrom_len, bp=40):
    cl = chrom_len.get(iv[0], 10**9)
    return (iv[0], max(0, iv[1] - bp), min(cl, iv[2] + bp)) + iv[3:]


def _confirm_sites(merged, combined, min_support, cut_events):
    """扫描线累加每个 merged 区间内 combined 的 tag 计数, 仅保留满足
    cut_events 个方向 (count>=min_support) 的位点。

    与原 O(merged*combined) 双重循环严格等价: 重叠条件为
    civ[1] < end and civ[2] > start。merged 已按 (chrom,start) 排序且互不
    重叠 (start/end 递增), combined 按 (chrom,start) 排序; 用按 civ.end 的
    最小堆维护当前窗口 (已入堆 => start<end; 未弹出 => end>start)。
    """
    import heapq
    comb_by_chrom = {}
    for civ in combined:
        comb_by_chrom.setdefault(civ[0], []).append(civ)
    for lst in comb_by_chrom.values():
        lst.sort(key=lambda x: (x[1], x[2]))
    merged_by_chrom = {}
    for idx, iv in enumerate(merged):
        merged_by_chrom.setdefault(iv[0], []).append(idx)
    confirmed = [None] * len(merged)
    for chrom, m_idx in merged_by_chrom.items():
        comb = comb_by_chrom.get(chrom, [])
        ci = 0
        heap = []  # (civ_end, civ_count, civ_tag)
        for idx in m_idx:
            iv = merged[idx]
            start, end = iv[1], iv[2]
            mid = (start + end) // 2
            while ci < len(comb) and comb[ci][1] < end:
                heapq.heappush(heap, (comb[ci][2], comb[ci][4], comb[ci][5]))
                ci += 1
            while heap and heap[0][0] <= start:
                heapq.heappop(heap)
            tag_counts = {"plus_plus": 0, "plus_minus": 0,
                          "minus_plus": 0, "minus_minus": 0}
            for _, cnt, tag in heap:
                if tag in tag_counts:
                    tag_counts[tag] += cnt
            if sum(1 for v in tag_counts.values() if v >= min_support) >= cut_events:
                # 保留 cluster 的 [start,end] (而非中点), 使 ext 窗口覆盖整个 cluster,
                # 避免宽 cluster 的 protospacer 落在 mid±150 之外导致 on-target 漏检。
                confirmed[idx] = (chrom, start, end, iv[3],
                                  tag_counts["plus_plus"], tag_counts["plus_minus"],
                                  tag_counts["minus_plus"], tag_counts["minus_minus"])
    return [iv for iv in confirmed if iv is not None]


def find_offtargets(
    plus_plus, plus_minus, minus_plus, minus_minus,
    control, blacklist, grna_seq, ref_fa, chromsize,
    prefix, outdir,
    min_support=1, cut_events=2, max_mismatch=6, pam="",
):
    import pysam
    from Bio import SeqIO
    from Bio.Align import PairwiseAligner

    outdir.mkdir(parents=True, exist_ok=True)

    chrom_len = {}
    with open(chromsize) as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 2:
                chrom_len[p[0]] = int(p[1])

    labels = {"plus_plus": plus_plus, "plus_minus": plus_minus,
              "minus_plus": minus_plus, "minus_minus": minus_minus}
    all_by_tag = {k: list(_read_bed(p, k)) for k, p in labels.items()}
    total = sum(len(v) for v in all_by_tag.values())
    if total == 0:
        logger.warning("{}: no sites", prefix)
        return None
    logger.info("{}: {} total sites", prefix, total)

    combined = []
    for tag, ivs in all_by_tag.items():
        for iv in ivs:
            combined.append(iv)
    merged = _merge(combined, 10)

    confirmed = _confirm_sites(merged, combined, min_support, cut_events)

    if not confirmed:
        logger.warning("{}: no confirmed sites (need {} orients >= {} reads)", prefix, cut_events, min_support)
        return None
    logger.info("{}: {} confirmed", prefix, len(confirmed))

    confirmed_raw = outdir / f"{prefix}.all.sites.merged.confirmed.raw"
    with open(confirmed_raw, "w") as f:
        for iv in confirmed:
            f.write("\t".join(map(str, iv)) + "\n")

    def _nchrom(c: str) -> str:
        return c[3:] if c.startswith("chr") else c

    # 预读 blacklist/control 区间(归一化 chrom)。稍后仅用于过滤 off-target(others);
    # on-target(perfect)豁免过滤, 避免真靶点被黑名单/对照误杀。
    excl = {}
    for src in (blacklist, control):
        for e in (list(_read_bed(src)) if src and Path(src).exists() else []):
            excl.setdefault(_nchrom(e[0]), []).append((int(e[1]), int(e[2])))
    for v in excl.values():
        v.sort()

    def _in_excl(chrom, start, end):
        for ss, ee in excl.get(_nchrom(chrom), ()):
            if start < ee and end > ss:
                return True
        return False

    extended = [_slop(iv, chrom_len, 150) for iv in confirmed]
    ext_bed = outdir / f"{prefix}.ext.bed"
    with open(ext_bed, "w") as f:
        for iv in extended:
            f.write(f"{iv[0]}\t{iv[1]}\t{iv[2]}\t{iv[3]}\n")

    ext_fa = outdir / f"{prefix}.ext.fa"
    with pysam.FastaFile(str(ref_fa)) as fa, open(ext_fa, "w") as f:
        for iv in extended:
            seq = fa.fetch(iv[0], iv[1], iv[2]).upper()
            if seq:
                # iv[4:] are [plus_plus, plus_minus, minus_plus, minus_minus] counts
                total_count = sum(iv[4:8]) if len(iv) >= 8 else 0
                f.write(f">{iv[0]}:{iv[1]}-{iv[2]}|cnt={total_count}|{iv[3]}\n{seq}\n")
    if ext_fa.stat().st_size == 0:
        return None

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 5
    aligner.mismatch_score = -4
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    # Semi-global: gRNA must align end-to-end, target ends are free.
    aligner.open_left_insertion_score = 0
    aligner.extend_left_insertion_score = 0
    aligner.open_right_insertion_score = 0
    aligner.extend_right_insertion_score = 0

    results = []
    records = list(SeqIO.parse(ext_fa, "fasta"))
    pbar = tqdm(total=len(records) * 2, unit="align", desc="  Smith-Waterman", leave=False)
    for record in records:
        for strand_label, tseq in [
            ("+", str(record.seq).upper()),
            ("-", str(record.seq.reverse_complement()).upper()),
        ]:
            alns = aligner.align(grna_seq.upper(), tseq)
            if alns:
                best = alns[0]
                # Build full-length aligned strings
                grna_aln, genome_aln = _build_alignment(grna_seq.upper(), tseq, aligner)
                if not grna_aln:
                    pbar.update(1)
                    continue
                mm, disp, t_start, t_end = _score_alignment(grna_seq, pam, grna_aln, genome_aln)
                if mm <= max_mismatch:
                    results.append((record.id, best.score, mm, strand_label, genome_aln, grna_aln, disp, t_start, t_end))
            pbar.update(1)
    pbar.close()

    if not results:
        logger.warning("{}: no off-targets after SW", prefix)
        return None

    bed = outdir / f"{prefix}.parsing_water_for_visualization.offtarget.bed"

    # Build entries with read count, aligned sequence, ref seq
    bed_entries = []
    for rid, sc, mm, strand, genome_aln, grna_aln, disp, t_start, t_end in results:
        if ":" in rid:
            chrom = rid.split(":")[0]
            rest = rid.split(":", 1)[1]
            coords = rest.split("|")[0] if "|" in rest else rest
            if "-" in coords:
                start_str, end_str = coords.split("-")
                ext_start = int(start_str)
                ext_end = int(end_str)
            else:
                ext_start, ext_end = 0, 0
        else:
            chrom, ext_start, ext_end = rid, 0, 0
        # gRNA match genomic coordinates (within the slopped ext region), including PAM
        pam_len = len(pam)
        if strand == "+":
            start = ext_start + t_start
            end = ext_start + t_end + pam_len
        else:
            start = ext_end - t_end - pam_len
            end = ext_end - t_start
        # Parse count from rid header: ...|cnt=N|...
        read_count = 1
        if "|cnt=" in rid:
            cnt_part = rid.split("|cnt=")[1]
            cnt_val = cnt_part.split("|")[0] if "|" in cnt_part else cnt_part
            try:
                read_count = int(cnt_val)
            except ValueError:
                pass
        bed_entries.append((chrom, start, end, rid, sc, mm, strand, read_count, genome_aln, grna_aln, disp))

    # Sort by reads desc, then score desc
    bed_entries.sort(key=lambda x: (-x[7], -x[4]))

    # Move perfect match (all dots) to top
    perfect = [e for e in bed_entries if e[10].replace(".", "") == ""]
    # on-target(perfect)豁免 blacklist/control; 仅对 off-target(others)过滤
    others = [e for e in bed_entries
              if e[10].replace(".", "") != "" and not _in_excl(e[0], e[1], e[2])]

    # Merge all perfect (on-target) matches into a single reference row
    if perfect:
        ref = list(perfect[0])
        ref[7] = sum(e[7] for e in perfect)          # sum read counts
        ref[3] = ",".join(e[3] for e in perfect)     # join ids/coords
        ref[1] = min(e[1] for e in perfect)          # min start
        ref[2] = max(e[2] for e in perfect)          # max end
        bed_entries = [tuple(ref)] + others
    else:
        bed_entries = others

    # Deduplicate by position: keep higher score strand
    seen_pos = set()
    deduped = []
    for e in bed_entries:
        pos_key = (e[0], e[1], e[2])
        if pos_key not in seen_pos:
            seen_pos.add(pos_key)
            deduped.append(e)
    bed_entries = deduped

    ref_seq = grna_seq.upper()
    with open(bed, "w") as f:
        f.write("#chrom\tstart\tend\tid\treads\taligned_seq\tref_seq\tstrand\n")
        for chrom, start, end, rid, sc, mm, strand, read_count, genome_aln, grna_aln, disp in bed_entries:
            f.write(f"{chrom}\t{start}\t{end}\t{rid}\t{read_count}\t{disp}\t{ref_seq}\t{strand}\n")

    logger.info("{}: {} off-targets", prefix, len(results))
    return bed


def _score_alignment(grna: str, pam: str, grna_aln: str, genome_aln: str) -> tuple[int, str]:
    """Compute mismatch count (incl. internal indels) and display string.

    gRNA must align end-to-end (semi-global). Only the region covered by gRNA
    characters contributes to the score; target prefixes/suffixes are free.
    Display: '.' = match, base = mismatch, '-' = gap (insertion/deletion).

    The PAM bases (gRNA 3' flank) are appended to the display string (matched
    as '.', otherwise the actual base) but never counted in the mismatch score.
    """
    gstart = 0
    while gstart < len(grna_aln) and grna_aln[gstart] == "-":
        gstart += 1
    gend = len(grna_aln) - 1
    while gend >= 0 and grna_aln[gend] == "-":
        gend -= 1
    mm = 0
    disp = []
    ai = gstart
    for _ in range(len(grna)):
        while ai <= gend and grna_aln[ai] == "-":
            ai += 1
            mm += 1  # insertion into target
        if ai > gend:
            disp.append("-")
            mm += 1
            continue
        gc = grna_aln[ai]
        tc = genome_aln[ai]
        if tc == "-":
            disp.append("-")
            mm += 1  # deletion from target
        elif gc == tc:
            disp.append(".")
        else:
            disp.append(tc)
            mm += 1
        ai += 1

    # Append PAM bases (gRNA 3' flank) — not counted in mm
    if pam:
        for i in range(len(pam)):
            while ai < len(genome_aln) and genome_aln[ai] == "-":
                ai += 1
            if ai >= len(genome_aln):
                disp.append("-")
                continue
            tc = genome_aln[ai]
            gc = pam[i]
            if tc == "-":
                disp.append("-")
            elif gc == tc:
                disp.append(".")
            else:
                disp.append(tc)
            ai += 1

    # gRNA match interval within the target sequence
    t_start = sum(1 for ch in genome_aln[:gstart] if ch != "-")
    t_end = sum(1 for ch in genome_aln[:gend + 1] if ch != "-")
    return mm, "".join(disp), t_start, t_end


def _build_alignment(query: str, target: str, aligner) -> tuple[str, str]:
    """Build gapped alignment strings of equal length using coordinates.

    Unlike ``Alignment.aligned`` (whose blocks may not correspond character-wise
    when gaps split one sequence into many blocks), ``Alignment.coordinates``
    gives a per-column track: each column advances query and/or target, so we can
    reconstruct correctly aligned, equal-length strings.
    """
    alns = aligner.align(query, target)
    if not alns:
        return "", ""
    best = alns[0]
    coords = best.coordinates
    q_pos = coords[0]
    t_pos = coords[1]
    q_parts, t_parts = [], []
    for k in range(len(q_pos) - 1):
        qs, qe = q_pos[k], q_pos[k + 1]
        ts, te = t_pos[k], t_pos[k + 1]
        ql, tl = qe - qs, te - ts
        if ql == 0:
            q_parts.append("-" * tl)
            t_parts.append(target[ts:te])
        elif tl == 0:
            q_parts.append(query[qs:qe])
            t_parts.append("-" * ql)
        else:
            q_parts.append(query[qs:qe])
            t_parts.append(target[ts:te])
    return "".join(q_parts), "".join(t_parts)
