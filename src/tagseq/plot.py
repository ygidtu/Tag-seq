"""Off-target alignment PDF/SVG plots — matching vis.py reference style."""

from __future__ import annotations
from pathlib import Path
from loguru import logger

BASE_COLORS = {
    "G": "#F5F500", "A": "#FF5454", "T": "#00D118", "C": "#26A8FF",
    "N": "#B3B3B3",
}

BS = 15  # box size matching vis.py


def plot_statistics(cfg) -> Path | None:
    try:
        import matplotlib
        matplotlib.use("pdf")
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available")
        return None
    from .report import generate_report
    reports = generate_report(cfg.prefix, cfg.outdir)
    if not reports:
        return None
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    colors = ["#6baed6", "#74c476", "#fd8d3c", "#9852a2", "#e6550d"]
    labels = ["Raw", "ODN+", "Trim", "Aligned", "Unique"]
    for idx, sr in enumerate(reports):
        ax = axes[idx]
        vals = [sr.total_reads, sr.odn_passed, sr.trim_passed,
                sr.alignment.total_aligned, sr.dedup_unique]
        ax.bar(labels, vals, color=colors, edgecolor="white", width=0.6)
        ax.set_title(sr.name, fontsize=13, fontweight="bold")
        ax.set_ylabel("Reads")
        ax.set_yscale("log")
        for b, v in zip(ax.patches, vals):
            ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.1,
                    f"{v:,}", ha="center", va="bottom", fontsize=8, rotation=45)
    plt.tight_layout()
    out = cfg.outdir / "statistics.pdf"
    plt.savefig(out, dpi=150)
    plt.close()
    logger.info("Statistics -> {}", out)
    return out


def plot_mismatch_profile(cfg) -> Path | None:
    """Generate mismatch position frequency bar chart."""
    try:
        import matplotlib
        matplotlib.use("pdf")
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    from Bio.Align import PairwiseAligner
    from .offtarget import _read_bed

    prefix = cfg.prefix
    grna_file = cfg.grna_file
    if not grna_file or not grna_file.exists():
        return None

    grnas = []
    with open(grna_file) as f:
        for line in f:
            p = line.strip().split()
            if len(p) >= 2:
                grnas.append((p[0], p[1].upper()))

    all_mm_positions = []

    for gid, grna_seq in grnas:
        target_dir = cfg.find_target_dir(gid)
        bed_file = target_dir / f"{gid}.parsing_water_for_visualization.offtarget.bed"
        ext_fa = target_dir / f"{gid}.ext.fa"
        if not bed_file.exists() or not ext_fa.exists():
            continue

        seq_map = {}
        with open(ext_fa) as f:
            lines = f.readlines()
        for i, line in enumerate(lines):
            if line.startswith(">"):
                key = line[1:].strip()
                seq = ""
                for j in range(i + 1, len(lines)):
                    if lines[j].startswith(">"): break
                    seq += lines[j].strip()
                seq_map[key] = seq.upper()

        aligner = PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 5
        aligner.mismatch_score = -4
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -0.5

        with open(bed_file) as f:
            for line in f:
                if line.startswith("#"): continue
                p = line.strip().split("\t")
                if len(p) < 5: continue
                tseq = seq_map.get(p[3], "")
                if not tseq: continue
                strand = p[6] if len(p) > 6 else "+"
                if strand == "-": tseq = str(tseq)[::-1]
                grna_aln, genome_aln = _build_alignment(grna_seq, tseq, aligner)
                if not grna_aln: continue
                ai = 0
                for ri in range(len(grna_seq)):
                    while ai < len(grna_aln) and grna_aln[ai] == "-":
                        ai += 1
                    if ai >= len(grna_aln): break
                    if ai < len(genome_aln) and genome_aln[ai] != "-" and grna_aln[ai] != genome_aln[ai]:
                        all_mm_positions.append(ri + 1)
                    ai += 1

    if not all_mm_positions:
        return None

    from collections import Counter
    pos_counts = Counter(all_mm_positions)
    positions = sorted(pos_counts.keys())
    counts = [pos_counts[p] for p in positions]
    total = sum(counts)

    fig, ax = plt.subplots(figsize=(max(8, len(grna_seq) * 0.4), 4))
    ax.bar(positions, counts, color="#e6550d", edgecolor="white", width=0.7)
    ax.set_xlabel("gRNA Position")
    ax.set_ylabel("Mismatch Count")
    ax.set_title(f"Mismatch Position Distribution (n={total})")
    ax.set_xticks(range(1, len(grna_seq) + 1))
    for p, c in zip(positions, counts):
        ax.text(p, c + 0.3, str(c), ha="center", va="bottom", fontsize=8)
    plt.tight_layout()
    out = cfg.outdir / "mismatch_profile.pdf"
    plt.savefig(out, dpi=150)
    plt.close()
    logger.info("Mismatch profile -> {}", out)
    return out


def _build_alignment(query: str, target: str, aligner) -> tuple[str, str]:
    """Build gapped alignment strings of equal length."""
    alns = aligner.align(query, target)
    if not alns:
        return "", ""
    best = alns[0]
    coords = best.aligned
    q_blocks = coords[0]
    t_blocks = coords[1]
    q_parts, t_parts = [], []
    qi = ti = 0
    for (qs, qe), (ts, te) in zip(q_blocks, t_blocks):
        if qs > qi:
            q_parts.append(query[qi:qs])
            t_parts.append("-" * (qs - qi))
        if ts > ti:
            q_parts.append("-" * (ts - ti))
            t_parts.append(target[ti:ts])
        q_parts.append(query[qs:qe])
        t_parts.append(target[ts:te])
        qi = qe
        ti = te
    if qi < len(query):
        q_parts.append(query[qi:])
        t_parts.append("-" * (len(query) - qi))
    if ti < len(target):
        q_parts.append("-" * (len(target) - ti))
        t_parts.append(target[ti:])
    return "".join(q_parts), "".join(t_parts)


def plot_offtarget_alignment(
    offtarget_bed: Path, ext_fa: Path,
    grna_seq: str, pam: str,
    out_pdf: Path, title: str = "",
    max_rows: int = 50,
    min_reads: int = 1,
) -> Path | None:
    import matplotlib
    matplotlib.use("pdf")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.font_manager import FontProperties
    from Bio.Align import PairwiseAligner

    if not offtarget_bed.exists() or offtarget_bed.stat().st_size == 0:
        logger.warning("Off-target BED missing: {}", offtarget_bed)
        return None

    # Read BED (format: chrom start end id reads aligned_seq ref_seq strand)
    all_rows = []
    with open(offtarget_bed) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.strip().split("\t")
            if len(p) < 8:
                continue
            chrom = p[0]
            start = int(p[1])
            end = int(p[2])
            read_count = int(p[4])
            aligned_seq = p[5]
            ref_seq = p[6]
            strand = p[7]
            sw_score = float(p[4])  # use reads as score for reference detection
            pos_key = (chrom, start, end)
            all_rows.append({
                "id": p[3], "score": read_count,
                "mm": 0, "strand": strand,
                "aligned_seq": aligned_seq,
                "ref_seq": ref_seq,
                "coord": f"{chrom}:{start}-{end}",
                "pos_key": pos_key,
            })

    if not all_rows:
        logger.warning("No off-target data")
        return None

    # Reference = site with SW score=100 (perfect match) — using reads as proxy
    # Perfect match = aligned_seq has no mismatches (all dots, no gaps)
    perfect = [r for r in all_rows if r["aligned_seq"].replace(".", "") == ""]
    if perfect:
        ref_site = perfect[0]
        ref_count = int(ref_site["score"])
        ref_coord = ref_site["coord"]
        seen = {ref_site["pos_key"]}
        rows = [ref_site]
    else:
        ref_count = 0
        ref_coord = ""
        seen = set()
        rows = []

    # Deduplicate by position: keep first occurrence (BED order)
    for r in all_rows:
        if r["pos_key"] not in seen:
            seen.add(r["pos_key"])
            rows.append(r)

    # Apply min_reads filter (keep reference regardless if exists)
    if perfect:
        rows = [rows[0]] + [r for r in rows[1:] if r["score"] >= min_reads]
    else:
        rows = [r for r in rows if r["score"] >= min_reads]
    rows = rows[:max_rows]

    # Align
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 5
    aligner.mismatch_score = -4
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    ref = grna_seq.upper()
    ref_len = len(ref)
    n = len(rows)
    font = FontProperties(family="monospace", size=11)
    font_small = FontProperties(family="monospace", size=8)
    font_title = FontProperties(family="monospace", size=13, weight="bold")

    # Layout — y increases downward
    margin_top = 50
    x0 = 20
    y_title = 20
    y_ticks = 38
    y_ref = 48
    row_h = BS + 3
    right_col = x0 + (ref_len + 3) * BS

    fig_w = (ref_len + 10) * BS / 72 + 1.5
    total_h = margin_top + 20 + (n + 2) * row_h
    fig_h = total_h / 72 + 0.5

    fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_h))
    ax.set_xlim(0, (ref_len + 10) * BS + 40)
    ax.set_ylim(total_h, 0)
    ax.set_aspect("equal")
    ax.axis("off")

    def center_x(i):
        return x0 + i * BS + BS / 2

    def center_y(y):
        return y + BS / 2

    # Title
    if title:
        ax.text(x0, y_title, title, fontproperties=font_title)

    # Position ticks (every 10th position)
    for pos in range(1, ref_len + 1):
        if pos == 1 or pos == ref_len or pos % 10 == 0:
            ax.text(center_x(pos - 1), y_ticks, str(pos),
                    ha="center", va="center", fontproperties=font_small)

    # PAM labels
    pam_len = len(pam)
    pam_start = ref_len - pam_len
    for pi in range(pam_len):
        ax.text(center_x(pam_start + pi), y_ticks - 12, pam[pi],
                ha="center", va="center",
                fontproperties=FontProperties(family="monospace", size=8, weight="bold"),
                color="#555")

    # Reference row: colored boxes for each base
    y = y_ref
    for i in range(ref_len):
        c = ref[i]
        color = "#D0D0D0" if i >= pam_start else BASE_COLORS.get(c, "#B3B3B3")
        ax.add_patch(mpatches.FancyBboxPatch(
            (x0 + i * BS, y), BS, BS, boxstyle="round,pad=0", facecolor=color, edgecolor="#999"))
        ax.text(center_x(i), center_y(y), c, ha="center", va="center", fontproperties=font, color="black")

    # Build aligned entries from pre-computed BED data (no re-alignment needed)
    ref_seq = grna_seq.upper() if not perfect else perfect[0]["ref_seq"]
    aligned = []
    for r in rows:
        aligned_seq = r["aligned_seq"]
        genome_aln = "".join(r["aligned_seq"])
        coord = r["coord"]
        aligned.append((r["id"], r["score"], r["mm"], r["strand"],
                        aligned_seq, genome_aln, coord))
    if not aligned:
        logger.warning("No alignments for {}", title)
        return None

    # Separate reference from off-targets
    if perfect:
        ref_entry = aligned[0]
        off_entries = aligned[1:]
    else:
        ref_entry = None
        off_entries = aligned

    total_off_reads = sum(sc for _, sc, _, _, _, _, _ in off_entries)
    ax.text(right_col, center_y(y_ref), "Sites", va="center",
            fontproperties=FontProperties(family="monospace", size=9, weight="bold"), color="#333")

    # Aligned rows — one column per reference base
    for j, (sid, score, mm, strand, aligned_seq, genome_aln, coord) in enumerate(off_entries):
        y = y_ref + (j + 1) * row_h

        for ri, ch in enumerate(aligned_seq):
            if ri >= ref_len:
                break
            x = x0 + ri * BS
            cx, cy = x + BS / 2, y + BS / 2
            if ch == "-" or ch == ".":
                ax.add_patch(mpatches.Circle((cx, cy), BS * 0.12, facecolor="#333", edgecolor="none"))
            else:
                color = BASE_COLORS.get(ch, "#B3B3B3")
                ax.add_patch(mpatches.FancyBboxPatch((x, y), BS, BS, boxstyle="round,pad=0", facecolor=color, edgecolor="#999"))
                ax.text(cx, cy, ch, ha="center", va="center", fontproperties=font, color="black")

        # Right side: count (pct) + coord
        pct = score / total_off_reads * 100 if total_off_reads > 0 else 0
        label = f"{int(score)} ({pct:.2f}%)  {coord}"
        ax.text(right_col, center_y(y), label, va="center", fontproperties=font_small, color="#333")

    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    plt.savefig(str(out_pdf), dpi=200, format="pdf", bbox_inches="tight")
    out_svg = out_pdf.with_suffix(".svg")
    plt.savefig(str(out_svg), dpi=200, format="svg", bbox_inches="tight")
    plt.close()
    logger.info("Alignment -> {} (PDF+SVG, {} sites)", out_pdf, n)
    return out_pdf
