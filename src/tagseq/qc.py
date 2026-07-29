"""QC statistics collection and HTML report generation."""

from __future__ import annotations
from pathlib import Path
from dataclasses import dataclass, field
from .config import reverse_complement


@dataclass
class SampleQC:
    name: str = ""
    # Raw
    raw_reads: int = 0
    # ODN
    odn_input: int = 0
    odn_passed: int = 0
    # Trim
    trim_input: int = 0
    trim_passed: int = 0
    # UMI
    umi_input: int = 0
    # STAR
    star_input: int = 0
    star_unique: int = 0
    star_multi: int = 0
    star_too_many: int = 0
    star_unmapped: int = 0
    # Dedup
    dedup_pre: int = 0
    dedup_post: int = 0
    # Targets
    clusters: int = 0
    plus_reads: int = 0
    minus_reads: int = 0


def collect_qc(cfg) -> list[SampleQC]:
    """Collect QC stats from pipeline output files."""
    from .alignment import parse_star_log
    import re

    reports = []
    for sample_name, tag in cfg.samples:
        qc = SampleQC(name=sample_name)
        dd = cfg.data_dir(sample_name)
        ad = cfg.align_dir(sample_name)
        td = cfg.target_dir(sample_name)

        # ODN stat gives us total reads (avoids counting huge FASTQ)
        stat_file = dd / f"{sample_name}.rmODN.stat"
        if stat_file.exists():
            with open(stat_file) as f:
                for line in f:
                    m = re.search(r"Raw flagment count:\s+(\d+)", line)
                    if m: qc.odn_input = int(m.group(1))
                    qc.raw_reads = qc.odn_input
                    m = re.search(r"Read count with ODN:\s+(\d+)", line)
                    if m: qc.odn_passed = int(m.group(1))

        # Trim (cutadapt report)
        trim_rpt = dd / f"{sample_name}.trim_report.txt"
        if trim_rpt.exists():
            with open(trim_rpt) as f:
                text = f.read()
            m = re.search(r"Total read pairs processed:\s+([\d,]+)", text)
            if m: qc.trim_input = int(m.group(1).replace(",", ""))
            m = re.search(r"Pairs written \(passing filters\):\s+([\d,]+)", text)
            if m: qc.trim_passed = int(m.group(1).replace(",", ""))

        # STAR
        star_log = ad / f"{sample_name}.Log.final.out"
        stats = parse_star_log(star_log)
        qc.star_input = stats.input_reads
        qc.star_unique = stats.unique_mapped
        qc.star_multi = stats.multi_mapped
        qc.star_too_many = stats.too_many_loci
        qc.star_unmapped = max(0, stats.input_reads - stats.total_aligned)

        # Dedup — count mapped reads in both BAMs using same method
        pre_bam = ad / f"{sample_name}.Aligned.sortedByCoord.out.bam"
        dedup_bam = ad / f"{sample_name}.Aligned.sortedByCoord.out.dedup.bam"
        import pysam
        if pre_bam.exists():
            with pysam.AlignmentFile(str(pre_bam), "rb") as bam:
                qc.dedup_pre = sum(1 for r in bam if not r.is_unmapped)
        if dedup_bam.exists():
            with pysam.AlignmentFile(str(dedup_bam), "rb") as bam:
                qc.dedup_post = sum(1 for r in bam if not r.is_unmapped)

        # Targets
        plus_bed = td / f"{sample_name}.plus.bed"
        minus_bed = td / f"{sample_name}.minus.bed"
        if plus_bed.exists():
            with open(plus_bed) as f:
                qc.plus_reads = sum(1 for _ in f)
        if minus_bed.exists():
            with open(minus_bed) as f:
                qc.minus_reads = sum(1 for _ in f)

        reports.append(qc)
    return reports


def _pct(a: int, b: int) -> float:
    return a / b * 100 if b > 0 else 0.0


def generate_html(cfg) -> Path | None:
    """Generate QC HTML report."""
    reports = collect_qc(cfg)

    if not reports:
        return None

    rows_html = ""
    for qc in reports:
        odn_pct = _pct(qc.odn_passed, qc.odn_input)
        trim_pct = _pct(qc.trim_passed, qc.trim_input)
        align_pct = _pct(qc.star_unique, qc.star_input)
        dedup_pct = _pct(qc.dedup_post, qc.dedup_pre) if qc.dedup_pre > 0 else 0
        reads_kept = _pct(qc.dedup_post, qc.raw_reads)

        def bar(pct, color="#4daf4a"):
            w = max(pct, 1)
            return f'<div style="background:#eee;border-radius:4px;width:200px;height:18px;display:inline-block;vertical-align:middle;margin-right:8px"><div style="background:{color};width:{w:.1f}%;height:18px;border-radius:4px"></div></div>'

        rows_html += f"""
<tr style="border-bottom:2px solid #ddd">
  <td colspan="2" style="padding:12px 8px;font-weight:bold;font-size:15px;background:#f8f9fa">{qc.name}</td>
</tr>
<tr>
  <td style="padding:4px 8px 4px 24px;color:#555">Raw FASTQ</td>
  <td style="padding:4px 8px;text-align:right">{qc.raw_reads:,}</td>
  <td style="padding:4px 8px">{bar(100)}100%</td>
  <td style="padding:4px 8px;color:#999">—</td>
</tr>
<tr>
  <td style="padding:4px 8px 4px 24px;color:#555">ODN filter</td>
  <td style="padding:4px 8px;text-align:right">{qc.odn_passed:,} / {qc.odn_input:,}</td>
  <td style="padding:4px 8px">{bar(odn_pct, "#377eb8")}{odn_pct:.2f}%</td>
  <td style="padding:4px 8px">Tag in R2 (exact match)</td>
</tr>
<tr>
  <td style="padding:4px 8px 4px 24px;color:#555">Adapter trim</td>
  <td style="padding:4px 8px;text-align:right">{qc.trim_passed:,} / {qc.trim_input:,}</td>
  <td style="padding:4px 8px">{bar(trim_pct, "#ff7f00")}{trim_pct:.2f}%</td>
  <td style="padding:4px 8px">cutadapt, minlen={cfg.minlen}</td>
</tr>
<tr>
  <td style="padding:4px 8px 4px 24px;color:#555">STAR alignment</td>
  <td style="padding:4px 8px;text-align:right">{qc.star_unique:,} / {qc.star_input:,}</td>
  <td style="padding:4px 8px">{bar(align_pct, "#e41a1c")}{align_pct:.2f}%</td>
  <td style="padding:4px 8px">Unique mapped</td>
</tr>
<tr>
  <td style="padding:4px 8px 4px 24px;color:#555">UMI dedup</td>
  <td style="padding:4px 8px;text-align:right">{qc.dedup_post:,} / {qc.dedup_pre:,}</td>
  <td style="padding:4px 8px">{bar(dedup_pct, "#984ea3")}{dedup_pct:.2f}%</td>
  <td style="padding:4px 8px">Mapped reads kept after UMI dedup</td>
</tr>
<tr style="border-top:1px dashed #ccc">
  <td style="padding:6px 8px 6px 24px;font-weight:bold;color:#333">Reads kept (end-to-end)</td>
  <td style="padding:6px 8px;text-align:right;font-weight:bold">{qc.dedup_post:,} / {qc.raw_reads:,}</td>
  <td style="padding:6px 8px">{bar(reads_kept, "#e6550d")}{reads_kept:.4f}%</td>
  <td style="padding:6px 8px;color:#999">Raw → Unique</td>
</tr>
<tr>
  <td style="padding:4px 8px 4px 24px;color:#555">Clusters (plus)</td>
  <td style="padding:4px 8px;text-align:right">{qc.plus_reads:,}</td>
  <td colspan="2" style="padding:4px 8px"></td>
</tr>
<tr>
  <td style="padding:4px 8px 4px 24px;color:#555">Clusters (minus)</td>
  <td style="padding:4px 8px;text-align:right">{qc.minus_reads:,}</td>
  <td colspan="2" style="padding:4px 8px"></td>
</tr>"""

    grna_section = ""
    if cfg.grna_file and cfg.grna_file.exists():
        from .pipeline import _read_grnas
        grnas = _read_grnas(cfg.grna_file)
        grna_rows = ""
        for gid, gseq, gpam in grnas:
            tdir = cfg.find_target_dir(gid)
            bed_file = tdir / f"{gid}.parsing_water_for_visualization.offtarget.bed"
            n_off = 0
            n_perfect = 0
            if bed_file.exists():
                with open(bed_file) as f:
                    for line in f:
                        if line.startswith("#"): continue
                        n_off += 1
                        p = line.strip().split("\t")
                        if len(p) >= 6 and p[5].replace(".", "") == "":
                            n_perfect += 1
            grna_rows += f"<tr><td>{gid}</td><td>{gseq}</td><td>{gpam}</td><td style='text-align:right'>{n_off}</td><td style='text-align:right'>{n_perfect}</td></tr>\n"
        if grna_rows:
            grna_section = f"""
<h3 style="color:#333;margin-top:24px">gRNA Off-target Summary</h3>
<table style="width:100%;border-collapse:collapse;font-size:13px">
<tr style="background:#f8f9fa;font-weight:bold"><td style="padding:8px">gRNA</td><td style="padding:8px">Sequence</td><td style="padding:8px">PAM</td><td style="padding:8px;text-align:right">Sites</td><td style="padding:8px;text-align:right">Perfect</td></tr>
{grna_rows}</table>"""

    html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>Tag-seq QC Report</title></head>
<body style="font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;max-width:960px;margin:20px auto;padding:0 20px;color:#333;line-height:1.5">

<h1 style="border-bottom:2px solid #4daf4a;padding-bottom:8px">Tag-seq QC Report</h1>
<p style="color:#666;font-size:14px">Sample: <strong>{cfg.prefix}</strong> | Genome: {cfg.genome} | Date: generated on-the-fly</p>

<h2 style="color:#333;margin-top:24px">Pipeline Statistics</h2>
<table style="width:100%;border-collapse:collapse;font-size:13px">
<tr style="background:#4daf4a;color:white">
  <th style="padding:10px 8px;text-align:left;width:28%">Step</th>
  <th style="padding:10px 8px;text-align:right;width:22%">Count</th>
  <th style="padding:10px 8px;text-align:left;width:30%">Survival</th>
  <th style="padding:10px 8px;text-align:left;width:20%">Note</th>
</tr>
{rows_html}
</table>

{grna_section}

<h2 style="color:#333;margin-top:24px">Configuration</h2>
<table style="width:100%;border-collapse:collapse;font-size:12px;color:#555">
<tr><td style="padding:4px 8px">MinSupportReadCount</td><td style="padding:4px 8px">{cfg.min_support_readcount}</td></tr>
<tr><td style="padding:4px 8px">MinCuttingEventCount</td><td style="padding:4px 8px">{cfg.min_cutting_event_count}</td></tr>
<tr><td style="padding:4px 8px">MaxMismatch</td><td style="padding:4px 8px">{cfg.max_mismatch}</td></tr>
<tr><td style="padding:4px 8px">MINLEN</td><td style="padding:4px 8px">{cfg.minlen}</td></tr>
<tr><td style="padding:4px 8px">THREAD</td><td style="padding:4px 8px">{cfg.threads}</td></tr>
</table>

<p style="color:#999;font-size:11px;margin-top:30px;border-top:1px solid #eee;padding-top:10px">
Tag-seq Pipeline v2 | Generated by tagseq.qc module</p>
</body></html>"""

    out = cfg.outdir / "qc_report.html"
    with open(out, "w") as f:
        f.write(html)
    return out
