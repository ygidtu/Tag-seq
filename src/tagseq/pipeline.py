"""Pipeline orchestrator — one function per step, no duplication."""

from __future__ import annotations

from loguru import logger
import subprocess
from pathlib import Path

from .alignment import run_star
from .config import warn_tag_orientation
from .exceptions import ExternalToolError
from .fastq import detect_phred
from .offtarget import find_offtargets
from .report import generate_report, print_report
from .targets import detect_targets
from .trim import remove_odn
from .types import PipelineConfig
from .umi import dedup_umi, extract_umi


# ── helpers ──


def _check_tools() -> None:
    """Verify external binaries exist. Only STAR is strictly required."""
    for cmd in ["STAR", "cutadapt"]:
        if not subprocess.run(["which", cmd], capture_output=True).returncode == 0:
            raise ExternalToolError(f"Required tool not found: {cmd}")


def proximal_path(sample_dir: Path, sample_name: str, strand: str) -> Path:
    return sample_dir / f"{sample_name}.{strand}.proximal"


def _run_cutadapt(r1: Path, r2: Path, cfg: PipelineConfig, sample: str) -> None:
    """Trim adapters with cutadapt (supports 3' and 5' adapters)."""
    dd = cfg.data_dir(sample)
    r1o = dd / f"{sample}.Trim.R1.fq.gz"
    r2o = dd / f"{sample}.Trim.R2.fq.gz"

    if not r1.exists() or not r2.exists():
        raise ExternalToolError(f"Input files for cutadapt not found: {r1}, {r2}")

    cmd = [
        "cutadapt",
        "--quality-cutoff", "5",
        "--minimum-length", str(cfg.minlen),
        "--trim-n",
        "-j", str(cfg.threads),
        "-o", str(r1o), "-p", str(r2o),
        str(r1), str(r2),
    ]
    if cfg.adapter and cfg.adapter.exists():
        with open(cfg.adapter) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                p = line.split("\t")
                if len(p) < 2:
                    continue
                if p[0].startswith("5p:") or p[0].startswith("5P:"):
                    cmd += ["-g", p[0].split(":", 1)[1], "-G", p[1]]
                else:
                    cmd += ["-a", p[0], "-A", p[1]]
    log = dd / f"{sample}.trim_report.txt"
    r = subprocess.run(cmd, stdout=open(log, "w"), stderr=subprocess.STDOUT, timeout=600)
    if r.returncode != 0:
        raise ExternalToolError(f"cutadapt failed (exit={r.returncode})")


# ── step implementations ──

def parse_mixed(r1_path: Path, r2_path: Path, outdir: Path, prefix: str, cfg: PipelineConfig) -> tuple[Path, Path]:
    """混测文库解析: 逐对 R1/R2, 按 L/R 判定提取 UMI/barcode, 剥离两端人工序列;
    无法判定(ambiguous)的 read 丢弃。UMI 与文库类型写入 read 名供后续 dedup/追溯。"""
    from collections import Counter
    from .library import classify_extract, parse_r2
    from .fastq import open_fastq, open_fastq_write
    spec = cfg.lib_spec
    r1_out = outdir / f"{prefix}.mixed.R1.fq.gz"
    r2_out = outdir / f"{prefix}.mixed.R2.fq.gz"
    st = Counter()
    with open_fastq(r1_path) as f1, open_fastq(r2_path) as f2, \
         open_fastq_write(r1_out) as o1, open_fastq_write(r2_out) as o2:
        while True:
            h1 = f1.readline()
            if not h1: break
            s1 = f1.readline().rstrip("\n"); p1 = f1.readline(); q1 = f1.readline().rstrip("\n")
            h2 = f2.readline(); s2 = f2.readline().rstrip("\n"); p2 = f2.readline(); q2 = f2.readline().rstrip("\n")
            st["total"] += 1
            ce = classify_extract(s1, spec, cfg.forward_tag, cfg.reverse_tag)
            if ce is None:
                st["ambiguous"] += 1
                continue
            r2p = parse_r2(s2, spec)
            import re as _re
            base = _re.sub(r"/[12]$", "", h1.rstrip("\n").split()[0])
            hdr = f"{base}:UMI_{ce['umi']}:LIB_{ce['type']}"
            o1.write(f"{hdr}\n{s1[ce['g_start']:ce['g_end']]}\n+\n{q1[ce['g_start']:ce['g_end']]}\n")
            o2.write(f"{hdr}\n{s2[r2p['g_start']:r2p['g_end']]}\n+\n{q2[r2p['g_start']:r2p['g_end']]}\n")
            st[ce["type"]] += 1
    logger.info("parse_mixed: total={total} L={L} R={R} ambiguous={ambiguous} (kept {kept})",
                total=st["total"], L=st["L"], R=st["R"], ambiguous=st["ambiguous"],
                kept=st["L"] + st["R"])
    return r1_out, r2_out


def step_create_makefile(cfg: PipelineConfig) -> None:
    cfg.outdir.mkdir(parents=True, exist_ok=True)
    warns = []
    for sample_name, tag in cfg.samples:
        warn_tag_orientation(cfg.lib_r2, tag)
        d = cfg.data_dir(sample_name)
        d.mkdir(parents=True, exist_ok=True)
        if not cfg.lib_r1.exists():
            warns.append(f"R1 not found: {cfg.lib_r1}")
        if not cfg.lib_r2.exists():
            warns.append(f"R2 not found: {cfg.lib_r2}")
    for w in warns:
        logger.warning(w)
    logger.info("Directories ready under {}", cfg.outdir)


def step_align(cfg: PipelineConfig) -> None:
    _check_tools()
    import shutil
    samples = cfg.samples

    # 比对入口口径由 ODN_MODE 决定 (默认 rmODN = 与 TAG-1/3/5 既有结果一致)
    first_name, _ = samples[0]
    dd0 = cfg.data_dir(first_name)
    dd0.mkdir(parents=True, exist_ok=True)

    if cfg.odn_mode == "mixed":
        logger.info("ODN_MODE=mixed: parse_mixed 结构解析 (5' primer 锚点判 L/R, 丢弃 ambiguous) — 实验口径")
        pre = "mixed"
        r1p, r2p = parse_mixed(cfg.lib_r1, cfg.lib_r2, dd0, cfg.prefix, cfg)
    else:
        logger.info("ODN_MODE=rmODN: remove_odn 过滤 (R1/R2 任一端含 tag 才保留) — 既有 TAG-1/3/5 口径")
        pre = "rmODN"
        r1p, r2p, _, _ = remove_odn(
            cfg.lib_r1, cfg.lib_r2, cfg.forward_tag, cfg.reverse_tag, dd0, first_name)

    for sample_name, tag in samples:
        logger.info("=== {} ===", sample_name)
        dd = cfg.data_dir(sample_name)
        dd.mkdir(parents=True, exist_ok=True)

        # 1. 过滤/解析结果复用（正负文库同一份 clean reads）
        r1 = dd / f"{sample_name}.{pre}.R1.fq.gz"
        r2 = dd / f"{sample_name}.{pre}.R2.fq.gz"
        if r1 != r1p or r2 != r2p:
            shutil.copy(r1p, r1)
            shutil.copy(r2p, r2)

        # 2. Cutadapt（修剪残留 3' 接头）
        _run_cutadapt(r1, r2, cfg, sample_name)
        r1t = dd / f"{sample_name}.Trim.R1.fq.gz"
        r2t = dd / f"{sample_name}.Trim.R2.fq.gz"

        # 3. UMI: rmODN 口径下在此按 offset 提取 (mixed 口径已在 read 名中带 UMI)
        if cfg.odn_mode != "mixed":
            r1t, r2t = extract_umi(r1t, r2t, dd, sample_name,
                                   umi_len=cfg.umi_len, umi_offset=cfg.umi_offset,
                                   umi_prefix=cfg.umi_prefix)

        # 4. STAR
        bam = run_star(r1t, r2t, cfg.index, cfg.outdir, sample_name, cfg.threads)

        # 5. Dedup
        dedup_bam = dedup_umi(bam)

        # 6. Targets
        detect_targets(dedup_bam, cfg.outdir / sample_name, sample_name, cfg.chromsize)

    # Report
    reports = generate_report(cfg.prefix, cfg.outdir)
    print_report(reports)


def step_find_target(cfg: PipelineConfig) -> None:
    if not cfg.grna_file:
        logger.error("No gRNA file configured")
        return
    grnas = _read_grnas(cfg.grna_file)
    logger.info("Loaded {} gRNAs from {}", len(grnas), cfg.grna_file)

    for gid, gseq, gpam in grnas:
        tdir = cfg.find_target_dir(gid)
        dm = tdir / f"{gid}.parsing_water_for_visualization.offtarget.bed"
        if dm.exists():
            logger.info("{}: already done", gid)
            continue

        sample_dir_plus = cfg.outdir / f"{cfg.prefix}_plus" / "02potentialTargets"
        sample_dir_minus = cfg.outdir / f"{cfg.prefix}_minus" / "02potentialTargets"
        pname = f"{cfg.prefix}_plus"
        mname = f"{cfg.prefix}_minus"

        pp = proximal_path(sample_dir_plus, pname, "plus")
        pm = proximal_path(sample_dir_plus, pname, "minus")
        mp = proximal_path(sample_dir_minus, mname, "plus")
        mm = proximal_path(sample_dir_minus, mname, "minus")

        logger.info("Off-target: {} (PAM={})", gid, gpam)
        find_offtargets(
            pp, pm, mp, mm,
            control=cfg.ctrl,
            blacklist=cfg.blacklist,
            grna_seq=gseq,
            ref_fa=cfg.ref,
            chromsize=cfg.chromsize,
            prefix=gid,
            outdir=tdir,
            min_support=cfg.min_support_readcount,
            cut_events=cfg.min_cutting_event_count,
            max_mismatch=cfg.max_mismatch,
            pam=gpam,
        )


def _read_grnas(path: Path) -> list[tuple[str, str, str]]:
    out = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            p = line.split()
            if len(p) >= 3:
                out.append((p[0], p[1], p[2]))
    return out
