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

    # 正负文库使用相同的 R1/R2，ODN 去除只需读取一次
    first_name, _ = samples[0]
    dd0 = cfg.data_dir(first_name)
    dd0.mkdir(parents=True, exist_ok=True)
    r1_odn, r2_odn, _, _ = remove_odn(
        cfg.lib_r1, cfg.lib_r2, cfg.forward_tag, cfg.reverse_tag, dd0, first_name)

    for sample_name, tag in samples:
        logger.info("=== {} ===", sample_name)
        dd = cfg.data_dir(sample_name)
        dd.mkdir(parents=True, exist_ok=True)

        # 1. ODN removal 结果（正负文库复用同一份）
        r1 = dd / f"{sample_name}.rmODN.R1.fq.gz"
        r2 = dd / f"{sample_name}.rmODN.R2.fq.gz"
        if r1 != r1_odn or r2 != r2_odn:
            shutil.copy(r1_odn, r1)
            shutil.copy(r2_odn, r2)

        # 2. Cutadapt（修剪 3' 接头，保留 5' primer/UMI 供 UMI 提取）
        _run_cutadapt(r1, r2, cfg, sample_name)
        r1t = dd / f"{sample_name}.Trim.R1.fq.gz"
        r2t = dd / f"{sample_name}.Trim.R2.fq.gz"

        # 3. UMI 提取（改造后文库 UMI 位于 primer 之后）
        r1u, r2u = extract_umi(r1t, r2t, dd, sample_name,
                               umi_len=cfg.umi_len, umi_offset=cfg.umi_offset,
                               umi_prefix=cfg.umi_prefix)

        # 4. STAR
        bam = run_star(r1u, r2u, cfg.index, cfg.outdir, sample_name, cfg.threads)

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
