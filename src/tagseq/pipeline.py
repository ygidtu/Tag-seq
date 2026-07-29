"""Pipeline orchestrator — one function per step, no duplication."""

from __future__ import annotations

import logging
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

logger = logging.getLogger(__name__)

# ── helpers ──


def _check_tools() -> None:
    """Verify external binaries exist. Only STAR is strictly required."""
    for cmd in ["STAR", "cutadapt"]:
        if not subprocess.run(["which", cmd], capture_output=True).returncode == 0:
            raise ExternalToolError(f"Required tool not found: {cmd}")


def _run_cutadapt(r1: Path, r2: Path, cfg: PipelineConfig, sample: str) -> None:
    """Trim adapters with cutadapt."""
    dd = cfg.data_dir(sample)
    r1o = dd / f"{sample}.Trim.R1.fq"
    r2o = dd / f"{sample}.Trim.R2.fq"
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
                p = line.strip().split("\t")
                if len(p) >= 2:
                    cmd += ["-a", p[0], "-A", p[1]]
    log = dd / f"{sample}.trim_report.txt"
    r = subprocess.run(cmd, stdout=open(log, "w"), stderr=subprocess.STDOUT)
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
    logger.info("Directories ready under %s", cfg.outdir)


def step_align(cfg: PipelineConfig) -> None:
    _check_tools()
    for sample_name, tag in cfg.samples:
        logger.info("=== %s ===", sample_name)
        dd = cfg.data_dir(sample_name)
        dd.mkdir(parents=True, exist_ok=True)

        # 1. ODN removal
        r1, r2, _ = remove_odn(cfg.lib_r1, cfg.lib_r2, tag, dd, sample_name)

        # 2. Cutadapt
        _run_cutadapt(r1, r2, cfg, sample_name)
        r1t = dd / f"{sample_name}.Trim.R1.fq"
        r2t = dd / f"{sample_name}.Trim.R2.fq"

        # 3. UMI extraction
        r1u, r2u = extract_umi(r1t, r2t, dd, sample_name)

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
    logger.info("Loaded %d gRNAs from %s", len(grnas), cfg.grna_file)

    for gid, gseq, gpam in grnas:
        tdir = cfg.find_target_dir(gid)
        dm = tdir / f"{gid}.parsing_water_for_visualization.offtarget.bed"
        if dm.exists():
            logger.info("%s: already done", gid)
            continue

        def pf(name):
            return cfg.outdir / f"{name}/02potentialTargets/{name}"

        pp = pf(f"{cfg.prefix}_plus") / "plus.proximal"
        pm = pf(f"{cfg.prefix}_plus") / "minus.proximal"
        mp = pf(f"{cfg.prefix}_minus") / "plus.proximal"
        mm = pf(f"{cfg.prefix}_minus") / "minus.proximal"

        logger.info("Off-target: %s (PAM=%s)", gid, gpam)
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
