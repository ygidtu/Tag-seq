"""CLI — click-based command line interface."""

from __future__ import annotations

from loguru import logger
import sys

import click

from .config import load_config, validate_config
from .exceptions import TagseqError
from .pipeline import step_align, step_create_makefile, step_find_target, _read_grnas
from .report import generate_report, print_report
from .plot import plot_statistics, plot_offtarget_alignment, plot_mismatch_profile
from .qc import generate_html



def _setup_logging(v: bool) -> None:
    logger.remove()
    level = "DEBUG" if v else "INFO"
    logger.add(sys.stderr, level=level)


@click.group()
@click.option("-c", "--config", required=True, help="Config file")
@click.option("-v", "--verbose", is_flag=True)
@click.option("--resume", is_flag=True, help="Skip completed steps")
@click.pass_context
def cli(ctx, config, verbose, resume):
    _setup_logging(verbose)
    ctx.ensure_object(dict)
    ctx.obj["cfg"] = load_config(config)
    issues = validate_config(ctx.obj["cfg"])
    for iss in issues:
        logger.error("Config: {}", iss)
    if issues:
        raise click.Abort()
    ctx.obj["resume"] = resume


@cli.command()
@click.pass_context
def create_makefile(ctx):
    """Prepare output directories."""
    step_create_makefile(ctx.obj["cfg"])


@cli.command()
@click.pass_context
def align(ctx):
    """Run QC → trimming → alignment → target detection."""
    step_align(ctx.obj["cfg"])


@cli.command(name="find-target")
@click.pass_context
def find_target(ctx):
    """Identify off-target sites per gRNA."""
    step_find_target(ctx.obj["cfg"])


@cli.command()
@click.pass_context
def report(ctx):
    """Show statistics only."""
    cfg = ctx.obj["cfg"]
    print_report(generate_report(cfg.prefix, cfg.outdir))


@cli.command()
@click.option("--min-reads", default=1, type=int, help="Minimum read count filter")
@click.pass_context
def all(ctx, min_reads):
    """Full pipeline + plot."""
    ctx.invoke(create_makefile)
    ctx.invoke(align)
    ctx.invoke(find_target)
    ctx.obj["min_reads"] = min_reads
    ctx.invoke(plot)


@cli.command()
@click.option("--min-reads", default=1, type=int, help="Minimum read count filter (default: 1)")
@click.pass_context
def plot(ctx, min_reads):
    """Generate plots from existing results (skip analysis)."""
    cfg = ctx.obj["cfg"]

    # Stats chart
    plot_statistics(cfg)

    # Mismatch profile
    plot_mismatch_profile(cfg)

    # QC report
    qc_path = generate_html(cfg)
    if qc_path:
        logger.info("QC report -> {}", qc_path)

    # Off-target alignment plots for each gRNA
    if cfg.grna_file and cfg.grna_file.exists():
        grnas = _read_grnas(cfg.grna_file)
        for gid, gseq, gpam in grnas:
            target_dir = cfg.find_target_dir(gid)
            offtarget_bed = target_dir / f"{gid}.parsing_water_for_visualization.offtarget.bed"
            ext_fa = target_dir / f"{gid}.ext.fa"
            if offtarget_bed.exists():
                out_pdf = target_dir / f"{gid}_offtargets.pdf"
                plot_offtarget_alignment(offtarget_bed, ext_fa, gseq, gpam or "NGG", out_pdf, title=gid, min_reads=min_reads)
            else:
                logger.warning("No off-target BED for {}, run 'find-target' first", gid)
    else:
        logger.warning("No gRNA file configured, skipping off-target plots")

    # QC report
    qc_path = generate_html(ctx.obj["cfg"])
    if qc_path:
        logger.info("QC report -> {}", qc_path)


@cli.command()
@click.pass_context
def qc(ctx):
    """Generate QC HTML report from existing results."""
    qc_path = generate_html(ctx.obj["cfg"])
    if qc_path:
        logger.info("QC report -> {}", qc_path)


def main():
    cli()


if __name__ == "__main__":
    main()
