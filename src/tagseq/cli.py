"""CLI — click-based command line interface."""

from __future__ import annotations

import logging
import sys

import click

from .config import load_config
from .exceptions import TagseqError
from .pipeline import step_align, step_create_makefile, step_find_target
from .report import generate_report, print_report

logger = logging.getLogger("tagseq")


def _setup_logging(v: bool) -> None:
    logging.basicConfig(
        level=logging.DEBUG if v else logging.INFO,
        format="%(levelname)s %(message)s",
        stream=sys.stderr,
    )


@click.group()
@click.option("-c", "--config", required=True, help="Config file")
@click.option("-v", "--verbose", is_flag=True)
@click.option("--resume", is_flag=True, help="Skip completed steps")
@click.pass_context
def cli(ctx, config, verbose, resume):
    _setup_logging(verbose)
    ctx.ensure_object(dict)
    ctx.obj["cfg"] = load_config(config)
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
@click.pass_context
def all(ctx):
    """Full pipeline."""
    ctx.invoke(create_makefile)
    ctx.invoke(align)
    ctx.invoke(find_target)


def main():
    cli()


if __name__ == "__main__":
    main()
