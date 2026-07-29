"""Config file parser — supports new (LIB_R1/LIB_R2) and legacy format."""

from __future__ import annotations

import gzip
import warnings
from pathlib import Path
from typing import Optional

from .exceptions import ConfigError
from .types import PipelineConfig

# ── helpers ──

def reverse_complement(seq: str) -> str:
    table = str.maketrans("ATCGatcg", "TAGCtagc")
    return seq.translate(table)[::-1]


def _resolve(path: str) -> Path:
    p = Path(path)
    return p if p.is_absolute() else (Path.cwd() / p).resolve()


def _read_optional(raw: dict, key: str, default: str = "") -> str:
    return raw.get(key, raw.get(key.upper(), raw.get(key.lower(), default)))


def _read_path(raw: dict, key: str) -> Optional[Path]:
    val = _read_optional(raw, key)
    if not val or val.lower() == "none":
        return None
    return _resolve(val)


def _read_int(raw: dict, key: str, default: int = 0) -> int:
    try:
        return int(_read_optional(raw, key, str(default)))
    except (ValueError, TypeError):
        return default


# ── backward-compatible key aliases ──
_OLD_KEYS = {
    "FORWARD_LIB_TAG": "FORWARD_TAG",
    "REVERSE_LIB_TAG": "REVERSE_TAG",
}


def _normalise(raw: dict) -> dict:
    """Map legacy keys to new names."""
    for old, new in _OLD_KEYS.items():
        if old in raw and new not in raw:
            raw[new] = raw[old]
    return raw


# ── public API ──

def load_config(path: str) -> PipelineConfig:
    raw: dict[str, str] = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            key, _, val = line.partition("\t")
            if not val:
                key, _, val = line.partition(" ")
            raw[key.strip()] = val.strip()

    raw = _normalise(raw)

    # Detect which format
    has_new = "LIB_R1" in raw
    has_old = "FORWARD_LIB_R1" in raw

    if has_new:
        lib_r1 = _read_path(raw, "LIB_R1")
        lib_r2 = _read_path(raw, "LIB_R2")
    elif has_old:
        # Both orientations use the same lib in practice
        lib_r1 = _read_path(raw, "FORWARD_LIB_R1")
        lib_r2 = _read_path(raw, "FORWARD_LIB_R2")
    else:
        raise ConfigError("Missing LIB_R1 / FORWARD_LIB_R1")

    if _read_optional(raw, "FORWARD_TAG"):
        fwd_tag = _read_optional(raw, "FORWARD_TAG")
        rev_tag = _read_optional(raw, "REVERSE_TAG")
    else:
        # Legacy: FORWARD_LIB_TAG was already the rev-comp
        fwd_tag = _read_optional(raw, "FORWARD_LIB_TAG", "")
        rev_tag = _read_optional(raw, "REVERSE_LIB_TAG", "")

    if not fwd_tag or not rev_tag or not lib_r1 or not lib_r2:
        raise ConfigError(
            "Missing required keys: LIB_R1/LIB_R2 + FORWARD_TAG/REVERSE_TAG "
            "(or legacy FORWARD_LIB_R1/FORWARD_LIB_R2)"
        )

    # Prefix / outdir
    prefix = raw.get("PREFIX", raw.get("prefix", "SAMPLE"))
    outdir = _read_path(raw, "OUTDIR") or Path.cwd() / "outdir" / prefix

    # Files that must exist
    ctrl = _read_path(raw, "CTRL")
    adapter = _read_path(raw, "ADAPTER") or Path.cwd() / "bin/adapters.txt"
    bin_dir = _read_path(raw, "BIN")

    # Reference
    genome = _read_optional(raw, "GENOME", "hg38")
    index = _read_path(raw, "INDEX")
    ref = _read_path(raw, "REF")
    chromsize = _read_path(raw, "CHROMSIZE")

    if not index or not ref or not chromsize:
        raise ConfigError("Missing INDEX / REF / CHROMSIZE")

    # Optional
    grna_file = _read_path(raw, "GRNA")
    blacklist = _read_path(raw, "BLACKLIST")
    if not blacklist and bin_dir:
        blacklist = bin_dir.parent / "data" / f"{genome}.blacklist.bed"

    return PipelineConfig(
        prefix=prefix,
        outdir=outdir,
        lib_r1=lib_r1,
        lib_r2=lib_r2,
        forward_tag=fwd_tag,
        reverse_tag=rev_tag,
        index=index,
        ref=ref,
        chromsize=chromsize,
        genome=genome,
        minlen=_read_int(raw, "MINLEN", 50),
        readlen=_read_int(raw, "READLEN", 150),
        maxins=_read_int(raw, "MAXINS", 1000),
        threads=_read_int(raw, "THREAD", 4),
        min_support_readcount=_read_int(raw, "MinSupportReadCount", 1),
        min_cutting_event_count=_read_int(raw, "MinCuttingEventCount", 2),
        max_mismatch=_read_int(raw, "MaxMismatch", 6),
        max_gap=_read_int(raw, "MaxGap", 2),
        max_gap_mismatch=_read_int(raw, "MaxGapMismatch", 4),
        grna_file=grna_file,
        ctrl=ctrl,
        adapter=adapter,
        blacklist=blacklist if blacklist and blacklist.exists() else None,
        bin_dir=bin_dir,
    )


def warn_tag_orientation(r2_path: Path, tag: str, n_check: int = 1000) -> None:
    """Sample R2 to see if tag or rev-comp is more common (catch orientation bug)."""
    if not r2_path.exists():
        return
    rev = reverse_complement(tag)
    tag_cnt = rev_cnt = 0
    opener = gzip.open if str(r2_path).endswith(".gz") else open
    with opener(r2_path, "rt") as f:
        for _ in range(n_check):
            try:
                next(f); seq = next(f); next(f); next(f)
            except StopIteration:
                break
            if tag in seq:
                tag_cnt += 1
            elif rev in seq:
                rev_cnt += 1
    if rev_cnt > tag_cnt:
        warnings.warn(
            f"Tag orientation may be reversed. "
            f"Rev-comp found {rev_cnt}/{n_check} vs tag {tag_cnt}/{n_check}. "
            f"Check FORWARD_TAG / REVERSE_TAG in config."
        )
