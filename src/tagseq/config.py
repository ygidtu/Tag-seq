from __future__ import annotations
import gzip
from loguru import logger
from pathlib import Path
from typing import Optional
from .exceptions import ConfigError
from .types import PipelineConfig, LibrarySpec
from .library import validate_spec

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
    if not val or val.lower() == "none": return None
    return _resolve(val)

def _read_int(raw: dict, key: str, default: int = 0) -> int:
    try: return int(_read_optional(raw, key, str(default)))
    except: return default

def _read_bool(raw: dict, key: str, default: bool = False) -> bool:
    v = _read_optional(raw, key, "").strip().lower()
    if v in ("1", "true", "yes", "y", "on"): return True
    if v in ("0", "false", "no", "n", "off"): return False
    return default

def load_config(path: str) -> PipelineConfig:
    raw = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"): continue
            key, _, val = line.partition("\t")
            if not val: key, _, val = line.partition(" ")
            raw[key.strip()] = val.strip()
    has_new = "LIB_R1" in raw
    has_old = "FORWARD_LIB_R1" in raw
    if has_new:
        lib_r1 = _read_path(raw, "LIB_R1")
        lib_r2 = _read_path(raw, "LIB_R2")
    elif has_old:
        lib_r1 = _read_path(raw, "FORWARD_LIB_R1")
        lib_r2 = _read_path(raw, "FORWARD_LIB_R2")
    else:
        raise ConfigError("Missing LIB_R1 / FORWARD_LIB_R1")
    fwd_tag = _read_optional(raw, "FORWARD_TAG") or _read_optional(raw, "FORWARD_LIB_TAG", "")
    rev_tag = _read_optional(raw, "REVERSE_TAG") or _read_optional(raw, "REVERSE_LIB_TAG", "")
    if not fwd_tag or not rev_tag or not lib_r1 or not lib_r2:
        raise ConfigError("Missing LIB_R1/LIB_R2 + FORWARD_TAG/REVERSE_TAG")
    prefix = raw.get("PREFIX", raw.get("prefix", "SAMPLE"))
    outdir = _read_path(raw, "OUTDIR") or Path.cwd() / "outdir" / prefix
    ctrl = _read_path(raw, "CTRL")
    adapter = _read_path(raw, "ADAPTER") or Path.cwd() / "bin/adapters.txt"
    bin_dir = _read_path(raw, "BIN")
    genome = _read_optional(raw, "GENOME", "hg38")
    index = _read_path(raw, "INDEX")
    ref = _read_path(raw, "REF")
    chromsize = _read_path(raw, "CHROMSIZE")
    if not index or not ref or not chromsize:
        raise ConfigError("Missing INDEX / REF / CHROMSIZE")
    grna_file = _read_path(raw, "GRNA")
    blacklist = _read_path(raw, "BLACKLIST")
    if not blacklist and bin_dir:
        blacklist = bin_dir.parent / "data" / f"{genome}.blacklist.bed"

    # 文库类型与 5' 结构 (见 library.py / TAG-seq信息.xlsx)
    lib_type = _read_optional(raw, "LIB_TYPE", "legacy").strip().lower()
    if lib_type not in ("legacy", "modified"):
        raise ConfigError(f"Unknown LIB_TYPE: {lib_type} (expect 'legacy' or 'modified')")
    # 比对入口的过滤口径; 默认 rmODN 与 TAG-1/3/5 既有结果一致
    odn_mode = _read_optional(raw, "ODN_MODE", "rmODN").strip().lower()
    if odn_mode not in ("rmodn", "mixed"):
        raise ConfigError(f"Unknown ODN_MODE: {odn_mode} (expect 'rmODN' or 'mixed')")
    odn_mode = "rmODN" if odn_mode == "rmodn" else "mixed"
    if lib_type == "modified":
        d = dict(primer="TAGCACCACGGATGG", umi_len=8, umi_pattern="NNWNNWNN",
                 barcode="TTTAGCTC", barcode_len=8, junction="ACTACTAATACGACT",
                 tag_len=25, strip_5p=True)
    else:
        d = dict(primer="", umi_len=8, umi_pattern="", barcode="", barcode_len=8,
                 junction="", tag_len=25, strip_5p=False)
    primer = _read_optional(raw, "UNIVERSAL_PRIMER") or _read_optional(raw, "UMI_PREFIX") or d["primer"]
    spec = LibrarySpec(
        lib_type=lib_type,
        primer=primer,
        umi_len=_read_int(raw, "UMI_LEN", d["umi_len"]),
        umi_pattern=_read_optional(raw, "UMI_PATTERN", d["umi_pattern"]),
        barcode=_read_optional(raw, "BARCODE", d["barcode"]),
        barcode_len=_read_int(raw, "BARCODE_LEN", d["barcode_len"]),
        junction=_read_optional(raw, "JUNCTION", d["junction"]),
        tag_len=_read_int(raw, "TAG_LEN", d["tag_len"]),
        strip_5p=_read_bool(raw, "STRIP_5P", d["strip_5p"]),
        strict=_read_bool(raw, "STRICT", False),
    )
    spec_errs = validate_spec(spec)
    if spec_errs:
        raise ConfigError("LibrarySpec 配置错误: " + "; ".join(spec_errs))

    return PipelineConfig(
        prefix=prefix, outdir=outdir, lib_r1=lib_r1, lib_r2=lib_r2,
        forward_tag=fwd_tag, reverse_tag=rev_tag,
        index=index, ref=ref, chromsize=chromsize, genome=genome,
        minlen=_read_int(raw, "MINLEN", 50),
        readlen=_read_int(raw, "READLEN", 150),
        maxins=_read_int(raw, "MAXINS", 1000),
        threads=_read_int(raw, "THREAD", 4),
        umi_len=spec.umi_len,
        umi_offset=len(spec.primer),
        umi_prefix=spec.primer,
        lib_type=lib_type,
        odn_mode=odn_mode,
        min_support_readcount=_read_int(raw, "MinSupportReadCount", 1),
        min_cutting_event_count=_read_int(raw, "MinCuttingEventCount", 2),
        max_mismatch=_read_int(raw, "MaxMismatch", 6),
        max_gap=_read_int(raw, "MaxGap", 2),
        max_gap_mismatch=_read_int(raw, "MaxGapMismatch", 4),
        grna_file=grna_file, ctrl=ctrl, adapter=adapter, bin_dir=bin_dir,
        blacklist=blacklist if blacklist and blacklist.exists() else None,
        lib_spec=spec,
    )

def validate_config(cfg: PipelineConfig) -> list[str]:
    """Return list of warnings/errors about config. Empty = all good."""
    issues = []
    pairs = [
        ("LIB_R1", cfg.lib_r1),
        ("LIB_R2", cfg.lib_r2),
        ("ADAPTER", cfg.adapter),
        ("INDEX (Genome)", cfg.index / "Genome"),
        ("REF", cfg.ref),
        ("CHROMSIZE", cfg.chromsize),
    ]
    for name, path in pairs:
        if not path.exists():
            issues.append(f"{name} not found: {path}")
    if cfg.grna_file and not cfg.grna_file.exists():
        issues.append(f"GRNA file not found: {cfg.grna_file}")
    if cfg.ctrl and cfg.ctrl != Path("none") and not cfg.ctrl.exists():
        issues.append(f"CTRL file not found: {cfg.ctrl}")
    return issues

def warn_tag_orientation(r2_path: Path, tag: str, n_check: int = 1000) -> None:
    if not r2_path.exists(): return
    rev = reverse_complement(tag)
    tag_cnt = rev_cnt = 0
    opener = gzip.open if str(r2_path).endswith(".gz") else open
    with opener(r2_path, "rt") as f:
        for _ in range(n_check):
            try:
                next(f); seq = next(f); next(f); next(f)
            except StopIteration: break
            if tag in seq: tag_cnt += 1
            elif rev in seq: rev_cnt += 1
    if rev_cnt > tag_cnt:
        logger.warning("Tag orientation reversed: rev-comp {}/{} vs tag {}/{}", rev_cnt, n_check, tag_cnt, n_check)
