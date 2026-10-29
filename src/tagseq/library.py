"""文库结构解析与校验 (依据 TAG-seq信息.xlsx)。

原版文库 (legacy):   P5 + UMI + barcode + adapter_junction + genome + Tag + P7
改造后文库 (modified): Universal_primer + UMI + barcode + adapter_junction + genome + Tag
  Universal primer = TAGCACCACGGATGG (15)
  UMI              = 8  (模式 NNWNNWNN, W=A/T)
  barcode          = TTTAGCTC (8)
  adapter junction = ACTACTAATACGACT (15)
  Tag              = 27 (位于 R1 3' 端 / R2 5' 端)

本模块把"5' 结构"配置化, 供 extract_umi 提取 UMI/barcode 并剥离人工序列,
同时提供序列级校验 (primer 前缀 / UMI 简并模式 / barcode / junction)。
"""
from __future__ import annotations
import re
from dataclasses import dataclass, field
from typing import Optional

from .types import LibrarySpec

# IUPAC 简并码 -> 字符类
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]", "K": "[GT]", "M": "[AC]",
    "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}


def pattern_to_regex(pat: str) -> Optional[re.Pattern]:
    """将 IUPAC 简并模式编译为整串匹配的正则; 空模式返回 None。"""
    if not pat:
        return None
    try:
        return re.compile("".join(IUPAC.get(c.upper(), re.escape(c)) for c in pat) + r"\Z")
    except re.error:
        return None


@dataclass
class ReadParse:
    umi: str
    barcode: str
    insert: str                       # 处理后用于比对的序列
    primer_ok: bool
    barcode_ok: Optional[bool] = None  # None = 未设 barcode, 不校验
    umi_ok: Optional[bool] = None      # None = 未设 umi_pattern, 不校验
    junction_ok: Optional[bool] = None
    fp_off: int = 0                    # 5' 结构总长 (strip_5p 时 insert 起点)
    umi_start: int = 0                 # UMI 在 seq 中的起点
    umi_end: int = 0                   # UMI 在 seq 中的终点
    issues: list[str] = field(default_factory=list)


def validate_spec(spec: LibrarySpec) -> list[str]:
    """校验 LibrarySpec 自身字段完整性/自洽性, 返回错误列表 (空=通过)。"""
    errs: list[str] = []
    if spec.umi_len < 0:
        errs.append(f"umi_len 非法: {spec.umi_len}")
    if spec.barcode_len < 0:
        errs.append(f"barcode_len 非法: {spec.barcode_len}")
    if spec.umi_pattern and len(spec.umi_pattern) != spec.umi_len:
        errs.append(f"umi_pattern 长度 {len(spec.umi_pattern)} != umi_len {spec.umi_len}")
    if spec.umi_pattern and pattern_to_regex(spec.umi_pattern) is None:
        errs.append(f"umi_pattern 含非法字符: {spec.umi_pattern}")
    if spec.barcode and len(spec.barcode) != spec.barcode_len:
        errs.append(f"barcode 长度 {len(spec.barcode)} != barcode_len {spec.barcode_len}")
    for nm, s in (("primer", spec.primer), ("junction", spec.junction)):
        if s and not set(s.upper()) <= set("ACGTN"):
            errs.append(f"{nm} 含非碱基字符: {s}")
    return errs


def parse_5p(seq: str, spec: LibrarySpec) -> ReadParse:
    """按 spec 解析 read 的 5' 结构。

    - spec.primer 非空: 要求 seq 以 primer 开头, 否则 primer_ok=False (不剥离)。
    - spec.primer 为空 (legacy): 视为无 primer, 从 offset 0 解析 UMI。
    - strip_5p=True: insert 剥离 primer+UMI+barcode+junction (改造后文库)。
    - strip_5p=False: insert 仅去 UMI, 保留其余 (原版文库)。
    """
    issues: list[str] = []
    pl = len(spec.primer)
    if spec.primer:
        if not seq.startswith(spec.primer):
            issues.append("primer_mismatch")
            return ReadParse("", "", seq, False, issues=issues)
    us, ue = pl, pl + spec.umi_len
    umi = seq[us:ue]
    rx = pattern_to_regex(spec.umi_pattern)
    umi_ok = (rx.match(umi) is not None) if rx else None
    if rx and not umi_ok:
        issues.append(f"umi_pattern_violation:{umi}")

    off = ue
    bc = seq[off:off + spec.barcode_len] if spec.barcode_len else ""
    off += spec.barcode_len
    barcode_ok = None
    if spec.barcode:
        barcode_ok = (bc == spec.barcode)
        if not barcode_ok:
            issues.append(f"barcode_mismatch:{bc}")

    junction_ok = None
    if spec.junction:
        j = seq[off:off + len(spec.junction)]
        junction_ok = (j == spec.junction)
        if not junction_ok:
            issues.append(f"junction_mismatch:{j}")
        off += len(spec.junction)

    insert = seq[off:] if spec.strip_5p else (seq[:us] + seq[ue:])
    return ReadParse(umi, bc, insert, True, barcode_ok, umi_ok, junction_ok, off, us, ue, issues)


_COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


def classify_extract(r1: str, spec: LibrarySpec, fwd_tag: str, rev_tag: str) -> Optional[dict]:
    """混测 R1 文库判定与字段提取 (以 5' 主锚点判型, 不强制 3' 读穿)。

    L 文库 R1: UP + UMI + BC + J_L + genome (+ Tag-L, 仅短 insert 读穿时可见)
    R 文库 R1: Tag-R + genome (+ J_R + BC + UMI + UP_R, 仅短 insert 读穿时可见)
    返回 {type, umi, barcode, tag, g_start, g_end}; 两个 5' 主锚点都不/都命中 -> None(ambiguous)。
    genome = r1[g_start:g_end]; 3' 元件读不到时置空并降置信, 不丢弃整条 read。
    """
    UP, JL = spec.primer, spec.junction
    ul, bl = spec.umi_len, spec.barcode_len
    pl, jl = len(UP), len(JL)
    five = pl + ul + bl + jl
    upr = rc(UP)
    is_L = bool(UP) and r1.startswith(UP)
    is_R = bool(rev_tag) and r1.startswith(rev_tag)
    if is_L == is_R:
        return None
    if is_L:
        umi = r1[pl:pl + ul]
        bc = r1[pl + ul:pl + ul + bl]
        tpos = r1.rfind(fwd_tag)
        if tpos >= five:
            tag = r1[tpos:tpos + len(fwd_tag)]; g_end = tpos
        else:
            tag = ""; g_end = len(r1)
        return {"type": "L", "umi": umi, "barcode": bc, "tag": tag,
                "g_start": five, "g_end": g_end}
    g_start = len(rev_tag)
    if r1.endswith(upr):
        e = len(r1)
        umi = r1[e - ul - len(upr):e - len(upr)]
        bc = r1[e - ul - len(upr) - bl:e - ul - len(upr)]
        g_end = e - (jl + bl + ul + len(upr))
    else:
        umi, bc = "", ""
        g_end = len(r1)
    return {"type": "R", "umi": umi, "barcode": bc, "tag": r1[:g_start],
            "g_start": g_start, "g_end": g_end}


def parse_r2(r2: str, spec: LibrarySpec) -> dict:
    """R2 结构: Tag + genome + J_R + BC + UMI + UP_R (tag 开头, UP_R 结尾)。
    返回 {g_start, g_end, umi}; 若 UP_R 缺失则保留全长(g_start=0,g_end=len)。"""
    upr, jr = rc(spec.primer), rc(spec.junction)
    suffix = len(jr) + spec.barcode_len + spec.umi_len + len(upr)
    if upr and r2.endswith(upr) and len(r2) >= spec.tag_len + suffix:
        end = len(r2)
        return {"g_start": spec.tag_len, "g_end": end - suffix,
                "umi": r2[end - spec.umi_len - len(upr):end - len(upr)]}
    return {"g_start": 0, "g_end": len(r2), "umi": ""}
