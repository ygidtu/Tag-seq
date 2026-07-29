#!/bin/bash
# ==========================================
# Tag-seq 环境检查脚本
# 检测所有依赖工具和参考基因组文件
# ==========================================

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BIN_DIR="$SCRIPT_DIR/bin"
DATA_DIR="$SCRIPT_DIR/data"

echo "========================================="
echo "  Tag-seq Environment Check"
echo "========================================="

# ========== 工具检查 ==========
echo ""
echo "[Tools] Checking required software..."

check_tool() {
    local name=$1
    local path=$2
    if command -v "$path" &>/dev/null; then
        echo "  [OK] $name: $path"
        return 0
    elif [ -x "$path" ] || [ -f "$path" ]; then
        echo "  [OK] $name: $path"
        return 0
    else
        echo "  [FAIL] $name: $path NOT FOUND"
        return 1
    fi
}

check_tool "perl" "$(which perl)"
check_tool "python3" "$(which python3)"
check_tool "STAR" "/usr/local/bin/STAR"
check_tool "AdapterRemoval" "/usr/local/bin/AdapterRemoval"
check_tool "fastqc" "/usr/bin/fastqc"
check_tool "samtools" "/usr/bin/samtools"
check_tool "bedtools" "/usr/bin/bedtools"
check_tool "picard" "/usr/bin/picard"
check_tool "water (emboss)" "/usr/bin/water"
check_tool "bedops" "/home/zym/.local/bin/bedops"
check_tool "umi_tools" "/home/zym/.local/umitools/bin/umi_tools"
check_tool "rsvg-convert" "$(which rsvg-convert 2>/dev/null || echo 'not found')"

# Check python packages
echo ""
echo "[Python packages] Checking..."
python3 -c "import svgwrite" 2>/dev/null && echo "  [OK] svgwrite" || echo "  [FAIL] svgwrite (install: pip install svgwrite)"
python3 -c "import argparse" 2>/dev/null && echo "  [OK] argparse" || echo "  [FAIL] argparse"

# ========== 参考基因组检查 ==========
echo ""
echo "[Reference] Checking genome files..."

HG38_FA="/rawdata/reference/hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
HG38_INDEX="/rawdata/reference/hg38/star/2.7.10b"
HG38_CHROM="/rawdata/reference/hg38/hg38.chrom.sizes"

[ -f "$HG38_FA" ] && echo "  [OK] hg38 FASTA: $HG38_FA" || echo "  [FAIL] hg38 FASTA: $HG38_FA"
[ -d "$HG38_INDEX" ] && echo "  [OK] hg38 STAR index: $HG38_INDEX" || echo "  [FAIL] hg38 STAR index: $HG38_INDEX"
[ -f "$HG38_CHROM" ] && echo "  [OK] hg38 chrom.sizes: $HG38_CHROM" || echo "  [FAIL] hg38 chrom.sizes: $HG38_CHROM"

MM10_FA="/rawdata/reference/mm10/mm10.fa"
MM10_INDEX="/rawdata/reference/mm10/star_index"
MM10_CHROM="/rawdata/reference/mm10/mm10.chrom.sizes"

[ -f "$MM10_FA" ] && echo "  [OK] mm10 FASTA: $MM10_FA" || echo "  [WARN] mm10 FASTA: $MM10_FA (not needed for hg38 runs)"
[ -d "$MM10_INDEX" ] && echo "  [OK] mm10 STAR index: $MM10_INDEX" || echo "  [WARN] mm10 STAR index: $MM10_INDEX"
[ -f "$MM10_CHROM" ] && echo "  [OK] mm10 chrom.sizes: $MM10_CHROM" || echo "  [WARN] mm10 chrom.sizes: $MM10_CHROM"

# ========== 其他数据文件 ==========
echo ""
echo "[Data] Checking data files..."

ADAPTER="$BIN_DIR/adapters.txt"
[ -f "$ADAPTER" ] && echo "  [OK] Adapter list: $ADAPTER" || echo "  [FAIL] Adapter list: $ADAPTER"

BLACKLIST_HG38="$DATA_DIR/hg38.blacklist.bed"
BLACKLIST_MM10="$DATA_DIR/mm10.blacklist.bed"
[ -f "$BLACKLIST_HG38" ] && echo "  [OK] hg38 blacklist: $BLACKLIST_HG38" || echo "  [WARN] hg38 blacklist: $BLACKLIST_HG38"
[ -f "$BLACKLIST_MM10" ] && echo "  [OK] mm10 blacklist: $BLACKLIST_MM10" || echo "  [WARN] mm10 blacklist: $BLACKLIST_MM10"

CTRL_BED="$SCRIPT_DIR/test/control.bed"
[ -f "$CTRL_BED" ] && echo "  [OK] Ctrl BED: $CTRL_BED" || echo "  [WARN] Ctrl BED: $CTRL_BED (create blank if not needed)"

# ========== 脚本语法检查 ==========
echo ""
echo "[Scripts] Checking Perl syntax..."

for pl in "$BIN_DIR"/*.pl; do
    name=$(basename "$pl")
    perl -c "$pl" >/dev/null 2>&1 && echo "  [OK] $name" || echo "  [FAIL] $name (syntax error)"
done

echo ""
echo "========================================="
echo "  Check complete."
echo "========================================="
