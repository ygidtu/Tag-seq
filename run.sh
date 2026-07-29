#!/usr/bin/env bash
# ==========================================
# Tag-seq 分析流程 — 完整运行脚本
# 使用方法: bash run.sh
# ==========================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CONFIG="$SCRIPT_DIR/configs/TAG-1.txt"
LOG_DIR="$SCRIPT_DIR/outdir/TAG-1"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# ── 激活环境 ──
if [ -d "$SCRIPT_DIR/.venv" ]; then
    source "$SCRIPT_DIR/.venv/bin/activate"
elif command -v tagseq &>/dev/null; then
    :  # tagseq 已在 PATH 中
else
    echo "[ERROR] 请先安装环境: cd $SCRIPT_DIR && uv venv && uv pip install -e ."
    exit 1
fi

mkdir -p "$LOG_DIR"

echo "========================================="
echo "  Tag-seq Pipeline Start: $(date)"
echo "  Config: $CONFIG"
echo "========================================="

# ── Step 1: 准备输出目录 ──
echo "[$(date '+%H:%M:%S')] Step 1/4: create-makefile"
tagseq -c "$CONFIG" create-makefile 2>&1 | tee "$LOG_DIR/create_makefile.log"

# ── Step 2: 比对与靶点检测 ──
echo "[$(date '+%H:%M:%S')] Step 2/4: align (ODN → cutadapt → UMI → STAR → dedup → targets)"
tagseq -c "$CONFIG" align 2>&1 | tee "$LOG_DIR/align.log"

# ── Step 3: 脱靶鉴定 ──
echo "[$(date '+%H:%M:%S')] Step 3/4: find-target"
tagseq -c "$CONFIG" find-target 2>&1 | tee "$LOG_DIR/find_target.log"

# ── Step 4: 统计报告 ──
echo "[$(date '+%H:%M:%S')] Step 4/4: report"
tagseq -c "$CONFIG" report 2>&1 | tee "$LOG_DIR/report.log"

echo "========================================="
echo "  Pipeline Complete: $(date)"
echo "  Results: $LOG_DIR"
echo "========================================="

# ── 输出关键结果汇总 ──
echo ""
echo "===== 结果摘要 ====="
grep -E "ODN|Trim|Align|Unique|Targets" "$LOG_DIR/report.log" 2>/dev/null || echo "  (查看 $LOG_DIR/report.log)"
