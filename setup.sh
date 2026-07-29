#!/usr/bin/env bash
# ==========================================
# Tag-seq 一键环境部署脚本
# 无需 root 权限，使用 Conda 安装所有依赖
# ==========================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ENV_NAME="tagseq"

echo "========================================="
echo "  Tag-seq Environment Setup"
echo "========================================="
echo ""

# ---- 检测 Conda ----
if command -v conda &>/dev/null; then
    CONDA_CMD="conda"
elif command -v mamba &>/dev/null; then
    CONDA_CMD="mamba"
elif [ -f "$HOME/miniconda3/bin/conda" ]; then
    echo "[INFO] Found conda at ~/miniconda3/bin/conda"
    CONDA_CMD="$HOME/miniconda3/bin/conda"
elif [ -f "$HOME/anaconda3/bin/conda" ]; then
    echo "[INFO] Found conda at ~/anaconda3/bin/conda"
    CONDA_CMD="$HOME/anaconda3/bin/conda"
else
    echo "[INFO] Conda not found. Installing Miniconda ..."
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh
    bash /tmp/miniconda.sh -b -p "$HOME/miniconda3"
    rm /tmp/miniconda.sh
    echo "[INFO] Miniconda installed."
    CONDA_CMD="$HOME/miniconda3/bin/conda"
    # 初始化 conda
    "$CONDA_CMD" init bash &>/dev/null || true
fi

export PATH="$HOME/miniconda3/bin:$PATH"

# ---- 安装 mamba（加速依赖解析） ----
if ! command -v mamba &>/dev/null; then
    echo "[INFO] Installing mamba (faster dependency solver) ..."
    "$CONDA_CMD" install -y -c conda-forge mamba
    CONDA_CMD="mamba"
fi

# ---- 创建 tagseq 环境 ----
echo "[INFO] Creating conda environment '$ENV_NAME' ..."
if "$CONDA_CMD" env list | grep -q "^$ENV_NAME "; then
    echo "[INFO] Environment '$ENV_NAME' already exists, updating ..."
    "$CONDA_CMD" env update -f "$SCRIPT_DIR/environment.yml" --prune
else
    "$CONDA_CMD" env create -f "$SCRIPT_DIR/environment.yml"
fi

echo "[INFO] Environment '$ENV_NAME' created successfully."

# ---- 创建激活脚本 ----
cat > "$SCRIPT_DIR/activate_env.sh" << 'ENVEOF'
#!/usr/bin/env bash
# 激活 Tag-seq 环境
# 用法: source activate_env.sh
ENV_DIR="$(cd "$(dirname "$0")" && pwd)"
export PATH="$HOME/miniconda3/envs/tagseq/bin:$ENV_DIR/bin:$PATH"
echo "[INFO] Tag-seq environment activated."
echo "       Run: perl bin/run_guideseq.pl config.txt all"
ENVEOF
chmod +x "$SCRIPT_DIR/activate_env.sh"

# ---- 验证安装 ----
echo ""
echo "========================================="
echo "  Verifying installation ..."
echo "========================================="
source "$SCRIPT_DIR/activate_env.sh" 2>/dev/null || true

# 使用 conda run 验证
verify_tool() {
    local name=$1
    local cmd=$2
    if "$CONDA_CMD" run -n "$ENV_NAME" which "$cmd" &>/dev/null; then
        echo "  [OK] $name"
    else
        echo "  [FAIL] $name (not found)"
    fi
}

verify_tool "perl"        "perl"
verify_tool "python3"     "python3"
verify_tool "STAR"        "STAR"
verify_tool "AdapterRemoval" "AdapterRemoval"
verify_tool "fastqc"      "fastqc"
verify_tool "samtools"    "samtools"
verify_tool "bedtools"    "bedtools"
verify_tool "picard"      "picard"
verify_tool "umi_tools"   "umi_tools"
verify_tool "water"       "water"
verify_tool "bedops"      "bedops"
verify_tool "bedGraphToBigWig" "bedGraphToBigWig"
verify_tool "rsvg-convert" "rsvg-convert"

"$CONDA_CMD" run -n "$ENV_NAME" python3 -c "import svgwrite" &>/dev/null && \
    echo "  [OK] svgwrite (Python)" || echo "  [FAIL] svgwrite (Python)"

# ---- 检查参考基因组 ----
echo ""
echo "[Reference genome]"
REF_DIR="/rawdata/reference"
if [ -d "$REF_DIR" ]; then
    echo "  [OK] Reference directory: $REF_DIR"
else
    echo "  [WARN] Reference directory not found: $REF_DIR"
    echo "  [HINT] Download reference genome and update INDEX/REF/CHROMSIZE in config"
fi

echo ""
echo "========================================="
echo "  Setup complete!"
echo "========================================="
echo ""
echo "  Usage:"
echo "    1. Activate environment:"
echo "       source activate_env.sh"
echo ""
echo "    2. Run analysis:"
echo "       perl bin/run_guideseq.pl config.txt all"
echo ""
