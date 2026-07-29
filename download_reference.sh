#!/usr/bin/env bash
# ==========================================
# 下载参考基因组 (hg38) 并构建 STAR 索引
# 无需 root 权限
# ==========================================
set -euo pipefail

REF_DIR="${1:-./reference/hg38}"
THREADS="${2:-8}"

echo "[INFO] Downloading hg38 reference to: $REF_DIR"
mkdir -p "$REF_DIR"

# 下载 FASTA
FASTA="$REF_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
if [ ! -f "$FASTA" ]; then
    echo "[INFO] Downloading hg38 FASTA (this may take a while)..."
    wget -q -O "$FASTA.gz" \
        https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    gunzip "$FASTA.gz"
    echo "[INFO] FASTA downloaded: $FASTA"
else
    echo "[INFO] FASTA already exists, skip download."
fi

# 构建 FASTA 索引
if [ ! -f "$FASTA.fai" ]; then
    echo "[INFO] Indexing FASTA ..."
    samtools faidx "$FASTA"
fi

# 生成 chrom.sizes
CHROM="$REF_DIR/hg38.chrom.sizes"
if [ ! -f "$CHROM" ]; then
    awk -v OFS='\t' '{print $1, $2}' "$FASTA.fai" > "$CHROM"
    echo "[INFO] Created: $CHROM"
fi

# 构建 STAR 索引
STAR_INDEX="$REF_DIR/star_index"
if [ ! -d "$STAR_INDEX" ]; then
    echo "[INFO] Building STAR index (this may take a while)..."
    mkdir -p "$STAR_INDEX"
    STAR --runMode genomeGenerate \
         --genomeDir "$STAR_INDEX" \
         --genomeFastaFiles "$FASTA" \
         --runThreadN "$THREADS" \
         --genomeSAindexNbases 14
    echo "[INFO] STAR index built: $STAR_INDEX"
else
    echo "[INFO] STAR index already exists, skip."
fi

echo ""
echo "========================================="
echo "  Reference genome ready!"
echo "========================================="
echo ""
echo "  Update your config.txt:"
echo "    INDEX       $STAR_INDEX"
echo "    GENOME      hg38"
echo "    REF         $FASTA"
echo "    CHROMSIZE   $CHROM"
echo ""
