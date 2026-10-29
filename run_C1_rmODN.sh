#!/bin/bash
# 仅跑 C1, rmODN 口径 (与 TAG-1/3/5 一致); 记录各阶段耗时
cd /rawdata1/projects/TagSeq_ZXY/Tag-seq
mkdir -p outdir
LOG=/rawdata1/projects/TagSeq_ZXY/Tag-seq/outdir/run_C1_rmODN_$(date +%Y%m%d).log
T0=$(date +%s)
echo "===== C1 (ODN_MODE=rmODN) start $(date '+%F %T') =====" | tee -a "$LOG"
.venv/bin/tagseq -c configs/C1.txt all >> "$LOG" 2>&1
rc=$?
T1=$(date +%s)
echo "===== C1 end $(date '+%F %T') rc=$rc elapsed=$(( (T1-T0)/60 ))min =====" | tee -a "$LOG"
