#!/bin/bash
# 三样本 mm=3 全量: K5/K8 全流程, C1 只重跑 find-target+plot (align 产物复用 rmODN 口径)
cd /rawdata1/projects/TagSeq_ZXY/Tag-seq
mkdir -p outdir
LOG=/rawdata1/projects/TagSeq_ZXY/Tag-seq/outdir/run_LMY_mm3_$(date +%Y%m%d).log
for s in K5 K8; do
  echo "===== [$s] all (mm=3) start $(date '+%F %T') =====" | tee -a "$LOG"
  .venv/bin/tagseq -c configs/$s.txt all >> "$LOG" 2>&1
  echo "===== [$s] end $(date '+%F %T') rc=$? =====" | tee -a "$LOG"
done
echo "===== [C1] find-target+plot (mm=3) start $(date '+%F %T') =====" | tee -a "$LOG"
.venv/bin/tagseq -c configs/C1.txt find-target >> "$LOG" 2>&1
.venv/bin/tagseq -c configs/C1.txt plot >> "$LOG" 2>&1
echo "===== [C1] end $(date '+%F %T') rc=$? =====" | tee -a "$LOG"
echo "ALLDONE $(date '+%F %T')" | tee -a "$LOG"
