"""检验: 各样本最强切割峰附近是否存在其 spacer (回答 'on-target 富集缺失' 是样本问题还是分析问题)"""
import glob
import pysam
from collections import defaultdict

REF = "/rawdata1/reference/hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
SPACER = {"C1": "TGTGGTAGCTGGGTGT", "K5": "AAGTGTCCATCAACAG", "K8": "AGGTGGAAAGAGAAAG"}
rc = lambda s: s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def top_peaks(sample, n=6):
    best = {}
    for f in glob.glob(f"outdir/{sample}/{sample}_*/02potentialTargets/*.proximal"):
        for line in open(f):
            p = line.split("\t")
            if len(p) < 6:
                continue
            key = (p[0], int(p[1]), int(p[2]))
            c = int(p[4])
            if c > best.get(key, 0):
                best[key] = c
    return sorted(best.items(), key=lambda kv: -kv[1])[:n]


def best_mm(seq, pat):
    out = (99, None, None)
    for i in range(len(seq) - len(pat) + 1):
        w = seq[i:i + len(pat)]
        d = sum(1 for a, b in zip(w, pat) if a != b)
        if d < out[0]:
            out = (d, i, w)
        if d == 0:
            break
    return out


fa = pysam.FastaFile(REF)
for s in ("C1", "K5", "K8"):
    sp = SPACER[s]
    print(f"=== {s}  spacer={sp}  rc={rc(sp)}")
    for (chrom, st, en), cnt in top_peaks(s):
        a, b = max(0, st - 300), en + 300
        seq = fa.fetch(chrom, a, b).upper()
        if len(seq) < len(sp):
            print(f"   chr{chrom}:{st}-{en} reads={cnt:4d}  区域过短")
            continue
        d1, _, w1 = best_mm(seq, sp)
        d2, _, w2 = best_mm(seq, rc(sp))
        print(f"   chr{chrom}:{st}-{en} reads={cnt:4d}  spacer最佳mm={d1} ({w1})  rc最佳mm={d2} ({w2})")
