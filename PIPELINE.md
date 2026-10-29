# Tag-seq 分析流程详解(算法级拆解)

> 配套 `README.md`。README 讲"怎么装、怎么跑";本文讲"每一步在做什么、为什么、如何影响脱靶效率"。
> 代码位置以 `src/tagseq/` 为准;参数以 `configs/TAG-1/3/5.txt` 当前值为准。

## 目录
1. 术语速查
2. 总览与模块映射
3. 阶段 A `create-makefile`
4. 阶段 B `align`
5. 阶段 C `find-target`
6. 阶段 D `report`
7. 阶段 E `plot`
8. 脱靶效率计算
9. 序列组分拆解
10. 关键参数与"旋钮"
11. 原版(Perl) vs 改造版差异
12. 已知问题 / 待确认

---

## 1. 术语速查
| 术语 | 含义 |
|---|---|
| reads | 测序得到的短 DNA 字符串(输入) |
| gRNA | 引导 Cas 切割的向导序列(脱靶比对的探针) |
| on-target | 与 gRNA 完全匹配(mm=0)的预期切割位点 |
| off-target | 与 gRNA 相似(1–`MaxMismatch` 错配)的误切位点 |
| UMI | 分子标签,用于把 PCR 重复归并为一条原始分子 |
| peak/cluster/proximal | 比对末端按距离聚合成的"候选切割峰" |
| confirmed | 通过方向/支持阈值的可信切割位点 |

---

## 2. 总览与模块映射
```
configs/TAG-*.txt ──config.py──► PipelineConfig(+LibrarySpec)
        │
cli.py ─► pipeline.py(编排)
   create-makefile → align → find-target → report → plot
   align:        parse_mixed → _run_cutadapt → run_star → dedup_umi → detect_targets
   find-target:  _read_grnas → _merge → _confirm_sites → _slop/ext → SW(_score_alignment) → 过滤 → 脱靶表
   统计/绘图:    report.py / qc.py / plot.py
```

---

## 3. 阶段 A `create-makefile`(`pipeline.step_create_makefile`)
- **输入**:config。
- **输出**:目录骨架 `outdir/<S>/{<S>_plus,<S>_minus}/{00datafilter,01alignment,02potentialTargets}`。
- **判定逻辑**:`warn_tag_orientation` 抽查 R2 中 tag 与其反向互补出现次数,反向更多则告警。
- **意义**:准备可复现输出结构;早期发现"tag 方向填反"这类会让全流程失效的输入错误。

---

## 4. 阶段 B `align`(`pipeline.step_align`)

### B1 文库解析 / ODN 过滤
当前代码用 `parse_mixed`(`pipeline.py` + `library.classify_extract`/`parse_r2`);生产 TAG 结果目前沿用旧 `remove_odn`(`trim.py`)产物。
- **输入**:`LIB_R1/R2`、`LibrarySpec`、`forward_tag`、`reverse_tag`。
- **输出**:`<S>.mixed.R1/R2.fq.gz`(仅基因组片段,read 名带 `:UMI_x:LIB_L|R`);旧法为 `<S>.rmODN.R1/R2`。
- **判定逻辑**:
  - `parse_mixed`:5′ 锚点判型——`startswith(UP)`→L,`startswith(Tag-R)`→R;都不/都命中→`ambiguous` 丢弃;读穿到 3′ tag/UP_R 则剥离。
  - `remove_odn`:R1 或 R2 任意位置含 4-tag(fwd/rev+rc)即保留,从 tag 后截断。
- **意义**:去除人工序列、提取 UMI、过滤背景;决定进入比对的有效 reads 量,直接影响峰与脱靶数。

### B2 `_run_cutadapt`
- **输入**:`mixed/rmODN.R1/R2`、`data/adapters.txt`。
- **输出**:`<S>.Trim.R1/R2.fq.gz` + `trim_report.txt`。
- **判定逻辑**:cutadapt 修 3′ 接头;`--quality-cutoff 5`、`--minimum-length MINLEN(50)`、`--trim-n`。
- **意义**:清除残留接头/低质量尾,避免接头被误判为切割信号。

### B3 `run_star`(`alignment.py`)
- **输入**:`Trim.R1/R2`、`INDEX`。
- **输出**:`...Aligned.sortedByCoord.out.bam` + `Log.final.out`。
- **判定逻辑**:STAR 双端比对,`--alignIntronMax 50`、`--outFilterScoreMinOverLread 0.5`。
- **意义**:把 reads 定位到参考坐标(染色体/位置/链),是切割位点的空间来源。

### B4 `dedup_umi`(`umi.py`)
- **输入**:比对 BAM。
- **输出**:`...dedup.bam` + `.flagstat`。
- **判定逻辑**:key=`(chrom, reference_start, UMI, is_read1)`,同 key 只留一条。
- **意义**:PCR 重复归一,防止某位点因扩增偏倚被虚高计数。

### B5 `detect_targets`(`targets.py`)
- **输入**:`dedup.bam`、`chromsize`。
- **输出**:`<S>.{plus,minus}.proximal`(每行=峰+重叠计数)、`plus.bed/minus.bed`、`clustered.region.bed`。
- **判定逻辑**:仅 `is_read2` 且 `MAPQ≥30`;取末端中点 `mid`;按 (chrom,mid) 排序后 `merge dist=10` 聚峰;`_overlap_count` 计每峰重叠 reads(第 5 列)。
- **意义**:把离散 read 末端聚成候选切割峰,峰支持 reads=信号强度,是脱靶分析的原料。

---

## 5. 阶段 C `find-target`(`pipeline.step_find_target` → `offtarget.find_offtargets`)

### C1 `_read_grnas`
- **输入**:`GRNA` 文件(行:`id  序列  PAM`)。**输出**:`[(id,seq,pam)]`。**判定**:空白切分≥3 字段。**意义**:提供脱靶比对探针。

### C2 读峰 + `_merge`
- **输入**:4 个 proximal。**输出**:`combined`、`merged`。**判定**:按 (chrom,start) 排序,距离≤10 合并、拼接 id、累加 count。**意义**:统一四方向证据、去重复峰。

### C3 `_confirm_sites`(位点确认)
- **输入**:`merged`、`combined`、`min_support(1)`、`cut_events(2)`。**输出**:`confirmed`(峰 + 四方向 reads)。
- **判定逻辑**:对每峰累加 4 方向 reads;达标方向数=count≥min_support 的个数;保留当且仅当达标方向数≥cut_events。
- **意义**:用"双向/多方向证据"剔除单侧噪声假峰,只让可信位点进入脱靶比对——直接决定脱靶候选数量。

### C4 `_slop` + `ext.fa`
- **输入**:`confirmed`、`REF`、`chromsize`。**输出**:`ext.bed`、`ext.fa`。
- **判定逻辑**:本版取峰 `[start-150, end+150]`(非中点),`pysam.FastaFile.fetch` 提序列。
- **意义**:为每峰准备含 gRNA 匹配区的局部序列;用峰范围可避免宽峰漏掉 protospacer(修复 on-target 丢失)。

### C5 Smith-Waterman + `_score_alignment`
- **输入**:gRNA、`ext.fa`(正/反链)。**输出**:`results`(命中 + 错配数 mm + 显示串 disp)。
- **判定逻辑**:`PairwiseAligner` 半全局(gRNA 全长对齐、target 端自由);`_score_alignment` 计 gRNA 区 mm(PAM 只展示不计分);保留 `mm≤max_mismatch(4)`。
- **意义**:mm=0→on-target,1–4→脱靶候选;`MaxMismatch` 是脱靶率核心旋钮。

### C6 blacklist/control 过滤(on-target 豁免)
- **输入**:`results`、`BLACKLIST`、`CTRL`(chrom 归一化)。**判定**:perfect 保留;非 perfect 命中黑名单/对照则删。**意义**:去已知假开放区/对照噪声,且不误杀真靶点。

### C7 导出 `offtarget.bed`
- **输出**:列 `chrom start end id reads disp gRNA strand`;首行=合并的 on-target。**判定**:perfect 合并、按 reads 降序、按位置去重。**意义**:最终脱靶清单,供脱靶效率与可视化。

---

## 6. 阶段 D `report`(`report.py`)
- **输入**:`rmODN.stat`、`trim_report`、`Log.final.out`、`flagstat`、`proximal`。
- **输出**:漏斗 Raw→ODN→Trim→Aligned→Unique→Targets。
- **判定逻辑**:`_count_min(col=5,minv=5)` 统计支持≥5 的峰数。
- **意义**:监控各步留存率,定位"reads 多但可靠少"发生在哪一步。

## 7. 阶段 E `plot`(`plot.py`/`qc.py`)
- **输出**:`statistics.pdf`、`mismatch_profile.pdf`、`qc_report.html`、`<S>_offtargets.pdf/.svg`。
- **意义**:可视化脱靶比对与错配分布,便于人工复核。

---

## 8. 脱靶效率计算(汇总脚本,非 pipeline 内置)
设 on-target=perfect 行(reads `R_on`),off-target=其余(`N_off` 个,reads `R_off`):
```
site_rate% = N_off / (1 + N_off) × 100        # 位点口径(对背景敏感, 常偏高)
read_rate% = R_off / (R_on + R_off) × 100     # reads 口径(反映真实切割负荷)
```
当前示例:read_rate TAG-1/3/5 ≈ 33%/48%/57%。缺阴性对照,绝对值仅供相对比较。

---

## 9. 序列组分拆解(改造后文库,5′→3′;第二行字符=组件来源)
L 文库 R1:
```
5'- TAGCACCACGGATGG  NNWNNWNN  TTTAGCTC  ACTACTAATACGACT  ACGTACGT...ACGT  ATCTCTGAGCCTTATGCGAAATGCG -3'
    PPPPPPPPPPPPPPP  MMMMMMMM  BBBBBBBB  JJJJJJJJJJJJJJJ  GGGGGGGGGGGGGGGG  TTTTTTTTTTTTTTTTTTTTTTTTT
```
R 文库 R1:
```
5'- CTTATGCGAAATGCGTGTTATCGCA  ACGTACGT...ACGT  AGTCGTATTAGTAGT  TTTAGCTC  NNWNNWNN  CCATCCGTGGTGCTA -3'
    TTTTTTTTTTTTTTTTTTTTTTTTT  GGGGGGGGGGGGGGGG  JJJJJJJJJJJJJJJ  BBBBBBBB  MMMMMMMM  PPPPPPPPPPPPPPP
```
`P`=Universal primer(15) `M`=UMI(8,`W`=A/T) `B`=barcode(8,`TTTAGCTC`) `J`=junction(15) `G`=genome(不定长) `T`=tag(25)。改造后文库去 P5/P7,改用 15bp Universal primer。

---

## 10. 关键参数与"旋钮"(`configs/TAG-*.txt`)
| 参数 | 当前值 | 作用步骤 | 对脱靶率影响 |
|---|---|---|---|
| MinSupportReadCount | 1 | C3 | 越小保留越多峰 |
| MinCuttingEventCount | 2 | C3 | 越大越严(需多方向证据) |
| MaxMismatch | 4 | C5 | 越大脱靶越多(核心) |
| ext 窗口 | 150(峰范围) | C4 | 影响 on-target 是否被覆盖 |
| BLACKLIST/CTRL | hg38/27行 | C6 | 去背景;on-target 豁免 |
| LIB_TYPE | modified | B1 | 决定文库解析 |

---

## 11. 原版(Perl) vs 改造版差异
| 维度 | 原版 | 改造版 |
|---|---|---|
| 运行 | `perl run_guideseq.pl`(docker `tag-seq:latest`) | `tagseq -c <cfg> <step>` |
| ODN | R2-only 单 tag 精确 | `remove_odn` 4-tag 任意 / `parse_mixed` 结构解析 |
| 脱靶比对 | EMBOSS `water` 局部 + 长度≥50% + 禁 gap | `PairwiseAligner` 半全局 + 允许 gap |
| 黑名单 | `bedtools intersect` 强制 | config `BLACKLIST` + chrom 归一化 + on-target 豁免 |
| 去重 | `umi_tools dedup` | 自实现 `dedup_umi` |
| 文库结构 | 写死脚本 | `LibrarySpec` 可配置 |

---

## 12. 已知问题 / 待确认
- `MaxGap/MaxGapMismatch` 在 `find_offtargets` 未生效(gap 不受控)。
- 无阴性对照 → 脱靶率绝对值不可单独解读;正式口径(site vs read)待定。
- 真实文库 read 方向存疑:~76% R1 以 primer 开头却仅 ~1% 含 tag;是否启用 `parse_mixed` 重跑 align 待确认。
- barcode 是否多样本复用(现固定 `TTTAGCTC`)。
- 代码改动未提交(git HEAD `8d0600d` + 本地 `library/config/umi/pipeline/offtarget/targets` 改动)。
