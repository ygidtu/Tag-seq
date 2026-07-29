# Tag-seq 脱靶分析流程

## 目录

- [1. 概述](#1-概述)
- [2. 系统要求](#2-系统要求)
- [3. 依赖安装](#3-依赖安装)
- [4. 数据准备](#4-数据准备)
- [5. 配置文件说明](#5-配置文件说明)
- [6. 运行流程](#6-运行流程)
- [7. 结果解读](#7-结果解读)
- [8. 问题排查](#8-问题排查)
- [9. 当前环境已知问题](#9-当前环境已知问题)

---

## 1. 概述

Tag-seq 是一种基于双端标签（Tag）和 UMI（Unique Molecular Identifier）的脱靶检测方法。其核心策略是：

- **双文库设计**：使用正向（forward）和反向（reverse）两种 Tag 引物构建两个独立的文库
- **UMI 去重**：通过 8bp UMI 精确去除 PCR 重复
- **峰值检测**：在参考基因组上检测切割位点富集区域
- **交叉验证**：要求正负链 reads 同时支持，或同链但两个方向引物文库均支持

### 整体分析流程图

```
FASTQ (R1/R2)
    │
    ├─ ① FastQC 原始质控
    ├─ ② remove_ODN — 去除 Tag 序列
    ├─ ③ AdapterRemoval — 修剪接头 / 低质量碱基
    ├─ ④ umi_tools extract — 提取 8bp UMI
    ├─ ⑤ STAR 比对 → 参考基因组 (hg19/mm10)
    ├─ ⑥ umi_tools dedup — UMI 去重
    ├─ ⑦ 提取单端 read → 候选切割位点
    ├─ ⑧ bedtools merge (10bp 窗口) → 聚类
    ├─ ⑨ detect.peaks.py — scipy 峰值检测
    ├─ ⑩ 4 个方向文件合并 → find_targetsite.v2.pl
    │       (plus.plus / plus.minus / minus.plus / minus.minus)
    ├─ ⑪ Smith-Waterman 比对 gRNA 序列
    └─ ⑫ 输出 off-target 列表 + 可视化
```

---

## 2. 系统要求

- **操作系统**：Linux (CentOS / Ubuntu 均测试通过)
- **内存**：≥ 32 GB RAM（处理 30M 片段约需 30 分钟）
- **存储**：参考基因组索引约需 30 GB，中间文件约为原始数据的 5-10 倍
- **软件**：Perl ≥ 5, Python ≥ 3.8, R ≥ 4.0

---

## 3. 依赖安装

### 3.1 核心软件

| 软件 | 版本要求 | 用途 | 安装方式 |
|------|---------|------|---------|
| STAR | ≥ 2.7 | 比对至参考基因组 | `apt install star` 或源码编译 |
| fastqc | ≥ 0.12 | FASTQ 质控 | `apt install fastqc` |
| AdapterRemoval | ≥ 2.3 | 接头修剪 | `apt install adapterremoval` 或源码编译 |
| samtools | ≥ 1.10 | SAM/BAM 处理 | `apt install samtools` |
| bedtools | ≥ 2.30 | BED 文件操作 | `apt install bedtools` |
| bedops | ≥ 2.4 | BED 文件操作（bedmap） | `apt install bedops` |
| umi_tools | ≥ 1.1 | UMI 提取与去重 | `conda install -c bioconda umi_tools` |
| bedGraphToBigWig | — | bedGraph → bigWig 转换 | `conda install -c bioconda ucsc-bedgraphtobigwig` |
| EMBOSS water | ≥ 6.6 | Smith-Waterman 比对 | `apt install emboss` |
| picard | ≥ 2.25 | BAM 处理 | `apt install picard-tools` |

### 3.2 Python 依赖

```bash
pip install numpy scipy pandas biopython loguru tqdm
```

### 3.3 R 依赖（仅可视化）

```R
install.packages("RIdeogram")
```

### 3.4 安装参考（本项目环境）

本项目环境中的安装方式记录在 `./README.md`：

```bash
# 1. 安装 bedGraphToBigWig
mamba create -p ~/.local/ucsc-bedgraphtobigwig -c bioconda -c conda-forge ucsc-bedgraphtobigwig
ln -sf ~/.local/bedgraphtobigwig/bin/bedGraphToBigWig Tag-seq/bin/bedGraphToBigWig

# 2. 安装 umi_tools
mamba create -p ~/.local/umitools -c bioconda -c conda-forge umi_tools

# 3. umi_tools 因在新版 pysam 下 sort 有问题，需 patch：
#    编辑 /path/to/umi_tools/dedup.py 第 372 行附近
#    将 pysam.sort 替换为直接调用 samtools sort
```

### 3.5 准备好参考基因组索引

```bash
# hg19
mkdir -p /path/to/ref/hg19 && cd /path/to/ref/hg19
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.chrom.sizes
gunzip hg19.fa.gz
samtools faidx hg19.fa
STAR --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles ./hg19.fa --runThreadN 32

# mm10
mkdir -p /path/to/ref/mm10 && cd /path/to/ref/mm10
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
# ... 重复相同步骤
```

---

## 4. 数据准备

### 4.1 输入文件结构

Tag-seq 需要两组双端 FASTQ 文件——对应正向（forward）和反向（reverse）Tag 引物文库：

```
data/
├── forward_R1.fq.gz    # 正向文库 R1
├── forward_R2.fq.gz    # 正向文库 R2
├── reverse_R1.fq.gz    # 反向文库 R1
└── reverse_R2.fq.gz    # 反向文库 R2
```

**注意**：同一样本的正反向文库可以使用相同的 FASTQ 文件（即两次 PCR 合并测序），但 Tag 序列不同。

### 4.2 sgRNA 列表文件

格式：TSV，每行一个 sgRNA，包含 ID、序列和 PAM 基序：

```
# sgRNA_ID  sequence          PAM
hBCL11A-Cas-g1  CTTCCTGGAGCCTGTGATAAAAGC  NTTM
Pcsk9-Cas9-g11  GGCGATGGTCTTGATGGCAGGG  NGG
```

### 4.3 对照文件 (control.bed)

用于过滤非特异性整合位点（如 AAVS1 位点），BED 格式，9 列。如果不需要对照，可在配置中设置为 `none`。

### 4.4 数据预处理（可选）

如果需要从原始测序数据中拆分 barcode 或合并文件，可使用本项目自带的辅助脚本：

```bash
# 按 barcode 拆分 FASTQ
python ./split_fq.py -i raw_R1.fq.gz -j raw_R2.fq.gz -b barcodes.txt -o output_dir/

# 随机抽取子样本测试
python ./sample_fq.py -1 input_R1.fq.gz -2 input_R2.fq.gz -o subsample -f 0.1
```

---

## 5. 配置文件说明

### 5.1 配置模板

| 参数 | 说明 | 示例值 |
|------|------|--------|
| `PREFIX` | 样本前缀（用于命名中间文件） | `HEK` / `B16` |
| `OUTDIR` | 输出目录 | `./outdir/HEK_A2` |
| `FORWARD_LIB_R1` | 正向文库 R1 FASTQ | `/path/to/A2_R1.fq` |
| `FORWARD_LIB_R2` | 正向文库 R2 FASTQ | `/path/to/A2_R2.fq` |
| `FORWARD_LIB_TAG` | 正向 Tag 序列（29 bp） | `TGCGATAACACGCATTTCGCATAAG` |
| `REVERSE_LIB_R1` | 反向文库 R1 FASTQ | `/path/to/A2_R1.fq` |
| `REVERSE_LIB_R2` | 反向文库 R2 FASTQ | `/path/to/A2_R2.fq` |
| `REVERSE_LIB_TAG` | 反向 Tag 序列（29 bp） | `ATCTCTGAGCCTTATGCGAAATGCG` |
| `CTRL` | 对照 BED 文件路径 | `./control.bed` 或 `none` |
| `GRNA` | sgRNA 列表文件路径 | `./sgrna.lst` |
| `MinSupportReadCount` | 候选位点最低 read 支持数 | `1` |
| `MinCuttingEventCount` | 最低切割事件计数 | `2` |
| `MaxMismatch` | off-target 检测最大错配数 | `6` |
| `MaxGap` | off-target 检测最大 gap | `2` |
| `MaxGapMismatch` | gap 区域内最大错配 | `4` |
| `MINLEN` | 修剪后最短 read 长度 | `50` |
| `READLEN` | 测序读长 | `150` |
| `MAXINS` | 最大插入片段长度 | `1000` |
| `THREAD` | 线程数 | `16` |
| `ADAPTER` | 接头序列文件 | `Tag-seq/bin/adapters.txt` |
| `INDEX` | STAR 索引目录 | `/path/to/hg19/star_index` |
| `GENOME` | 基因组名称 | `hg19` / `mm10` |
| `REF` | 参考基因组 FASTA | `/path/to/hg19.fa` |
| `CHROMSIZE` | 染色体大小文件 | `/path/to/hg19.chrom.sizes` |

上述参数中，**软件路径**部分需根据实际安装位置调整：

| 参数 | 默认值 |
|------|--------|
| `BIN` | Tag-seq 的 bin 目录 |
| `FASTQC` | `/usr/bin/fastqc` |
| `STAR` | `/usr/local/bin/STAR` |
| `BEDTOOLS` | `/usr/bin/bedtools` |
| `SAMTOOLS` | `/usr/bin/samtools` |
| `PICARD` | `/usr/bin/picard` |
| `AdapterRemoval` | `/usr/local/bin/AdapterRemoval` |
| `umi_tools` | `umi_tools`（需在 PATH 中） |
| `water` | `/usr/bin/water` |
| `bedops` | `/usr/bin/bedops` |

### 5.2 本项目已有配置文件

所有配置文件位于本仓库的 `configs/` 目录下：

| 配置文件 | 样本 | 基因组 | sgRNA |
|---------|------|--------|-------|
| `HEK_config.txt` | HEK293T-A2 | hg19 | `sgrnas/sgrna_list.txt` |
| `HEK_config_v2.txt` | HEK293T-A3 | hg19 | `sgrnas/sgrna_list_v2.txt` |
| `HEK_config_v3.txt` | HEK293T-sg1617 | hg19 | `sgrnas/sgrna_list_v3.txt` |
| `HEK_config_v4.txt` | HEK293T-A6 | hg19 | `sgrnas/sgrna_list_v4.txt` |
| `B16_config.txt` | B16-原始 | mm10 | `sgrnas/sgrna_list_mouse.txt` |
| `B16_config_A3.txt` | B16-A3 | mm10 | `sgrnas/sgrna_list_mouse.txt` |
| `B16_config_A11.txt` | B16-A11 | mm10 | `sgrnas/sgrna_list_mouse_A11.txt` |
| `B16_config_A12.txt` | B16-A12 | mm10 | `sgrnas/sgrna_list_mouse_A12.txt` |

---

## 6. 运行流程

### 6.1 快速运行（一键执行）

```bash
# 进入工作目录
cd /path/to/Cas12f/res/offtarget

# 确保配置文件中的路径正确
# 然后执行（建议使用 nohup 后台运行）
nohup perl Tag-seq/bin/run_guideseq.pl configs/HEK_config.txt all \
    > logs/HEK_A2.log 2> logs/HEK_A2.err &

nohup perl Tag-seq/bin/run_guideseq.pl configs/B16_config_A12.txt all \
    > logs/B16_A12.log 2> logs/B16_A12.err &
```

### 6.2 分步运行

如果需要在某一步中断后继续，可分段执行：

```bash
# Step 1: 创建配置文件和 Makefile
perl Tag-seq/bin/run_guideseq.pl configs/HEK_config.txt create_makefile

# Step 2: 比对与去重（this step runs make for _plus and _minus）
perl Tag-seq/bin/run_guideseq.pl configs/HEK_config.txt align

# Step 3: 检测 off-target（对每个 sgRNA 逐一运行 find_targetsite.v2.pl）
perl Tag-seq/bin/run_guideseq.pl configs/HEK_config.txt find_target
```

### 6.3 运行时间

- **30M 片段**：约 30 分钟（使用 16 线程）
- 主要耗时在 STAR 比对和 umi_tools dedup 阶段
- off-target 检测与 sgRNA 数量成正比

---

## 7. 结果解读

### 7.1 输出目录结构

运行完成后，在 `OUTDIR` 下生成以下目录：

```
outdir/PREFIX/
├── config.txt                         # 生成的配置文件副本
├── sample.lst                         # 样本列表
├── PREFIX.log / PREFIX.err            # 主日志
├── stat.txt                           # 统计汇总
├── PREFIX_plus/                       # 正向文库分析结果
│   ├── 00datafilter/                  # 质控结果
│   │   ├── PREFIX.R1_fastqc.html
│   │   ├── PREFIX.rmODN.R1.fq         # 去除 Tag 后
│   │   ├── PREFIX.Trim.R1.umis.fq     # 提取 UMI 后
│   │   └── ...
│   ├── 01alignment/                   # 比对结果
│   │   ├── PREFIX.Aligned.sortedByCoord.out.bam
│   │   └── PREFIX.Aligned.sortedByCoord.out.dedup.bam  # UMI 去重后
│   ├── 02potentialTargets/            # 候选切割位点
│   │   ├── PREFIX.plus.proximal.sorted    # 正链候选位点
│   │   └── PREFIX.minus.proximal.sorted   # 负链候选位点
│   └── 03visual/                      # 可视化
│       ├── PREFIX.plus.bw             # 正链 bigWig
│       └── PREFIX.minus.bw            # 负链 bigWig
├── PREFIX_minus/                      # 反向文库（结构同上）
│   └── ...
└── sgRNA_ID.find.target/              # 每个 sgRNA 的 off-target 结果
    ├── grna.fa                        # sgRNA 序列
    ├── sgRNA_ID.parsing_water_for_visualization.offtarget.bed  # 最终 off-target
    ├── sgRNA_ID.all.sites.merged.confirmed    # 所有确认位点
    └── sites/
        └── monitoring_offtargets.log  # 位点监控日志
```

### 7.2 核心结果文件

**off-target BED 文件**：`{sgRNA}.parsing_water_for_visualization.offtarget.bed`

格式（12 列 BED）：

```
chr1    10111   10112   sgRNA1_E_minus_minus_2_9,sgRNA1_E_plus_minus_1_13  0   29  0   12
```

- 列 1-3：染色体、起始、结束
- 列 4：位点 ID（包含来源文库、方向、read 计数信息）
- 列 5：score（0）
- 列 6：正链 plus 文库中正链 read 数
- 列 7：正链 plus 文库中负链 read 数
- 列 8：负链 minus 文库中正链 read 数
- 列 9：负链 minus 文库中负链 read 数

### 7.3 统计文件：`stat.txt`

包含每个样本的 QC 和比对统计，包括：
- 原始 reads 数
- 去除 ODN 后 reads 数
- 修剪后 reads 数
- 比对率
- UMI 去重后 unique reads 数

### 7.4 筛选标准

候选 off-target 的确认条件（任一满足即可标记为潜在 DSB）：

1. 同时在 **正链和负链** 检测到 reads 的位点
2. 同一条链但 **正向和反向两个 Tag 引物文库** 均检测到 reads 的位点
3. 使用 Smith-Waterman 将 gRNA 序列与位点侧翼序列比对，计算错配数

---

## 8. 问题排查

### 8.1 常见错误

| 错误信息 | 原因 | 解决方案 |
|---------|------|---------|
| `Can't locate Some/Module.pm` | 缺少 Perl 模块 | `cpan install Some::Module` |
| `STAR: command not found` | STAR 未安装或不在 PATH | 安装 STAR 或修正配置路径 |
| `umi_tools: command not found` | umi_tools 未安装 | `pip install umi_tools` 或修正路径 |
| `bedGraphToBigWig: not found` | UCSC 工具未安装 | 参照 3.4 安装并创建软链接 |
| `No module named 'tqdm'` | Python 包缺失 | `pip install tqdm` |
| `samtools sort: failed` | umi_tools 内部调用问题 | 修改 umi_tools dedup.py 改用系统 samtools |

### 8.2 日志检查

```bash
# 查看主要进度
tail -f logs/HEK_A2.log

# 查看错误
tail -f logs/HEK_A2.err

# 查看 make 执行日志
cat outdir/HEK_A2/HEK_plus/log
cat outdir/HEK_A2/HEK_plus/err

# 查看各步骤详细日志
cat outdir/HEK_A2/HEK_plus/00datafilter/*.log
cat outdir/HEK_A2/HEK_plus/01alignment/*.log
```

### 8.3 中间文件清理

运行完成后，可选择删除中间文件以释放空间：

```bash
# 删除中间 FASTQ（保留为数据分析输出所需即可自行决定）
rm -rf outdir/PREFIX/PREFIX_plus/00datafilter/*.fq
rm -rf outdir/PREFIX/PREFIX_minus/00datafilter/*.fq

# 可以删除的是输入文件的副本，例如，仅需保留 Aligned.sortedByCoord.out.dedup.bam
```

---

## 9. 当前环境已知问题

### 9.1 配置文件路径修复（必须修改后使用）

本项目已有的 HEK 配置文件存在 **过期路径**：

| 配置文件 | 问题参数 | 当前值 | 应修正为 |
|---------|---------|--------|---------|
| `HEK_config.txt` | `BIN` | `/home/user/Tag-seq/bin` | `Tag-seq/bin` 的绝对路径 |
| `HEK_config_v2.txt` | `BIN` | `/home/user/Tag-seq/bin` | 同上 |
| `HEK_config_v3.txt` | `BIN` | `/home/user/Tag-seq/bin` | 同上 |
| `HEK_config_v4.txt` | `BIN` | `/home/user/Tag-seq/bin` | 同上 |
| `HEK_config.txt` | `umi_tools` | `umi_tools`（bare） | 需写入完整路径或确保在 PATH 中 |
| `HEK_config_v2.txt` | `umi_tools` | `umi_tools`（bare） | 同上 |

**B16 配置文件路径正常**，可直接使用（路径如下）。但请注意这些路径是服务器绝对路径，**其他用户需全部替换为自己的路径**。

*以上路径问题不会影响项目的长期运行时效果，但在复制到另一台服务器时，需要重新适配全部路径。*

### 9.2 软件依赖检查

当前环境中已安装的依赖状态：

| 依赖 | 状态 | 路径/版本 |
|------|------|-----------|
| Perl | 已安装 | v5.40.1 |
| Python | 已安装 | v3.13.14 |
| R | 已安装 | v4.5.1 |
| STAR | 已安装 | `/usr/local/bin/STAR` |
| bedtools | 已安装 | `/usr/bin/bedtools` |
| samtools | 已安装 | `/usr/bin/samtools`, `/opt/samtools/1.23.1/samtools` |
| fastqc | 已安装 | `/usr/bin/fastqc` |
| bedGraphToBigWig | 已安装 | `Tag-seq/bin/bedGraphToBigWig`（symlink） |
| umi_tools | 已安装 | `/home/zym/.local/umitools/bin/umi_tools` |
| AdapterRemoval | 已安装 | `/usr/local/bin/AdapterRemoval` |
| water (EMBOSS) | 已安装 | `/usr/bin/water` |
| bedops | 已安装 | `/home/zym/.local/bin/bedops` |
| picard | 已安装 | `/usr/bin/picard` |
| numpy/scipy/pandas | 已安装 | 正常 |
| biopython/loguru | 已安装 | 正常 |
| tqdm | **未安装** | 运行 `analysis.py` 前需 `pip install tqdm` |

### 9.3 代码兼容性修复

以下是本项目对原始 Tag-seq 所做的主要修改：

1. **detect.peaks.py**：从 Python 2 改为 Python 3（`print` 改为函数，`numpy` 类型处理适配）
2. **umi_tools patch**：在 dedup.py 中将 `pysam.sort` 替换为直接调用系统 `samtools sort`
3. **bedGraphToBigWig**：通过 conda 安装并软链接至 `Tag-seq/bin/` 下

### 9.4 在新服务器上的部署步骤

1. 安装所有依赖软件（参考第 3 节）
2. 下载并构建参考基因组索引（hg19 和/或 mm10）
3. 复制 Tag-seq 目录至服务器
4. 创建配置文件，**所有路径改为服务器实际路径**
5. 安装 bedGraphToBigWig 并创建 `Tag-seq/bin/bedGraphToBigWig` 软链接
6. 安装 Python 依赖：`pip install numpy scipy pandas biopython loguru tqdm`
7. 运行测试：`perl Tag-seq/bin/run_guideseq.pl configs/B16_config_A12.txt all`
