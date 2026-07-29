# Tag-seq 脱靶分析流程

## 目录

- [1. 概述](#1-概述)
- [2. 安装](#2-安装)
- [3. 配置](#3-配置)
- [4. 运行](#4-运行)
- [5. 输出说明](#5-输出说明)
- [6. 问题排查](#6-问题排查)

---

## 1. 概述

Tag-seq 是一种基于双端标签（Tag）和 UMI（Unique Molecular Identifier）的脱靶检测方法。

### 分析流程

```
FASTQ (R1/R2)
    │
    ├─ ① ODN 去除 — 纯 Python (tagseq.trim)
    ├─ ② cutadapt — 接头修剪 / 质控
    ├─ ③ UMI 提取 — 纯 Python (tagseq.umi)
    ├─ ④ STAR 比对 → 参考基因组
    ├─ ⑤ UMI 去重 — pysam (tagseq.umi)
    ├─ ⑥ 候选切割位点检测 — 纯 Python (tagseq.targets)
    ├─ ⑦ Smith-Waterman 比对 — Bio.Align (tagseq.offtarget)
    └─ ⑧ 输出 off-target 列表
```

### 主要改进

- **纯 Python 实现**：35 个 Perl 脚本 → 12 个 Python 模块
- **依赖精简**：只需 STAR + cutadapt 两个外部工具
- **配置简化**：6 行必填（LIB_R1、LIB_R2、FORWARD_TAG、REVERSE_TAG、INDEX、REF、CHROMSIZE）
- **自动 .fq.gz 支持**：无需手动解压
- **uv 项目管理**：可复现安装，锁定依赖版本

---

## 2. 安装

### 2.1 系统要求

- **操作系统**：Linux
- **内存**：≥ 32 GB（处理 30M 片段）
- **存储**：参考基因组索引约 30 GB
- **Python**：≥ 3.10

### 2.2 一键安装（Conda）

```bash
# 自动安装 Miniconda + 所有依赖
bash setup.sh
source activate_env.sh
```

### 2.3 手动安装（uv）

```bash
# 安装 uv（如未安装）
curl -LsSf https://astral.sh/uv/install.sh | sh

# 创建虚拟环境并安装
cd Tag-seq
uv venv
uv pip install -e .

# 激活环境
source .venv/bin/activate
```

### 2.4 外部工具

| 工具 | 用途 | 安装 |
|------|------|------|
| STAR | 比对至参考基因组 | `conda install -c bioconda star` 或 `apt install star` |
| cutadapt | 接头修剪 | `pip install cutadapt` 或 `conda install -c bioconda cutadapt` |

### 2.5 参考基因组

```bash
# hg38 示例：从 Ensembl 下载并构建 STAR 索引
bash download_reference.sh ./reference/hg38 16
```

---

## 3. 配置

### 3.1 最小配置（6 行）

```
LIB_R1      ./raw_data/sample_R1.fq.gz
LIB_R2      ./raw_data/sample_R2.fq.gz
FORWARD_TAG CTTATGCGAAATGCGTGTTATCGCA
REVERSE_TAG CGCATTTCGCATAAGGCTCAGAGAT
INDEX       /path/to/star_index
REF         /path/to/genome.fa
CHROMSIZE   /path/to/genome.chrom.sizes
```

### 3.2 完整配置

参见 `configs/example_full.txt`。

### 3.3 Tag 序列说明

`FORWARD_TAG` / `REVERSE_TAG` 需填入**原始 tag 的反向互补序列**。脚本内部会再做一次反向互补，最终得到原始 tag 用于匹配 R2。

```
原始 tag:      TGCGATAACACGCATTTCGCATAAG
配置中填入:    CTTATGCGAAATGCGTGTTATCGCA（即其反向互补）
脚本处理后:    TGCGATAACACGCATTTCGCATAAG（恢复原始 tag，匹配 R2）
```

### 3.4 本项目已有配置文件

| 文件 | 说明 |
|------|------|
| `configs/example_minimal.txt` | 最小配置模板 |
| `configs/example_full.txt` | 全参数配置模板 |
| `configs/TAG-1.txt` | TAG-1 样本 (hg38) |
| `configs/test_py.txt` | 测试配置 (hg19) |

---

## 4. 运行

### 4.1 全流程

```bash
tagseq -c config.txt all
```

### 4.2 分步运行

```bash
# 创建输出目录
tagseq -c config.txt create-makefile

# 比对与靶点检测
tagseq -c config.txt align

# 脱靶鉴定
tagseq -c config.txt find-target

# 查看统计
tagseq -c config.txt report
```

### 4.3 续跑（跳过已完成步骤）

```bash
tagseq -c config.txt --resume all
```

---

## 5. 输出说明

### 5.1 目录结构

```
outdir/PREFIX/
├── PREFIX_plus/                     # 正向文库
│   ├── 00datafilter/
│   │   ├── PREFIX.rmODN.R1/R2.fq    # ODN 去除后
│   │   ├── PREFIX.Trim.R1/R2.fq     # 修剪后
│   │   └── PREFIX.Trim.R1/R2.umis.fq # UMI 提取后
│   ├── 01alignment/
│   │   ├── PREFIX.Aligned.sortedByCoord.out.bam
│   │   └── PREFIX.Aligned.sortedByCoord.out.dedup.bam
│   └── 02potentialTargets/
│       ├── PREFIX.plus.proximal       # 正链候选位点
│       └── PREFIX.minus.proximal      # 负链候选位点
├── PREFIX_minus/                     # 反向文库（结构同上）
└── gRNA_ID.find.target/
    └── gRNA_ID.parsing_water_for_visualization.offtarget.bed
```

### 5.2 统计报告

运行 `tagseq -c config.txt report` 输出：

```
  ODN pass:    7877 / 21098025 (0.04%)
  Trim:        36 / 7877 (0.46%)
  Aligned:     3 / 36 (8.33%)
  Unique:      3
  Dup rate:    0.0%
  Usable:      3
  Targets+:    0
  Targets-:    0
```

---

## 6. 问题排查

### 6.1 常见错误

| 错误 | 原因 | 解决 |
|------|------|------|
| `STAR: command not found` | STAR 未安装 | `conda install -c bioconda star` |
| `Tag orientation may be reversed` | FORWARD_TAG 方向错误 | 填入原始 tag 的反向互补 |
| `Read pairs do not match` | R1/R2 文件不配对 | 检查 FASTQ 文件是否配对 |
| `No module named 'pysam'` | Python 包缺失 | `uv pip install pysam` |

### 6.2 查看日志

```bash
# 详细输出
tagseq -v -c config.txt all

# 查看各步骤日志
cat outdir/PREFIX/PREFIX_plus/00datafilter/*.log
cat outdir/PREFIX/PREFIX_plus/01alignment/*.star.log
```
