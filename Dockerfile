FROM mambaorg/micromamba:2-debian12

# 避免交互式安装提示
ENV DEBIAN_FRONTEND=noninteractive \
    TZ=Asia/Shanghai \
    LC_ALL=C \
    # 优化 pip 和 conda 行为
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    MAMBA_ROOT_PREFIX=/opt/conda

# 设置工作目录
WORKDIR /tmp

# 切换到 root 用户进行系统配置
USER root

# ============================================
# 阶段 1: 系统依赖安装（合并 RUN 减少层数）
# ============================================
RUN set -eux; \
    # 修改 debian.sources 文件中的源为阿里云镜像
    if [ -f /etc/apt/sources.list.d/debian.sources ]; then \
        sed -i 's|deb.debian.org|mirrors.aliyun.com|g' /etc/apt/sources.list.d/debian.sources; \
        sed -i 's|security.debian.org|mirrors.aliyun.com|g' /etc/apt/sources.list.d/debian.sources; \
    fi; \
    # 基础工具安装（合并安装，清理缓存）
    apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        wget \
        git \
        build-essential \
        python3 \
        python3-pip \
        python3-dev \
        bzip2 \
        cmake \
        default-jdk \
        ncurses-dev \
        tzdata \
        unzip \
        zlib1g \
        zlib1g-dev \
        libnss-sss \
        fastqc \
        emboss \
        adapterremoval \
        bedtools \
        samtools \
        bedops \
        librsvg2-bin \
        r-base \
        r-base-dev \
        && \
    # 清理 apt 缓存（减小镜像体积）
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/cache/apt/* /tmp/* /var/tmp/*

# ============================================
# 阶段 2: STAR 安装（使用多阶段构建思想，但单阶段优化）
# ============================================
ENV STAR_VERSION=2.7.10b
RUN set -eux; \
    wget -q --timeout=30 --tries=3 \
        https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz && \
    tar -xzf ${STAR_VERSION}.tar.gz && \
    cp STAR-${STAR_VERSION}/bin/Linux_x86_64_static/STAR /usr/local/bin/ && \
    chmod +x /usr/local/bin/STAR && \
    rm -rf STAR-${STAR_VERSION} ${STAR_VERSION}.tar.gz

# ============================================
# 阶段 3: R 包安装（使用阿里云镜像加速）
# ============================================
RUN set -eux; \
    Rscript -e " \
        options(repos = c(CRAN = 'https://mirrors.aliyun.com/CRAN/')); \
        install.packages(c('ggplot2', 'reshape2', 'RColorBrewer'), \
        dependencies = TRUE, \
        Ncpus = $(nproc) \
    )" && \
    # 清理 R 缓存
    rm -rf /tmp/Rtmp*

# ============================================
# 阶段 4: Conda/Mamba 环境创建（优化层数和缓存）
# ============================================
# 创建 ucsc-bedgraphtobigwig 环境
RUN micromamba create -y -p /opt/ucsc-bedgraphtobigwig \
        -c bioconda -c conda-forge \
        ucsc-bedgraphtobigwig && \
    micromamba clean -afy

# 创建 umitools 环境
RUN micromamba create -y -p /opt/umitools \
        -c bioconda -c conda-forge \
        umi_tools && \
    micromamba clean -afy

# ============================================
# 阶段 6: 安装 Python 包到 base 环境
# ============================================
RUN set -eux; \
    micromamba create -y -p /opt/venv -c conda-forge \
        svgwrite \
        scipy \
        numpy \
        biopython && \
    micromamba clean -afy


# ============================================
# 阶段 7: 应用代码复制和软链接创建
# ============================================
COPY ./ /opt/Tag-seq

# ============================================
# 阶段 5: 修复 umi_tools 中的 pysam.sort 问题
# ============================================
RUN mv /opt/Tag-seq/data/dedup.py /opt/umitools/lib/python3.12/site-packages/umi_tools/dedup.py



# 创建软链接（修复路径问题）
RUN ln -sf /opt/ucsc-bedgraphtobigwig/bin/bedGraphToBigWig /opt/Tag-seq/bin/bedGraphToBigWig && \
    # 设置可执行权限
    chmod +x /opt/Tag-seq/bin/*.pl /opt/Tag-seq/bin/*.py 2>/dev/null || true && \
    # 创建必要的目录
    mkdir -p /opt/Tag-seq/temp /opt/Tag-seq/results

# ============================================
# 阶段 8: 环境变量配置
# ============================================
ENV PATH="/opt/Tag-seq/bin:/opt/venv/bin:${PATH}"

# 设置工作目录
WORKDIR /workspace

# 切换回默认用户以提高安全性
USER mambauser


# 健康检查（可选）
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD which samtools && which STAR && which bedGraphToBigWig || exit 1

# 默认命令
ENTRYPOINT []
CMD []