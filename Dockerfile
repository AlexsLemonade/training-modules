FROM rocker/tidyverse:4.0.3
LABEL maintainer="ccdl@alexslemonade.org"
WORKDIR /rocker-build/

RUN apt-get update && apt-get install -y --no-install-recommends apt-utils
RUN apt-get install dialog apt-utils -y

RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    build-essential \
    libxml2 \
    libxml2-dev \
    libsqlite-dev \
    libmariadbd-dev \
    libmariadbclient-dev \
    libpq-dev \
    libssh2-1-dev \
    pandoc \
    libmagick++-dev \
    time \
    python3-pip \
    git \
    gcc \
    make \
    g++ \
    libboost-all-dev \
    liblzma-dev \
    libbz2-dev \
    ca-certificates \
    zlib1g-dev \
    curl \
    unzip \
    autoconf \
    libglpk-dev

# FastQC
RUN apt update && apt install -y fastqc

# fastp 
ENV FASTP_VERSION 0.20.1
RUN git clone https://github.com/OpenGene/fastp.git
RUN cd fastp && \
    git checkout tags/v${FASTP_VERSION} -b v${FASTP_VERSION} && \
    make && \
    sudo make install

# MultiQC
ENV MULTIQC_VERSION 1.9
RUN pip3 install multiqc==${MULTIQC_VERSION}

# Snakemake
ENV SNAKEMAKE_VERSION 5.19.3
RUN pip3 install snakemake==${SNAKEMAKE_VERSION}

# Salmon
ENV SALMON_VERSION 1.4.0
WORKDIR /usr/local/src
RUN wget --quiet https://github.com/COMBINE-lab/salmon/releases/download/v${SALMON_VERSION}/salmon-${SALMON_VERSION}_linux_x86_64.tar.gz
RUN tar xzf salmon-${SALMON_VERSION}_linux_x86_64.tar.gz && \
    rm -f salmon-${SALMON_VERSION}_linux_x86_64.tar.gz && \
    ln -s /usr/local/src/salmon-latest_linux_x86_64/bin/salmon /usr/local/bin/salmon

# renv
ENV RENV_VERSION 0.12.5-2
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

# CRAN packages
ENV CRAN https://packagemanager.rstudio.com/cran/__linux__/focal/2021-01-26
RUN install2.r --error --skipinstalled -r ${CRAN} \
    caTools \
    colorspace \
    fastqcr \
    foreach \
    ggpubr \
    glmnet \
    gplots \
    gtools \
    hexbin \
    iterators \
    igraph \
    optparse \
    pheatmap \
    remotes \
    rjson \
    Rtsne \
    umap \
    BiocManager
RUN rm -rf /tmp/downloaded_packages


ENV BIOC_VERSION 3.12
RUN R -e "options(repos = c(CRAN = '${CRAN}'), warn = 2); \
    BiocManager::install(c( \
        'alevinQC', \
        'AnnotationHub', \
        'ComplexHeatmap', \
        'ConsensusClusterPlus', \
        'DESeq2', \
        'ensembldb', \
        'GEOquery', \
        'org.Cf.eg.db',\
        'org.Dr.eg.db', \
        'org.Hs.eg.db', \
        'org.Mm.eg.db', \
        'qvalue', \
        'scater', \
        'scran',  \
        'tximport', \
        'vsn'), \
    version = '${BIOC_VERSION}', update = FALSE)"

# # Use renv for R packages
# ENV RENV_VERSION 0.12.5-2
# RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
# RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

# WORKDIR /renv
# COPY renv.lock renv.lock
# RUN R -e 'renv::consent(provided = TRUE)'
# RUN R -e 'renv::restore()'

WORKDIR /home/rstudio

