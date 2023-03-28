FROM rocker/tidyverse:4.2.2
LABEL maintainer="ccdl@alexslemonade.org"
WORKDIR /rocker-build/

RUN apt-get update && apt-get install -y --no-install-recommends apt-utils
RUN apt-get install dialog apt-utils -y

RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    autoconf \
    build-essential \
    ca-certificates \
    curl \
    g++ \
    gcc \
    git \
    groff \
    less \
    libboost-all-dev \
    libbz2-dev \
    libfftw3-dev \
    libglpk-dev \
    liblzma-dev \
    libmagick++-dev \
    libmariadb-dev-compat \
    libmariadbd-dev \
    libpq-dev \
    libproj-dev \
    libsqlite-dev \
    libssh2-1-dev \
    libxml2 \
    libxml2-dev \
    make \
    pandoc \
    python3-pip \
    time \
    unzip \
    zlib1g-dev

# AWS
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install

# FastQC
RUN apt-get update && apt-get install -y fastqc

# fastp
ENV FASTP_VERSION 0.20.1
RUN git clone https://github.com/OpenGene/fastp.git
RUN cd fastp && \
    git checkout tags/v${FASTP_VERSION} -b v${FASTP_VERSION} && \
    make && \
    make install

# MultiQC
ENV MULTIQC_VERSION 1.9
RUN pip3 install multiqc==${MULTIQC_VERSION}

# Snakemake
ENV SNAKEMAKE_VERSION 7.25.0
RUN pip3 install snakemake==${SNAKEMAKE_VERSION}

# Salmon
ENV SALMON_VERSION 1.10.0
WORKDIR /usr/local/src
RUN wget --quiet https://github.com/COMBINE-lab/salmon/releases/download/v${SALMON_VERSION}/salmon-${SALMON_VERSION}_linux_x86_64.tar.gz
RUN tar xzf salmon-${SALMON_VERSION}_linux_x86_64.tar.gz && \
    rm -f salmon-${SALMON_VERSION}_linux_x86_64.tar.gz && \
    ln -s /usr/local/src/salmon-latest_linux_x86_64/bin/salmon /usr/local/bin/salmon

# Use renv for R packages
ENV RENV_CONFIG_CACHE_ENABLED FALSE
RUN R -e "install.packages(c('renv', 'yaml', 'BiocManager'))"

WORKDIR /usr/local/renv
COPY renv.lock renv.lock
RUN R -e 'renv::consent(provided = TRUE); renv::restore()' && rm -rf ~/.local/share/renv


WORKDIR /home/rstudio
