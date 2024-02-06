# Build salmon from source in a separate image
FROM ubuntu:22.04 as build

ENV PACKAGES gcc g++ make cmake curl unzip ca-certificates \
    libboost-all-dev liblzma-dev libbz2-dev libcurl4-openssl-dev libdeflate-dev libisal-dev zlib1g-dev
RUN apt-get update -qq
RUN apt-get install -y --no-install-recommends ${PACKAGES}
WORKDIR /usr/local/src

# Build salmon
ENV SALMON_VERSION 1.10.1
RUN curl -LO https://github.com/COMBINE-lab/salmon/archive/refs/tags/v${SALMON_VERSION}.tar.gz
RUN tar xzf v${SALMON_VERSION}.tar.gz
RUN mkdir salmon-${SALMON_VERSION}/build
RUN cd salmon-${SALMON_VERSION}/build && \
    cmake -DCMAKE_INSTALL_PREFIX=/usr/local/salmon .. && \
    make && make install

# Build fastp
ENV FASTP_VERSION 0.23.4
RUN curl -LO https://github.com/OpenGene/fastp/archive/refs/tags/v${FASTP_VERSION}.tar.gz
RUN tar xzf v${FASTP_VERSION}.tar.gz
RUN cd fastp-${FASTP_VERSION} && \
    make && make install

# Main image with Biocconductor and other tools
FROM bioconductor/bioconductor_docker:RELEASE_3_18
LABEL maintainer="ccdl@alexslemonade.org"

WORKDIR /rocker-build/

RUN apt-get update -qq
RUN apt-get install -y --no-install-recommends \
    less

# copy salmon and fastp binaries from the build image
COPY --from=build /usr/local/salmon/ /usr/local/
COPY --from=build /usr/local/bin/fastp /usr/local/bin/fastp

# AWS
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install

# FastQC
RUN apt-get install -y --no-install-recommends fastqc

# MultiQC
ENV MULTIQC_VERSION 1.19
RUN pip install multiqc==${MULTIQC_VERSION}

# Snakemake
ENV SNAKEMAKE_VERSION 7.32.4
RUN pip install snakemake==${SNAKEMAKE_VERSION}

# Use renv for R packages
ENV RENV_CONFIG_CACHE_ENABLED FALSE
WORKDIR /usr/local/renv
COPY renv.lock renv.lock
RUN R -e "install.packages('renv')"
RUN MAKEFLAGS=-j$(nproc) R -e "renv::restore()" && rm -rf ~/.local/share/renv

WORKDIR /home/rstudio
