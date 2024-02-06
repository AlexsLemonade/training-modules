FROM bioconductor/bioconductor_docker:RELEASE_3_18
LABEL maintainer="ccdl@alexslemonade.org"

WORKDIR /rocker-build/

RUN apt-get update -qq

RUN apt-get install -y --no-install-recommends \
    less \
    libboost-all-dev
#     build-essential \
#     groff \
#     libboost-all-dev \
#     libglpk-dev \
#     libmariadb-dev-compat \
#     libmariadbd-dev \
#     libsqlite-dev \


# AWS
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install

# FastQC
RUN apt-get install -y --no-install-recommends fastqc

# fastp
ENV FASTP_VERSION 0.23.4
RUN curl -o /usr/local/bin/fastp http://opengene.org/fastp/fastp.${FASTP_VERSION} && \
    chmod +x /usr/local/bin/fastp

# MultiQC
ENV MULTIQC_VERSION 1.19
RUN pip install multiqc==${MULTIQC_VERSION}

# Snakemake
ENV SNAKEMAKE_VERSION 7.32.4
RUN pip install snakemake==${SNAKEMAKE_VERSION}

# Salmon
ENV SALMON_VERSION 1.10.1
WORKDIR /usr/local/src
RUN curl -LO https://github.com/COMBINE-lab/salmon/archive/refs/tags/v${SALMON_VERSION}.tar.gz && \
    tar xzf v${SALMON_VERSION}.tar.gz && \
    cd salmon-${SALMON_VERSION} && \
    cmake -DCMAKE_INSTALL_PREFIX=/usr/local/bin && \
    make
RUN make install

# Use renv for R packages
ENV RENV_CONFIG_CACHE_ENABLED FALSE
WORKDIR /usr/local/renv
COPY renv.lock renv.lock
RUN R -e "install.packages('renv')"
RUN R -e "renv::restore()" && rm -rf ~/.local/share/renv


WORKDIR /home/rstudio
