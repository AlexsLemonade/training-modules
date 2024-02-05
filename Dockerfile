FROM bioconductor/bioconductor_docker:RELEASE_3_18
LABEL maintainer="ccdl@alexslemonade.org"

WORKDIR /rocker-build/

RUN apt-get update -qq

RUN apt-get install -y --no-install-recommends \
    less
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
ENV SALMON_VERSION 1.10.0
WORKDIR /usr/local/src
RUN wget --quiet https://github.com/COMBINE-lab/salmon/releases/download/v${SALMON_VERSION}/salmon-${SALMON_VERSION}_linux_x86_64.tar.gz
RUN tar xzf salmon-${SALMON_VERSION}_linux_x86_64.tar.gz && \
    rm -f salmon-${SALMON_VERSION}_linux_x86_64.tar.gz && \
    ln -s /usr/local/src/salmon-latest_linux_x86_64/bin/salmon /usr/local/bin/salmon

# Use renv for R packages
ENV RENV_CONFIG_CACHE_ENABLED FALSE
RUN R -e "install.packages('renv')"

WORKDIR /usr/local/renv
COPY renv.lock renv.lock
RUN R -e 'renv::consent(provided = TRUE); renv::restore()' && rm -rf ~/.local/share/renv


WORKDIR /home/rstudio
