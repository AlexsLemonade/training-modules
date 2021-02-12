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

# Need this package to make plots colorblind friendly
RUN Rscript -e "remotes::install_github('clauswilke/colorblindr', ref = '1ac3d4d62dad047b68bb66c06cee927a4517d678', dependencies = TRUE)"

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

# Use renv for R packages
ENV RENV_VERSION 0.12.5-2
RUN R -e "install.packages('remotes')"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

WORKDIR /usr/local/renv
COPY renv.lock renv.lock
RUN R -e 'renv::consent(provided = TRUE)'
RUN R -e 'renv::restore()'

WORKDIR /home/rstudio
