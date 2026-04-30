# Build awscli, fastp, and salmon from source in a separate image
# matching base image from https://github.com/rocker-org/rocker-versioned2/blob/master/dockerfiles/r-ver_4.5.2.Dockerfile
FROM docker.io/library/ubuntu:noble AS build

# Build dependencies
RUN apt-get update -qq \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    ca-certificates \
    cmake \
    curl \
    g++ \
    gcc \
    libboost-all-dev \
    libbz2-dev \
    libcurl4-openssl-dev \
    libdeflate-dev \
    libisal-dev \
    liblzma-dev \
    libzstd-dev \
    make \
    pkg-config \
    unzip \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /usr/local/src

# Get AWS CLI
RUN curl -o awscliv2.zip "https://awscli.amazonaws.com/awscli-exe-linux-$(arch).zip"
RUN unzip awscliv2.zip
RUN ./aws/install

# Build salmon
ARG SALMON_VERSION=1.10.3
RUN curl -LO https://github.com/COMBINE-lab/salmon/archive/refs/tags/v${SALMON_VERSION}.tar.gz
RUN tar xzf v${SALMON_VERSION}.tar.gz
RUN mkdir salmon-${SALMON_VERSION}/build
RUN cd salmon-${SALMON_VERSION}/build && \
    cmake -DCMAKE_INSTALL_PREFIX=/usr/local/salmon .. && \
    make && make install

# Build fastp
ARG FASTP_VERSION=0.23.4
RUN curl -LO https://github.com/OpenGene/fastp/archive/refs/tags/v${FASTP_VERSION}.tar.gz
RUN tar xzf v${FASTP_VERSION}.tar.gz
RUN cd fastp-${FASTP_VERSION} && \
    make && make install

# Main image with Biocconductor and other tools
FROM bioconductor/bioconductor_docker:3.22 AS final
LABEL maintainer="ccdl@alexslemonade.org"

WORKDIR /rocker-build/

# Additional dependencies for AWS runtime and handy tools
RUN apt-get update -qq \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    emacs \
    glibc-source \
    groff \
    htop \
    less \
    libisal2 \
    nano \
    vim \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# FastQC
RUN apt-get update -qq \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    fastqc \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Get rclone
RUN curl -L https://rclone.org/install.sh | bash

# Python packages
COPY requirements.txt requirements.txt
RUN pip install -r requirements.txt --break-system-packages

# Use renv for R packages
WORKDIR /usr/local/renv
COPY renv.lock renv.lock
ENV RENV_CONFIG_CACHE_ENABLED=FALSE
ENV RENV_CONFIG_INSTALL_STAGED=FALSE
RUN Rscript -e "install.packages('renv')"
ENV RENV_CONFIG_INSTALL_STAGED=FALSE
RUN Rscript - <<'RSCRIPT_EOF'
arch <- R.version[['arch']]
RUN Rscript -e "renv::restore()" \
# set up repos for both BioC and CRAN
repos <- c(
    BioCsoft = 'https://bioconductor.org/packages/3.22/bioc',
    BioCann = 'https://bioconductor.org/packages/3.22/data/annotation',
    BioCexp = 'https://bioconductor.org/packages/3.22/data/experiment',
    BioCworkflows = 'https://bioconductor.org/packages/3.22/workflows',
    BioCbooks = 'https://bioconductor.org/packages/3.22/books',
    CRAN = 'https://packagemanager.posit.co/cran/__linux__/noble/latest'
)
# add binary repo for Bioc on x86_64
if (arch == 'x86_64') {
    repos <- c(repos, BioCcontainers = 'https://bioconductor.org/packages/3.22/container-binaries/bioconductor_docker')
}
options(repos = repos)
install.packages('renv')
renv::restore(repos = repos)

# clean up
unlink("~/.cache/R/renv", recursive=TRUE)
unlink("/tmp/downloaded_packages", recursive=TRUE)
unlink(list.files(tempdir(), pattern = "Rtmp", full.names = TRUE), recursive=TRUE)
RSCRIPT_EOF


# copy aws, salmon, and fastp binaries from the build image
COPY --from=build /usr/local/aws-cli/ /usr/local/aws-cli/
RUN ln -s /usr/local/aws-cli/v2/current/bin/aws /usr/local/bin/aws
RUN ln -s /usr/local/aws-cli/v2/current/bin/aws_completer /usr/local/bin/aws_completer
COPY --from=build /usr/local/salmon/ /usr/local/
COPY --from=build /usr/local/bin/fastp /usr/local/bin/fastp

# Create the skel directory by copying in the repository contents as a template
# (limited by the .dockerignore file)
# Then run the setup script to update the default skel directory
ARG template_dir=/etc/skel-template/training-modules
COPY . ${template_dir}
RUN python3 ${template_dir}/scripts/setup-skel.py \
    --base-dir ${template_dir} \
    --skel-dir /etc/skel \
    --module-file ${template_dir}/current-modules.json

WORKDIR /home/rstudio

