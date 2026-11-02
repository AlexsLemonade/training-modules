# Build awscli, fastp, and salmon from source in a separate image
# matching base image from https://github.com/rocker-org/rocker-versioned2/blob/master/dockerfiles/r-ver_4.5.2.Dockerfile
FROM docker.io/library/ubuntu:noble AS build

# Build dependencies
RUN apt-get update -qq \
    && DEBIAN_FRONTEND=noninteractive \
    apt-get install -y --no-install-recommends \
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
    && apt-get clean

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
    && DEBIAN_FRONTEND=noninteractive \
    apt-get install -y --no-install-recommends \
    ca-certificates \
    emacs \
    glibc-source \
    groff \
    less \
    libisal2 \
    && apt-get clean

# FastQC
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    fastqc \
    && apt-get clean

# Get rclone
RUN curl -L https://rclone.org/install.sh | bash

# Python packages in a venv that is set as default for users
COPY requirements.txt /tmp/requirements.txt
RUN python3 -m venv /opt/venv \
    && /opt/venv/bin/pip install --no-cache-dir -r /tmp/requirements.txt \
    && rm /tmp/requirements.txt \
    && chgrp -R staff /opt/venv && chmod -R g+ws /opt/venv \
    && echo '[ "$(id -u)" -ne 0 ] && export PATH=/opt/venv/bin:$PATH' > /etc/profile.d/venv.sh \
    && cat >> "$(R RHOME)/etc/Rprofile.site" <<'EOF'
Sys.setenv(
    RETICULATE_PYTHON = "/opt/venv/bin/python",
    VIRTUAL_ENV_DISABLE_PROMPT = "1"
  )
EOF

# Use renv for R packages
WORKDIR /usr/local/renv
ENV RENV_CONFIG_CACHE_ENABLED=FALSE
RUN Rscript -e "install.packages('renv')"


COPY renv.lock renv.lock
RUN Rscript -e "options(pkgType='binary'); renv::restore(repos = c(CRAN = 'https://packagemanager.posit.co/cran/__linux__/noble/latest'))" \
    && rm -rf ~/.cache/R/renv \
    && rm -rf /tmp/downloaded_packages \
    && rm -rf /tmp/Rtmp*

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

