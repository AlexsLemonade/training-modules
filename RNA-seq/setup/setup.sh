#!/bin/bash

set -euo pipefail

jobs=8

# create gastric cancer files (keep fastq)
snakemake -j ${jobs} --notemp

# create NB cell line files(fastq will be deleted)
snakemake -j ${jobs} --configfile config-NB.yaml