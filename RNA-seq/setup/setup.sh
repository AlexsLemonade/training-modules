#!/bin/bash

# This script runs two snakemake workflows to populate the 
# `/shared/training-modules/RNA-seq` directory with preprocessed data
# for the notebooks in that module


# safety flags
set -euo pipefail
# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

jobs=12

# create gastric cancer files (keep fastq)
snakemake -j ${jobs} --notemp --configfile config-GC.yaml

# create NB cell line files (fastq will be deleted)
snakemake -j ${jobs} --configfile config-NB.yaml

# create leukemia files
sankemake -j ${jobs} --configfile config-leuk.yaml
