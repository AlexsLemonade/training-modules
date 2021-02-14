#! /bin/bash

# This script is used to sync files required for training module notebooks to S3

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"
# move back up to the training modules root
cd ..

bucket=s3://ccdatalab-training-data/training-modules

sync_locs=(
  intro-to-R-tidyverse/data
  RNA-seq/data/gastric-cancer/gastric-cancer_metadata.tsv
  RNA-seq/data/gastric-cancer/salmon_quant
  RNA-seq/data/NB-cell/NB-cell_metadata.tsv
  RNA-seq/data/NB-cell/tximport
  RNA-seq/index/Homo_sapiens/Homo_sapiens.GRCh38.95_tx2gene.tsv
  machine-learning/data/open-pbta/processed
)

for loc in ${link_locs[@]}
do
  # upload to S3, make public, ignore timestamps, ignore hidden files
  aws s3 sync ${loc} ${bucket}/${loc} --acl public-read --size-only --exclude ".*"
done
