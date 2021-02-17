#! /bin/bash

# This script is used to sync files required for training module notebooks to S3

set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"
# move back up to the training modules root
cd ..

bucket=s3://ccdl-training-data/training-modules

sync_dirs=(
  intro-to-R-tidyverse/data
  RNA-seq/data/gastric-cancer/salmon_quant
  RNA-seq/data/NB-cell/tximport
  scRNA-seq/data/tabula-muris/alevin-quant/10X_P4_3
  scRNA-seq/index/Mus_musculus
  machine-learning/data/open-pbta/processed
)

sync_files=(
  RNA-seq/data/gastric-cancer/gastric-cancer_metadata.tsv
  RNA-seq/data/NB-cell/NB-cell_metadata.tsv
  RNA-seq/index/Homo_sapiens/Homo_sapiens.GRCh38.95_tx2gene.tsv
  scRNA-seq/data/glioblastoma/preprocessed/tximport/count_matrix.tsv
  scRNA-seq/data/glioblastoma/preprocessed/darmanis_metadata.tsv
  scRNA-seq/data/tabula-muris/normalized/TM_normalized.rds
  pathway-analysis/data/medulloblastoma/medulloblastoma_vst_collapsed.tsv
)

for loc in ${sync_dirs[@]}
do
  # upload to S3, make public, ignore timestamps, ignore hidden files
  aws s3 sync ${loc} ${bucket}/${loc} \
   --acl public-read --size-only \
   --exclude ".*"
done

for loc in ${sync_files[@]}
do
  # upload to S3, make public
  aws s3 cp ${loc} ${bucket}/${loc} \
   --acl public-read
done
