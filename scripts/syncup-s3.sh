#! /bin/bash

# This script is used to sync files required for training module notebooks to S3
# To run this script, you must have write permission to the destination bucket 
# and have set up credentials with `aws configure`

# As new directories and files are required, they should be added to the
# the `sync_dirs` or `sync_files` arrays below, followed by running this script

set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"
# move back up to the training modules root
cd ..

# destination bucket
bucket=s3://ccdl-training-data/training-modules

# Directories and files are listed separately so that we can take advantage
# of the ability of `aws s3 sync` to avoid copying files inside directories that 
# are already present on S3.  
# Unfortunately, `sync` does not work on individual files, so we have to handle 
# them separately. While we could presumably check individual files for changes
# before upload as well, it probably isn't worth it.
sync_dirs=(
  intro-to-R-tidyverse/data
  RNA-seq/data/gastric-cancer/salmon_quant
  RNA-seq/data/NB-cell/tximport
  scRNA-seq/data/tabula-muris/alevin/10X_P4_3
  scRNA-seq/index/Mus_musculus
  machine-learning/data/open-pbta/processed
)

sync_files=(
  RNA-seq/data/gastric-cancer/gastric-cancer_metadata.tsv
  RNA-seq/data/NB-cell/NB-cell_metadata.tsv
  RNA-seq/data/leukemia/SRP049821_metadata.tsv
  RNA-seq/index/Homo_sapiens/Homo_sapiens.GRCh38.95_tx2gene.tsv
  scRNA-seq/data/glioblastoma/preprocessed/tximport/count_matrix.tsv
  scRNA-seq/data/glioblastoma/preprocessed/darmanis_metadata.tsv
  scRNA-seq/data/tabula-muris/normalized/TM_normalized.rds
  pathway-analysis/data/medulloblastoma/medulloblastoma_vst_collapsed.tsv
)

for loc in ${sync_dirs[@]}
do
  # upload directories to S3, make public, ignore timestamps, ignore hidden files
  aws s3 sync ${loc} ${bucket}/${loc} \
   --acl public-read --size-only \
   --exclude ".*"
done

for loc in ${sync_files[@]}
do
  # upload individual files to S3, make public
  aws s3 cp ${loc} ${bucket}/${loc} \
   --acl public-read
done
