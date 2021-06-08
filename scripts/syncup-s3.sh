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
  RNA-seq/data/NB-cell/txi
  RNA-seq/data/leukemia/txi
  RNA-seq/data/medulloblastoma/txi
  RNA-seq/data/open-pbta
  scRNA-seq/data/hodgkins/cellranger
  scRNA-seq/data/tabula-muris/alevin-quant/10X_P4_3
  scRNA-seq/data/tabula-muris/alevin-quant/10X_P7_12
  scRNA-seq/index/Mus_musculus
  machine-learning/data/open-pbta/processed
  pathway-analysis/data/leukemia
  pathway-analysis/data/medulloblastoma
  pathway-analysis/data/open-pbta
)

sync_files=(
  RNA-seq/data/gastric-cancer/gastric-cancer_metadata.tsv
  RNA-seq/data/NB-cell/NB-cell_metadata.tsv
  RNA-seq/data/leukemia/SRP049821_metadata.tsv
  RNA-seq/data/medulloblastoma/SRP150101_metadata.tsv
  scRNA-seq/data/glioblastoma/preprocessed/txi/count_matrix.tsv
  scRNA-seq/data/glioblastoma/preprocessed/darmanis_metadata.tsv
  scRNA-seq/data/hodgkins/hs_mitochondrial_genes.tsv
  scRNA-seq/data/tabula-muris/normalized/TM_normalized.rds
  scRNA-seq/data/tabula-muris/mm_mitochondrial_genes.tsv
  scRNA-seq/data/tabula-muris/mm_ensdb95_tx2gene.tsv
)

for loc in ${sync_dirs[@]}
do
  if [[ -d ${loc} ]]; then
    # upload directories to S3, make public, ignore timestamps, ignore hidden files
    aws s3 sync ${loc} ${bucket}/${loc} \
    --acl public-read --size-only \
    --exclude ".*"
  else
    echo "${loc} does not exist."
  fi
done

echo "Directories synced"

for loc in ${sync_files[@]}
do
  if [[ -f ${loc} ]]; then
    # upload individual files to S3, make public
    aws s3 cp ${loc} ${bucket}/${loc} \
    --acl public-read
  else
    echo "${loc} does not exist."
  fi
done

echo "Files uploaded"
