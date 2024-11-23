#!/bin/bash

set -euo pipefail

# Profile to use with the download script
PROFILE=${PROFILE:-openscpca}
# Release to download
RELEASE=${RELEASE:-2024-08-22}

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Set up directories
ewing_data_dir="../../data/ewing-sarcoma"
annotations_dir="${ewing_data_dir}/annotations"
processed_dir="${ewing_data_dir}/processed"

# Create directories if they don't exist yet
mkdir -p "${annotations_dir}"
mkdir -p "${processed_dir}"

# Get download data script from OpenScPCA for convenience
curl -O \
    -L https://raw.githubusercontent.com/AlexsLemonade/OpenScPCA-analysis/a3d8a2c9144e8edb3894a7beeb89cdc6c3e6d681/download-data.py
# Make executable
chmod +x download-data.py

# Download Ewing sarcoma tumor samples
./download-data.py \
    --samples 'SCPCS000490,SCPCS000492,SCPCS000493,SCPCS000494,SCPCS000495,SCPCS000496,SCPCS000749' \
    --format SCE \
    --release ${RELEASE} \
    --data-dir ${ewing_data_dir} \
    --profile ${PROFILE}

# # Download Ewing sarcoma metadata
./download-data.py \
    --projects SCPCP000015 \
    --metadata-only \
    --release ${RELEASE} \
    --data-dir ${ewing_data_dir} \
    --profile ${PROFILE}

# Remove existing files from processed directory
if [ -z "$( ls -A ${processed_dir} )" ]; then
   echo "No processed files yet!"
else
   rm -r ${processed_dir}/*
fi

# Move files from release folder
mv ${ewing_data_dir}/${RELEASE}/SCPCP000015/* ${processed_dir}
mv "${processed_dir}/single_cell_metadata.tsv" "${annotations_dir}/ewing_sarcoma_sample_metadata.tsv"

# Remove PDX samples from metadata
Rscript - << EOF

sample_metadata_df <- readr::read_tsv("${annotations_dir}/ewing_sarcoma_sample_metadata.tsv")
sample_metadata_df |>
    dplyr::filter(stringr::str_detect(sample_type, "xenograft", negate = TRUE)) |>
    readr::write_tsv("${annotations_dir}/ewing_sarcoma_sample_metadata.tsv")

EOF

# Clean up download data script
rm download-data.py
# Clean up the remnants of download structure
rm -r ${ewing_data_dir}/${RELEASE}
