#!/bin/bash

set -euo pipefail

# Profile to use with the download script
PROFILE=${PROFILE:-openscpca}
# Release to download
RELEASE=${RELEASE:-2024-11-25}

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Define ScPCA project id
project_id="SCPCP000006"

# Set up directories
wilms_data_dir="../../data/wilms-tumor"
processed_dir="${wilms_data_dir}/processed"

# Create directories if they don't exist yet
mkdir -p "${processed_dir}"

# Get download data script from OpenScPCA for convenience
curl -O \
    -L https://raw.githubusercontent.com/AlexsLemonade/OpenScPCA-analysis/a3d8a2c9144e8edb3894a7beeb89cdc6c3e6d681/download-data.py
# Make executable
chmod +x download-data.py

# Download Wilms tumor samples
./download-data.py \
    --samples 'SCPCS000203' \
    --format SCE \
    --release ${RELEASE} \
    --data-dir ${wilms_data_dir} \
    --profile ${PROFILE}

# Remove existing files from processed directory
if [ -z "$( ls -A ${processed_dir} )" ]; then
   echo "No processed files yet!"
else
   rm -r ${processed_dir}/*
fi

# Move processed files from release folder into processed_dir
mv ${wilms_data_dir}/${RELEASE}/${project_id}/* ${processed_dir}

# Clean up download data script
rm download-data.py
# Clean up the remnants of download structure
rm -r ${wilms_data_dir}/${RELEASE}
