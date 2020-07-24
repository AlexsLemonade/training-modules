#!/bin/bash
# Adapted from:
# https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/d23b5f07f6eadeeadf7a5061a374edf17e6bf279/download-data.sh

set -e
set -o pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Use the most recent release, as of writing this script, as the default
RELEASE=${OPENPBTA_RELEASE:-release-v16-20200320}

# Use the OpenPBTA bucket.
bucket_url=https://s3.amazonaws.com/kf-openaccess-us-east-1-prd-pbta/data

# All of the OpenPBTA data will be located in data/open-pbta/download folder
# (from root of module directory) - we won't include the release in the folder structure
data_dir=../data/open-pbta/download

# The md5sum file provides our single point of truth for which files are in a release.
curl --create-dirs $bucket_url/$RELEASE/md5sum.txt -o ${data_dir}/md5sum.txt

# Consider the filenames in the md5sum file and the release notes
FILES=(`tr -s ' ' < ${data_dir}/md5sum.txt | cut -d ' ' -f 2` release-notes.md)

# Download the items in FILES if newer than what's on server
for file in "${FILES[@]}"
do
  curl --create-dirs $bucket_url/$RELEASE/$file -o ${data_dir}/$file
done

cd ${data_dir}
echo "Checking MD5 hashes..."
md5sum -c md5sum.txt

