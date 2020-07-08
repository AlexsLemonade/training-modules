#!/bin/bash

set -e
set -o pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Use the OpenPBTA bucket as the default.
URL=${OPENPBTA_URL:-https://s3.amazonaws.com/kf-openaccess-us-east-1-prd-pbta/data}
RELEASE=${OPENPBTA_RELEASE:-release-v16-20200320}

# All of the OpenPBTA data will be located in data/open-pbta folder (from root
# of module directory) - we won't include the release in the folder structure
data_dir=../data/open-pbta

# The md5sum file provides our single point of truth for which files are in a release.
curl --create-dirs $URL/$RELEASE/md5sum.txt -o ${data_dir}/md5sum.txt -z ${data_dir}/md5sum.txt

# Consider the filenames in the md5sum file and the release notes
FILES=(`tr -s ' ' < ${data_dir}/md5sum.txt | cut -d ' ' -f 2` release-notes.md)

# Download the items in FILES if newer than what's on server
for file in "${FILES[@]}"
do
  curl --create-dirs $URL/$RELEASE/$file -o ${data_dir}/$file -z ${data_dir}/$file
done
