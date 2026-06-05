#!/bin/bash

# This script is to patch the RStudio server in the Docker image to fix hangs
# when dealing with large SCE files (hopefully)
#
# We will be patching a line from the function isSerializableImpl() in the file SessionEnvironment.R
# isSerializableImpl() in the file SessionEnvironment.R:
# This function has previously been implicated in a similar pegged CPU bug to the one we
# have encountered, which was fixed with the addition of igraph in the chunk below:
#
#   # Check for 'known-safe' object classes.
#   if (inherits(value, c("data.frame", "igraph")))
#      return(TRUE)
#
# There is only one instance of this line in the file, so we can be confident that we have patched the correct line.
# If the file changes, then we will fail on the MD5sum check.

set -euo pipefail

# Check that the file exists and has the expected hash
FILE="/usr/lib/rstudio-server/R/modules/SessionEnvironment.R"
MD5="69d4d879c71d13daf6ab4996288921ef"

if [[ ! -f "$FILE" ]]; then
    echo "Error: $FILE does not exist."
    exit 1
fi
if [[ "$(md5sum "$FILE" | awk '{print $1}')" != "$MD5" ]]; then
    echo "Error: $FILE has an unexpected hash."
    exit 1
fi

# Patch the file to add SummarizedExperiment (includes SingleCellExperiment) and Matrix
# to this list of types known to be serializable

sed -i.bak 's/inherits(value, c("data.frame", "igraph"))/inherits(value, c("data.frame", "igraph", "SummarizedExperiment", "Matrix"))/' "$FILE"
