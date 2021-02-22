#!/bin/bash
# Uses the index from AlexsLemonade/training-txome-prep

set -e
set -o pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"
data_dir=../data

mkdir -p ${data_dir}/quants
for srr_iter in {1652983..1653014}
do
  salmon --no-version-check --threads=8 quant \
  	-l A -i index/Mus_musculus/short_index/ \
  	-1 ${data_dir}/fastq/SRR${srr_iter}/SRR${srr_iter}_1.fastq.gz \
  	-2 ${data_dir}/fastq/SRR${srr_iter}/SRR${srr_iter}_2.fastq.gz \
  	-o ${data_dir}/quants/SRR${srr_iter} \
  	--validateMappings --rangeFactorizationBins 4 \
  	--seqBias --gcBias
done
