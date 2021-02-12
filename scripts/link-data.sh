#! /bin/bash

# This script is used to establish symlinks for training modules to use data
# stored in a shared directory.

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"
# move back up to the training modules root
cd ..

# location for the data modules
share_base=/shared/data
modules_base=${share_base}/training-modules


# RNA-seq module directories
mkdir -p RNA-seq/data/gastric-cancer
mkdir -p RNA-seq/data/gastric-cancer/salmon_quant
mkdir -p RNA-seq/data/NB-cell
mkdir -p RNA-seq/QC/gastric-cancer/fastp/
mkdir -p RNA-seq/QC/gastric-cancer/fastqc/
mkdir -p machine-learning/data



link_locs=(
  RNA-seq/data/gastric-cancer/gastric-cancer_metadata.tsv
  RNA-seq/data/gastric-cancer/fastq
  RNA-seq/data/gastric-cancer/salmon_quant/SRR585571
  RNA-seq/data/gastric-cancer/salmon_quant/SRR585572
  RNA-seq/data/gastric-cancer/salmon_quant/SRR585573
  RNA-seq/data/gastric-cancer/salmon_quant/SRR585574
  RNA-seq/data/gastric-cancer/salmon_quant/SRR585575
  RNA-seq/data/gastric-cancer/salmon_quant/SRR585576
  RNA-seq/data/gastric-cancer/salmon_quant/SRR585577
  RNA-seq/QC/gastric-cancer/fastp/SRR585571
  RNA-seq/QC/gastric-cancer/fastqc/SRR585571
  RNA-seq/data/NB-cell/NB-cell_metadata.tsv
  RNA-seq/data/NB-cell/salmon_quant
  machine-learning/data/open-pbta
)
for loc in ${link_locs[@]}
do
  # only make the links if replacing an old link or the file doesn't exist
  if [[ -L ${loc} || ! -e ${loc} ]]
  then
    ln -nsf ${modules_base}/${loc} ${loc}
  else
    echo "${loc} already exists and is not a link, delete or move it to create a link."
  fi
done

## link indexes
mkdir -p RNA-seq/index/Homo_sapiens/
hs_index_dest=RNA-seq/index/Homo_sapiens/Homo_sapiens.GRCh38.95_tx2gene.tsv
hs_index_source=${share_base}/reference/refgenie/hg38_cdna/salmon_index/short
if [[ -L ${hs_index_dest} || ! -e ${hs_index_dest} ]]
then
  ln -nsf ${hs_index_source} ${hs_index_dest} 
else
  echo "${hs_index_dest} already exists and is not a link, delete or move it to create a link."
fi

hs_tx2gene_dest=RNA-seq/index/Homo_sapiens/Homo_sapiens.GRCh38.95_tx2gene.tsv
hs_tx2gene_source=${share_base}/reference/tx2gene/Homo_sapiens.GRCh38.95_tx2gene.tsv
if [[ -L ${hs_tx2gene_dest} || ! -e ${hs_tx2gene_dest} ]]
then
  ln -sf ${hs_tx2gene_source} ${hs_tx2gene_dest}
else
  echo "${hs_tx2gene_dest} already exists and is not a link, delete or move it to create a link."

fi
