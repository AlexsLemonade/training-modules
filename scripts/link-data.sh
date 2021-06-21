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
mkdir -p RNA-seq/data/leukemia
mkdir -p RNA-seq/data/medulloblastoma
mkdir -p RNA-seq/data/zebrafish-cortisol
mkdir -p RNA-seq/QC/gastric-cancer/fastp
mkdir -p RNA-seq/QC/gastric-cancer/fastqc
mkdir -p RNA-seq/data/open-pbta

# scRNA-seq module directories
mkdir -p scRNA-seq/analysis/mouse-liver/markers
mkdir -p scRNA-seq/data/glioblastoma
mkdir -p scRNA-seq/data/tabula-muris/alevin-quant
mkdir -p scRNA-seq/data/tabula-muris/normalized
mkdir -p scRNA-seq/data/hodgkins
mkdir -p scRNA-seq/data/mouse-liver/normalized
mkdir -p scRNA-seq/gene-sets

# Machine learning module directory
mkdir -p machine-learning/data

# Pathway analysis module directory
mkdir -p pathway-analysis/data

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
  RNA-seq/data/leukemia/SRP049821_metadata.tsv
  RNA-seq/data/leukemia/txi
  RNA-seq/data/medulloblastoma/SRP150101_metadata.tsv
  RNA-seq/data/medulloblastoma/txi
  RNA-seq/data/zebrafish-cortisol/zebrafish-cortisol_metadata.tsv
  RNA-seq/data/zebrafish-cortisol/salmon_quant
  RNA-seq/data/open-pbta/pbta-histologies-subset.tsv 
  RNA-seq/data/open-pbta/pbta-rsem-expected_count-subset.rds
  scRNA-seq/analysis/mouse-liver/markers
  scRNA-seq/data/glioblastoma/preprocessed
  scRNA-seq/data/hodgkins/cellranger
  scRNA-seq/data/mouse-liver/normalized
  scRNA-seq/data/hodgkins/hs_mitochondrial_genes.tsv
  scRNA-seq/data/tabula-muris/alevin-quant/10X_P4_3
  scRNA-seq/data/tabula-muris/alevin-quant/10X_P7_12
  scRNA-seq/data/tabula-muris/fastq
  scRNA-seq/data/tabula-muris/normalized/TM_normalized.rds
  scRNA-seq/data/tabula-muris/mm_mitochondrial_genes.tsv
  scRNA-seq/data/tabula-muris/mm_ensdb95_tx2gene.tsv
  scRNA-seq/gene-sets
  machine-learning/data/open-pbta
  pathway-analysis/data/leukemia
  pathway-analysis/data/medulloblastoma
  pathway-analysis/data/open-pbta
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


# Link the human indices
mkdir -p RNA-seq/index/Homo_sapiens/
hs_index_dest=RNA-seq/index/Homo_sapiens/short_index
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

# Link mouse indices to single cell
mkdir -p scRNA-seq/index/Mus_musculus/
mm_index_dest=scRNA-seq/index/Mus_musculus/short_index
mm_index_source=${share_base}/reference/refgenie/mm10_cdna/salmon_index/short
if [[ -L ${mm_index_dest} || ! -e ${mm_index_dest} ]]
then
  ln -nsf ${mm_index_source} ${mm_index_dest} 
else
  echo "${mm_index_dest} already exists and is not a link, delete or move it to create a link."
fi

mm_tx2gene_dest=scRNA-seq/index/Mus_musculus/Mus_musculus.GRCm38.95.versioned_tx2gene.tsv
mm_tx2gene_source=${share_base}/reference/tx2gene/Mus_musculus.GRCm38.95.versioned_tx2gene.tsv
if [[ -L ${mm_tx2gene_dest} || ! -e ${mm_tx2gene_dest} ]]
then
  ln -sf ${mm_tx2gene_source} ${mm_tx2gene_dest}
else
  echo "${mm_tx2gene_dest} already exists and is not a link, delete or move it to create a link."

fi

# Link the zebrafish indices
mkdir -p RNA-seq/index/Danio_rerio/
dr_index_dest=RNA-seq/index/Danio_rerio/short_index
dr_index_source=${share_base}/reference/refgenie/z11_cdna/salmon_index/short
if [[ -L ${dr_index_dest} || ! -e ${dr_index_dest} ]]
then
  ln -nsf ${dr_index_source} ${dr_index_dest} 
else
  echo "${dr_index_dest} already exists and is not a link, delete or move it to create a link."
fi
