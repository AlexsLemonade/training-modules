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
mkdir -p scRNA-seq/data/glioblastoma/
mkdir -p scRNA-seq/data/tabula-muris/
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
  scRNA-seq/data/glioblastoma/preprocessed
  scRNA-seq/data/tabula-muris/fastq
  scRNA-seq/data/tabula-muris/normalized/TM_normalized.rds
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



# Link the indices

index_dir_names=(
  RNA-seq/index/Homo_sapiens
  RNA-seq/index/Danio_rerio
)

index_source_names=(
  ${share_base}/reference/refgenie/hg38_cdna/salmon_index/short
  ${share_base}/reference/refgenie/z11_cdna/salmon_index/short
)

# Run through this for each dir name
for index_dir in ${index_dir_names[@]}
do 
  # Make the directory itself
  mkdir -p ${index_dir}
  
  # Declare the short_index directory and source directory
  index_dest=${index_dir}/short_index
  index_source=${index_dir_names}
  
  if [[ -L ${index_dest} || ! -e ${index_dest} ]]
  then
    ln -nsf ${index_source} ${index_dest} 
  else
    echo "${index_dest} already exists and is not a link, delete or move it to create a link."
  fi
done

hs_tx2gene_dest=RNA-seq/index/Homo_sapiens/Homo_sapiens.GRCh38.95_tx2gene.tsv
hs_tx2gene_source=${share_base}/reference/tx2gene/Homo_sapiens.GRCh38.95_tx2gene.tsv
if [[ -L ${hs_tx2gene_dest} || ! -e ${hs_tx2gene_dest} ]]
then
  ln -sf ${hs_tx2gene_source} ${hs_tx2gene_dest}
else
  echo "${hs_tx2gene_dest} already exists and is not a link, delete or move it to create a link."
fi

# Mouse tx2gene for single cell
mkdir -p scRNA-seq/index/Mus_musculus/
mm_tx2gene_dest=scRNA-seq/index/Mus_musculus/Mus_musculus.GRCm38.95.versioned_tx2gene.tsv
mm_tx2gene_source=${share_base}/reference/tx2gene/Mus_musculus.GRCm38.95.versioned_tx2gene.tsv
if [[ -L ${mm_tx2gene_dest} || ! -e ${mm_tx2gene_dest} ]]
then
  ln -sf ${mm_tx2gene_source} ${mm_tx2gene_dest}
else
  echo "${mm_tx2gene_dest} already exists and is not a link, delete or move it to create a link."

fi
