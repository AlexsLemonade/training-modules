#! /bin/bash

share_base=/shared/data/training-modules/
# RNA-seq module
mkdir -p RNA-seq/data/gastric-cancer
mkdir -p RNA-seq/data/gastric-cancer/salmon_quant
mkdir -p RNA-seq/data/NB-cell
mkdir -p RNA-seq/QC/gastric-cancer/fastp/
mkdir -p RNA-seq/QC/gastric-cancer/fastqc/



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
  RNA-seq/QC/gastric-cancer/fastp/SRR585572
  RNA-seq/QC/gastric-cancer/fastp/SRR585573
  RNA-seq/QC/gastric-cancer/fastp/SRR585574
  RNA-seq/QC/gastric-cancer/fastp/SRR585575
  RNA-seq/QC/gastric-cancer/fastp/SRR585576
  RNA-seq/QC/gastric-cancer/fastp/SRR585577
  RNA-seq/QC/gastric-cancer/fastqc/SRR585571
  RNA-seq/QC/gastric-cancer/fastqc/SRR585572
  RNA-seq/QC/gastric-cancer/fastqc/SRR585573
  RNA-seq/QC/gastric-cancer/fastqc/SRR585574
  RNA-seq/QC/gastric-cancer/fastqc/SRR585575
  RNA-seq/QC/gastric-cancer/fastqc/SRR585576
  RNA-seq/QC/gastric-cancer/fastqc/SRR585577
  RNA-seq/data/NB-cell/NB-cell_metadata.tsv
  RNA-seq/data/NB-cell/salmon_quant
)
for loc in ${link_locs[@]}
do
  ln -nsf ${share_base}/${loc} ${loc}
done

## link indexes
mkdir -p RNA-seq/index/Homo_sapiens/
ln -nsf /shared/data/reference/refgenie/hg38_cdna/salmon_index/short RNA-seq/index/Homo_sapiens/short_index 
ln -sf /shared/data/reference/tx2gene/Homo_sapiens.GRCh38.95_tx2gene.tsv RNA-seq/index/Homo_sapiens/Homo_sapiens.GRCh38.95_tx2gene.tsv
