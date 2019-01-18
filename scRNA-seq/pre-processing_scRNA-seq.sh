#!/bin/bash
# C. Savonen
# CCDL for ALSF 2018

# Purpose: running the post-processing steps for single cell RNA-seq data.

# Data from: 

# Make a directory to store the data you download
mkdir 1kPBMC_data
mkdir ref_files

# Download the dataset
cd 1kPBMC_data
curl -O http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v2/pbmc_1k_v2_fastqs.tar
tar -xvf pbmc_1k_v2_fastqs.tar

#--------------Will run setup only if it hasn't been ran before----------------#
# Will check for genome index first before running
if [ ! -e ref_files/human_index ]; then

  # Get the human transcriptome
  curl ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz \
    -o ref_files/Homo_sapiens.GRCh38.cdna.all.fa.gz

  # Index the human transcriptome
  salmon --threads=16 --no-version-check index \
    -t ref_files/Homo_sapiens.GRCh38.cdna.all.fa.gz \
    -i ref_files/human_index \
    -k 23
fi

#http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v2/pbmc_1k_v2_molecule_info.h5
#--------------------------- Quantify samples-------------------------------#
# For each fastq file pair run salmon for quantfication
  cd pbmc_1k_v2_fastqs
  
  for f in `ls *R1_001.fastq.gz | sed 's/R1_001.fastq.gz//' `
    do
    echo "Processing sample ${f}"

    # Run Salmon
    salmon alevin -l ISR \
      -i ../ref_files/human_index \
      -1 ${f}R1_001.fastq.gz \
      -2 ${f}R2_001.fastq.gz\
      --chromium  \
      -p 10 \
      -o alevin_output_${f} \
      --tgMap ../genes_2_tx.tsv \
      –dumpCsvCounts
    done

# Run QC using alevin QC
R -e "alevinQCReport(baseDir = system.file("extdata/alevin_example", package = "alevinQC"),
                    sampleId = "testSample", 
                    outputFile = "alevinReport.html", 
                    outputFormat = "html_document",
                    outputDir = tempdir(), forceOverwrite = TRUE)"
