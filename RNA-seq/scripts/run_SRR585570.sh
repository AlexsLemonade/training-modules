#!/bin/bash

# CCDL 2021
# Run a single sample, SRR585570, through our bulk RNA-seq pipeline, which
# includes FastQC, light trimming with Fastp, and Salmon with mapping validation
# enabled

# FastQC -----------------------------------------------------------------------

# Create a directory to hold the FastQC output
mkdir -p QC/gastric-cancer/fastqc/SRR585570

# In the interest of time, we'll run one of the fastq files through FastQC
fastqc data/gastric-cancer/fastq/SRR585570/SRR585570_1.fastq.gz \
    -o QC/gastric-cancer/fastqc/SRR585570

# fastp ------------------------------------------------------------------------

# Create a directory to hold the trimmed fastq files
mkdir -p data/gastric-cancer/fastq-trimmed/SRR585570
# Create a directory to hold the QC output from Fastp
mkdir -p QC/gastric-cancer/fastp/SRR585570

# Run the adapter and quality trimming step -- also produces QC report
fastp -i data/gastric-cancer/fastq/SRR585570/SRR585570_1.fastq.gz \
    -I data/gastric-cancer/fastq/SRR585570/SRR585570_2.fastq.gz \
    -o data/gastric-cancer/fastq-trimmed/SRR585570/SRR585570_fastp_1.fastq.gz \
    -O data/gastric-cancer/fastq-trimmed/SRR585570/SRR585570_fastp_2.fastq.gz \
    --qualified_quality_phred 15 \
    --length_required 20 \
    --report_title "SRR585570" \
    --json QC/gastric-cancer/fastp/SRR585570/SRR585570_fastp.json \
    --html QC/gastric-cancer/fastp/SRR585570/SRR585570_fastp.html

# Quantification with Salmon ---------------------------------------------------

# We perform quantification on the files that have been trimmed
# and use the index generated with -k 23, as this may "improve sensitivity"
# per the Salmon documentation
salmon quant -i index/Homo_sapiens/short_index \
    -l A \
    -1 data/gastric-cancer/fastq-trimmed/SRR585570/SRR585570_fastp_1.fastq.gz \
    -2 data/gastric-cancer/fastq-trimmed/SRR585570/SRR585570_fastp_2.fastq.gz \
    -o data/gastric-cancer/salmon_quant/SRR585570 \
    --validateMappings --rangeFactorizationBins 4 \
    --gcBias --seqBias \
    --threads 4
