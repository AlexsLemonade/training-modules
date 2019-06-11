#!/bin/bash

# CCDL 2019
# Run a single sample, SRR585570, through our bulk RNA-seq pipeline, which
# includes FastQC, light trimming with Fastp, and Salmon with mapping validation
# enabled

# FastQC -----------------------------------------------------------------------

# In the interest of time, we'll run one of the fastq files through FastQC
fastqc data/fastq/gastric_cancer/SRR585570/SRR585570_1.fastq.gz \
	-o data/QC/gastric_cancer/

# Fastp ------------------------------------------------------------------------

# Create a directory to hold the JSON and HTML output from Fastp
mkdir fastp_output

# Run the adapter and quality trimming step -- also produces QC report
fastp -i data/fastq/gastric_cancer/SRR585570/SRR585570_1.fastq.gz \
	-I data/fastq/gastric_cancer/SRR585570/SRR585570_2.fastq.gz \
	-o data/fastq/gastric_cancer/SRR585570/SRR585570_fastp_1.fastq.gz \
	-O data/fastq/gastric_cancer/SRR585570/SRR585570_fastp_2.fastq.gz \
	--qualified_quality_phred 15 \
	--length_required 20 \
	--report_title "SRR585570" \
	--json fastp_output/SRR585570_fastp.json \
	--html fastp_output/SRR585570_fastp.html

# Quantification with Salmon ---------------------------------------------------

# We perform quantification on the files that have been trimmed
# and use the index generated with -k 23, as this may "improve sensitivity"
# per the Salmon documentation
salmon quant -i index/Homo_sapiens/short_index \
	-l A \
    -1 data/fastq/gastric_cancer/SRR585570/SRR585570_fastp_1.fastq.gz \
    -2 data/fastq/gastric_cancer/SRR585570/SRR585570_fastp_2.fastq.gz \
    -o data/quant/gastric_cancer/SRR585570 \
    --validateMappings --rangeFactorization 4 \
	--gcBias --seqBias
