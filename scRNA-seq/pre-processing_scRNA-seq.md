# Pre-processing single cell RNA-seq data

In this section, we will be running through the basics of pre-processing single-cell rna-seq data. 

Single cell RNA-seq data can be split into two major categories based on how the cells are sorted. 
Depending on how the cells are sorted and what technology is used, the pre-processing steps are a bit different and 
the biases to look out for in postt-processing also vary. 

For more information:
Kiselev at al have very [good tutorial for scRNA-seq processing in general](https://hemberg-lab.github.io/scRNA.seq.course/introduction-to-single-cell-rna-seq.html#experimental-methods). 

### 1) Full-length
Example: Smart-seq2
Bidrectional sequencing allows you to obtain the full length sequence for your reads. 
*Pros:* more complete coverage of transcripts, which is better for transcript discovery purposes. 
*Cons:* It is not as efficient as tag based. 

### 2) Tag - based 
Example: 10X Genomics
Individual cells are given barcodes to parse out the source of sequences
*Pros:* Can run 1000s or millions of cells at once. Highly efficient.
*Cons:* Sequencing is not bidirectional so data will have more intense 3' bias.

## Obtaining the data


## Processing scRNA-seq fastq files:

### Step 1: Index the human transcriptome with Salmon
Before you can quantify with salmon and alevin, we need a transcriptome to be indexed.
```
salmon --threads=16 --no-version-check index \
-t ref_files/Homo_sapiens.GRCh38.cdna.all.fa.gz \
-i ref_files/human_index \
-k 23
```

### Step 2: For each sample, run Salmon with Alevin for quantification
In this instance, files with `R1` contain the Unique Molecular Identifiers while
`R2` files contain the full reads for that sample.  
The file `genes_2_tx.tsv` is a gene to transcript key created based on the ensembl transcriptome. 

```
salmon alevin -l ISR \
    -i ../ref_files/human_index \
    -1 pbmc_1k_v2_S1_L001_R1_001.fastq.gz \
    -2 pbmc_1k_v2_S1_L001_R2_001.fastq.gz \
    --chromium  \
    -p 10 \
    -o alevin_output \
    --tgMap ../genes_2_tx.tsv
```

### Step 3: Run QC using alevin QC
In R, we can run initial quality checks for this mapping using [csoneson/alevinQC](https://github.com/csoneson/alevinQC) package.
```r
alevinQCReport(baseDir = system.file("extdata/alevin_example", package = "alevinQC"),
               sampleId = "testSample", 
               outputFile = "alevinReport.html", 
               outputFormat = "html_document",
               outputDir = tempdir(), forceOverwrite = TRUE)"
```
