# Pre-processing single cell RNA-seq data **CCDL 2019**

#### In this section, we will be running through the basics of pre-processing single-cell rna-seq data. 

For the purposes of this tutorial, we'll summarize single-cell technologies as 
being one of two groups, based on their capture methods and quantitative nature.
Depending on how the cells are sorted and what technology is used, the pre-processing steps are a bit different and the biases to look out for in post-processing also vary. 

For a more extensive run through on single-cell technologies, 
Kiselev at al have very [good tutorial for scRNA-seq](https://hemberg-lab.github.io/scRNA.seq.course/introduction-to-single-cell-rna-seq.html#experimental-methods). 
  
### 1) Full-length scRNA-seq
*Examples: Smart-seq2*  
Cells are physically separated generally into individual wells in a plate and 
often also sorted by other means (eg. Fluorescence Activated Cell Sorting). 
Each cell is then sequenced individual and has it's own fastq file (this will be two fastq files if this is paired-end sequencing.) 
The data pre-processing steps for these types of scRNA-seq data is/can be more similar to 
bulk rna-seq methods.

#### Pros:  
- Can be paired end sequencing which has less risk for 3' bias.  
- More complete coverage of transcripts, which may be better for transcript 
discovery purposes.   
  
#### Cons:  
- Is not very efficient (generally 96 cells per plate)  
- Takes much longer to run.  
- Is a lot more expensive.  
  
### 2) Tag-based scRNA-seq
*Examples: 10X Genomics, Drop-seq*  
Cells are separated by emulsion/droplets, and individual cells are given barcodes. 
Everything is then sequenced. 
These types of methods, because they are newer, are more likely to have Unique 
Molecular Identifiers (UMIs) which allow you to better control for PCR amplification errors and biases.
Individual samples have two fastq files: one for the cell barcodes
and another with the individual reads. 

#### Pros:  
- Can run potentially millions of cells at once.   
- A lot quicker computing wise.  
- Won't take up all your computer's storage.  
- A lot cheaper  
  
#### Cons:  
- Sequencing is not bidirectional so data will likely have more intense 3' bias.  
- Coverage of these technologies generally is not as deep.  

### Processing scRNA-seq data

#### Step 0: Obtaining the data
10X Genomics has plenty of [datasets](https://support.10xgenomics.com/single-cell-gene-expression/datasets) 
available to play around with. 
For this example, we will use the [1k pbmc](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_1k_v2) 
(peripheral blood mononuclear cells). 
We will use one of these samples' subsetted fastq file of this as an example.  
  
*Note*: depending on the state of the data you are working with, ie. if you have 
a `.bcl` file, you will need to use `CellRanger` with their `mkfastq` command to 
make fastq files from this.
However, because most publicly available data is in fastq format and any
core you might be working with will also likely provide you with fastq files, 
we are starting from this point.  

#### Step 1: Set up your directories
As described in a previous module, we'll make some directories from command line
for us to store these data in. 
``` 
mkdir data
mkdir ref_files
mkdir results
mkdir data/alevin_output
```

#### Step 2: Index the human transcriptome with Salmon
Before you can quantify with salmon and alevin, we need a transcriptome to be indexed.
You can use the same trancriptome index as was used for bulk-rna-seq, however,
due to the smaller pieces and amounts of single cell RNA-seq as opposed to bulk, 
you may want to build the index with a smaller `-k`. 
In this instance, we used a `-k` of 23 using the ensemble transcriptome. 
```
salmon --threads=16 --no-version-check index \
-t ref_files/Homo_sapiens.GRCh38.cdna.all.fa.gz \
-i ref_files/human_index \
-k 23
```

#### Step 3: For each sample, run Salmon with Alevin for quantification
From the command line, running Alevin is not too much different from running 
Salmon for bulk rna-seq. You'll recognize a lot of these options as the same.
In this instance, files with `R1` contain the barcodes for the cells as well as 
the Unique Molecular Identifiers while `R2` files contain the full reads for that sample.  
The file `genes_2_tx.tsv` is a gene to transcript key created based on the ensembl transcriptome. 
Note that for downstream analyses, we will also want to use the `--dumpCsvCounts` option 
as well as the `--dumpFeatures` option. 
  
`--dumpCsvCounts` will make the counts in a csv file and make it easier for us to 
import into R later.  
`--dumpFeatures` will print out information that we will need for quality checks
later on, including files with information on the UMI's and cell barcodes ("CB").  
  
See the [Alevin](https://salmon.readthedocs.io/en/latest/alevin.html) 
documentation and tutorials for a complete list of the Alevin options.  

```
salmon alevin -l ISR \
  -i ../ref_files/human_index \
  -1 pbmc_1k_v2_S1_L001_R1_subset.fastq.gz \
  -2 pbmc_1k_v2_S1_L001_R2_subset.fastq.gz \
  --chromium  \
  -p 10 \
  -o ../alevin_output/pbmc_1k_v2_S1_L001_subset \
  --tgMap ../genes_2_tx.tsv \
  --dumpCsvCounts \
  --dumpFeatures
```
