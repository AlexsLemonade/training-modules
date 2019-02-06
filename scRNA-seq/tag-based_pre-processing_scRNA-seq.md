# Pre-processing single cell RNA-seq data 

**CCDL 2019**

#### In this section, we will be running through the basics of pre-processing
#### single-cell RNA-seq data.

For the purposes of this tutorial, we'll summarize single-cell technologies as
being one of two groups, based on their capture methods and quantitative nature.
Depending on how the cells are sorted and what technology is used, the pre-processing steps are a bit different and the biases to look out for in post-processing also vary.

For a more extensive run through on single-cell technologies,
Kiselev et al have very [good tutorial for scRNA-seq](https://hemberg-lab.github.io/scRNA.seq.course/introduction-to-single-cell-rna-seq.html#experimental-methods).

### 1) Non-Tag-Based scRNA-seq  
*Example:* Smart-seq2 [(Picelli et al, 2014)](https://www.nature.com/articles/nprot.2014.006)   
Cells are physically separated generally into individual wells in a plate and
often also sorted by other means (eg. Fluorescence Activated Cell Sorting).
Each cell is then sequenced individual and has it's own fastq file (this will be two fastq files if this is paired-end sequencing.)
The data pre-processing steps for these types of scRNA-seq data is/can be more similar to
bulk RNA-seq methods.

#### Pros:  
- Can be paired end sequencing which has less risk for 3' bias.  
- More complete coverage of transcripts, which may be better for transcript
discovery purposes.   

#### Cons:  
- Is not very efficient (generally 96 cells per plate).  
- Takes much longer to run (Days/weeks depending on sample size).
- Is a lot more expensive.  

### 2) Tag-based scRNA-seq  
*Example:* 10X Genomics [(Zheng et al, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/28091601)  
Cells are separated by emulsion/droplets, and individual cells are given barcodes.
Everything is then sequenced.
These types of methods, because they are newer, are more likely to have
[Unique Molecular Identifiers (UMIs)](http://www.nature.com/doifinder/10.1038/nmeth.2772)
which allow you to better control for PCR amplification errors and biases.
Individual samples have two fastq files: one for the cell barcodes
and another with the individual reads.

#### Pros:  
- Can run potentially millions of cells at once.   
- A lot quicker computing wise.  
- Won't take up all your computer's storage.  
- A lot cheaper.  

#### Cons:  
- Sequencing is not bidirectional so data will likely have more intense 3' bias.  
- Coverage of these technologies generally is not as deep.  

*More sources on the comparisons and explanations of these technologies:*   
- [Zhang et al, 2018](https://doi.org/10.1016/j.molcel.2018.10.020)  
- [AlJanahi et al, 2018](https://doi.org/10.1016/j.omtm.2018.07.003)  
- [Angerer et al, 2017](http://dx.doi.org/10.1016/j.coisb.2017.07.004)  
- [Baran-Gale et al, 2018](https://doi.org/10.1093/bfgp/elx035)  

## Steps for Processing scRNA-seq Data:

### Step 0: Obtaining the data
10X Genomics has plenty of [datasets](https://support.10xgenomics.com/single-cell-gene-expression/datasets)
available to play around with.
For this example, we will use the [1k pbmc](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_1k_v2)
(peripheral blood mononuclear cells).
We will process a fastq file of this as an example.  
For practical purposes as far as time constraints, we have subset this fastq file
and have it ready for you to process but with only the first 3 million reads.

*Note*: depending on the state of the data you are working with, ie. if you have
a `.bcl` file, you will need to use `CellRanger` with their `mkfastq` command to
make fastq files from this.
However, because most publicly available data is in fastq format and any
core you might be working with will also likely provide you with fastq files,
we are starting from this point.  

### Step 1: Set up your output directory
As described in a previous module, we'll make a directory from command line
for us to store our output data in.
```
mkdir alevin_output  
```

### Step 2: Index the human transcriptome with Salmon
Before you can quantify with Salmon and
[Alevin](https://www.biorxiv.org/content/10.1101/335000v2), we need a transcriptome
to be indexed.
You can use the same trancriptome index as was used for bulk-RNA-seq, however,
due to the smaller pieces and amounts of single cell RNA-seq as opposed to bulk,
you may want to build the index with a smaller `-k`.
In this instance, we used a `-k` of 23 using the ensemble transcriptome.

In the interest of time, we have already run the command below and have the index
built for you.
But for your own reference, here is how you'd do it yourself:
```
salmon --threads=16 --no-version-check index \
  -t Homo_sapiens.GRCh38.cdna.all.fa.gz \
  -i human_index \
  -k 23
```

### Step 3: For each sample, run [Alevin](https://www.biorxiv.org/content/10.1101/335000v2)
for quantification
From the command line, running Alevin is not too much different from running
Salmon for bulk RNA-seq. You'll recognize a lot of these options as the same.
In this instance, files with `R1` contain the barcodes for the cells as well as
the Unique Molecular Identifiers while `R2` files contain the full reads for that sample.  

#### `-l`
The ISR library type is what is recommended for single cell data quant.

#### `--chromium`
Because we are using 10X chromium data, we have to use this flag. However,
Drop-seq data is also supported, and in this case you would use a `--dropseq`
flag instead of this.

#### `--tgMap`
This is needed to supply a trasncript to gene key that Alevin will use to
quantify the genes. For our example, we've premade the file `genes_2_tx.tsv` from
the ensembl transcriptome that we indexed above. The file has to be a .tsv file.

#### `--dumpCsvCounts`
Using this option will make the counts in a csv file and make it easier for us to
import into R later.

#### `--dumpFeatures`
This option will print out information that we will need for quality checks
later on, including files with information on the UMI's and cell barcodes ("CB").  

#### Running Salmon Alevin
Copy and paste this in your command line to run Alevin quantification
```
salmon alevin -l ISR \
  -i data/human_index \
  -1 data/pbmc_1k_v2_S1_L001_R1_subset.fastq.gz \
  -2 data/pbmc_1k_v2_S1_L001_R2_subset.fastq.gz \
  --chromium  \
  -p 10 \
  -o alevin_output \
  --tgMap data/genes_2_tx.tsv \
  --dumpCsvCounts \
  --dumpFeatures
```

See the [Alevin documentation](https://salmon.readthedocs.io/en/latest/alevin.html)
for a complete list of the Alevin options and see the
[Alevin tutorials](https://combine-lab.github.io/alevin-tutorial/2018/running-alevin/)
for example analyses.

### Step 4: Perform QC checks with alevinQC
In order to perform quality control checks, we will need to open R.
Alevin provides count data output for each transcript and cell. To read this 
data into R, we will import a function from the script "ReadAlevin.R" which is 
located in the "scripts" folder.
```r 
# Import the function to read alevin output data
source(file.path("scripts", "ReadAlevin.R"))

# Read in the data
alevin_file <- ReadAlevin("alevin_output")
```
Now that our data is imported into the R environment, we can run quality checks
using alevinQC package.
This will provide html output with graphs evaluating the data. Keep in mind that 
the data we have processed in this workshop is not the full dataset, and won't 
look as good as the full set.
```r 
# Produce a QC report
alevinQC::alevinQCReport(alevin_file,
                         sampleId = "pbmc_1k_v2_S1_L001_subset", 
                         outputFile = "pbmc_1k_v2_S1_L001_qc_report.html", 
                         outputDir = "data",
                         outputFormat = "html_document")
```
Now you can check out "pbmc_1k_v2_S1_L001_qc_report.html" in order to examine 
the quality of your data and performance of Alevin. 
