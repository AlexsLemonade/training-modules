# Pre-processing single-cell RNA-seq data 

**CCDL 2019**

#### In this section, we will be running through the basics of pre-processing
single-cell RNA-seq data.

We will be using a peripheral blood mononuclear cells (pbmc's) 10X Genomics [(Zheng et al, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/28091601) dataset as an example.
For 10X Genomics scRNA-seq data, cells are separated by emulsion/droplets, and 
individual cells are given barcodes.
These data also have
[Unique Molecular Identifiers (UMIs)](http://www.nature.com/doifinder/10.1038/nmeth.2772)
which allow us to better control for PCR amplification errors and biases.

Note: Raw single-cell RNA-seq data from non-tag-based methods, like the Smart-seq2 
dataset we were working with in the previous section, can be processed using 
Salmon, just like was done in the bulk RNA-seq module. 

## Steps for Processing scRNA-seq Data:

### Step 0: Obtaining the data

10X Genomics has plenty of [datasets](https://support.10xgenomics.com/single-cell-gene-expression/datasets)
available to play around with.
For this example, we will use the [1k pbmc](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_1k_v2).
We will process a fastq file of this as an example.  
For practical purposes as far as time constraints, we have subset this fastq file
to the first 3 million reads.

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
You can use the same transcriptome index as bulk RNA-seq, however,
due to the shorter read lengths as opposed to bulk, you may want to build the 
index with a smaller `-k`.
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
Salmon for bulk RNA-seq. 
You downloaded two files:
- `R1` contain the barcodes for cells as well as the UMIs
- `R2` files contain the full reads for that sample.  

You'll recognize a lot of these options as the same as regular `Salmon` such as
- `-l` to designate library type
- `-1` and `-2` for file input
- `-o` to designate a folder for output

### `Alevin` options summary:

#### `-l`
As mentioned, `-l` is for designating library type. For all single-cell quant, 
you will want to use the `ISR` library type. 
See [Salmon's documentation](https://salmon.readthedocs.io/en/latest/library_type.html)
for more information on fragment library types.

#### `--chromium`
Because we are using 10X chromium data, we have to use this flag. However,
Drop-seq data is also supported, and in this case you would use a `--dropseq`
flag instead of this.

#### `--tgMap`
This is needed to supply a transcript to gene key that Alevin will use to
quantify the genes. For our example, we've premade the file `genes_2_tx.tsv` from
the Ensembl transcriptome that we indexed above. The file has to be a tsv file.

#### `--dumpCsvCounts`
Using this option will make the counts in a csv file and make it easier for us to
import into R later.

#### `--dumpFeatures`
This option will print out information that we will need for quality checks
later on, including files with information on the UMIs and cell barcodes.

#### Running Salmon Alevin
Copy and paste this in your command line to run Alevin quantification.
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

### Step 4: Perform QC checks with `alevinQC`
In order to perform quality control checks, we will need to open R.
Alevin provides count data output for each transcript and cell. To read this 
data into R, we will import a function from the script `ReadAlevin.R` which is 
located in the `scripts` folder.
```r 
# Import the function to read alevin output data
source(file.path("scripts", "ReadAlevin.R"))

# Read in the data
alevin_file <- ReadAlevin("alevin_output")
```
Now that our data is imported into the R environment, we can run quality control
checks using `alevinQC` package.
This will provide html output with graphs evaluating the data. 

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

TODO: Make a link to a "good" and "bad" alevinQC output "For future reference"
