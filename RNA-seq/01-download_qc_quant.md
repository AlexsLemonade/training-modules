## Introduction to RNA-seq data processing

We will first learn how to process RNA-seq data at the command line using two samples that were assayed with paired-end sequencing. 

These samples come from a project ([`PRJNA178120`](https://www.ebi.ac.uk/ena/data/view/PRJNA178120)) that includes 8 samples from normal gastric tissue, gastric cancer cell lines and primary gastric tumor cell cultures.

We'll download data from [European Nucleotide Archive](https://www.ebi.ac.uk/ena), do some quality control checks, and estimate the transcript abundances for these two samples.

**Our objectives:**

* Get comfortable with the command line
* Learn about project organization
* Become familiar with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [Salmon](https://combine-lab.github.io/salmon/)

Later, we will use the full dataset (n = 8) to explore how to summarize estimates to the gene level and do some exploratory data analyses with data the course directors have processed ahead of time.

---

For these exercises, we'll want to set our **LOCAL FOLDER** in Kitematic to the `RNA-seq` folder. 

Copy and paste the text in the code blocks below into your `Terminal` window in RStudio. 
It should be in the lower left hand corner as a tab next to `Console`.

We'll first want to set our working directory to the top-level of the RNA-seq folder like so

```bash
cd kitematic
```

### Download data

First, let's create a directory for our raw sequencing data (`.fastq.gz`) files.

```bash
mkdir data/fastq/gastric_cancer -p
```

We'll use the `-P` option for `wget` to download the files into our new directory.

#### SRR585570

```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR585/SRR585570/SRR585570_1.fastq.gz -P data/fastq/gastric_cancer 
```
```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR585/SRR585570/SRR585570_2.fastq.gz -P data/fastq/gastric_cancer 
```

#### SRR585574

```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR585/SRR585574/SRR585574_1.fastq.gz -P data/fastq/gastric_cancer 
```
```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR585/SRR585574/SRR585574_2.fastq.gz -P data/fastq/gastric_cancer 
```

### Quality control with FastQC

We'll use [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for quality control in command line mode.
Here's a link to the FastQC documentation: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/

**Let's take a look at some example reports from the authors of FastQC:**

* [Good Illumina data example report](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html) from FastQC
* [Bad Illumina data example report](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html) from FastQC

FastQC runs a series of quality checks on sequencing data and provides an HTML report. As the authors point out in the [docs](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/2%20Basic%20Operations/2.2%20Evaluating%20Results.html):

> It is important to stress that although the analysis results appear to give a pass/fail result, these evaluations must be taken in the context of what you expect from your library.

The [documentation for individual modules/analyses](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/) in FastQC is a great resource!

To save time, we'll only run one `fastq` file for each sample.
Let's make a directory to hold our QC results.

```bash
mkdir QC/gastric_cancer/ -p
```

We'll use the `-o` option for `fastqc` to place the reports generated in the directory we just created.

#### SRR585570

```bash
fastqc data/fastq/gastric_cancer/SRR585570_1.fastq.gz -o QC/gastric_cancer/
```

#### SRR585574

```bash
fastqc data/fastq/gastric_cancer/SRR585574_1.fastq.gz -o QC/gastric_cancer/
```

Let's look at the results together.

### Quantification with Salmon

We'll use [Salmon](https://combine-lab.github.io/salmon/) for quantifying transcript expression ([documentation](http://salmon.readthedocs.io/en/latest/)). 
Salmon ([Patro, et al. _Nature Methods._ 2017.](https://doi.org/10.1038/nmeth.4197)) is fast and requires very little memory, which makes it a great choice for running on your laptop and for a project like refine.bio!
We can use the output for downstream analyses like differential expression analysis and clustering. 
Below, we'll walk through the arguments options we'll use when running `salmon quant`. 
We'll be using Salmon in quasi-mapping mode.

#### Transcriptome index (`-i`)

Salmon requires a set of transcripts (what we want to quantify) in the form of a transcriptome index built with `salmon index`.
Building an index can take a while (but you only have to do it once!), so we've built the one we will use today ahead of time. 
Before we use it, we'll take a moment to give a bit of background.

One of the things we learned by using FastQC is read length. 
`salmon index` has a parameter `-k` which sets the k-mer length.
The most important point is that the recommended value for _k_ depends on the read length, with _k_ = 31 is appropriate for 75bp or longer reads.
The index we'll use now was built with `-k 31` and can be found here:

```
index/Homo_sapiens/long_index
```

#### `--gcBias`

With this option enabled, Salmon will attempt to correct for fragment GC-bias. 
Regions with high or low GC content tend to be underrepresented in sequencing data.

It should be noted that this is only appropriate for use with paired-end reads, as fragment length can not be inferred from single-end reads (see [this Github issue](https://github.com/COMBINE-lab/salmon/issues/83)).

#### `--seqBias`

With this option enabled, Salmon will attempt to correct for the bias that occurs when using random hexamer priming (preferential sequencing of reads when certain motifs appear at the beginning).

#### `--biasSpeedSamp` 

We'll set this value to `5` (the default) to speed up the GC bias correction. 
This should have very little effect on the results according to the Salmon documentation.

#### `-o`

Output directory, `salmon quant` should create this for us if it doesn't exist yet.

#### `-l`

We'll use `-l A` to allow Salmon to automatically infer the library type based on a subset of reads, but you can also provide the [library type](http://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype) to Salmon with this argument. 

### Running `salmon quant`

Because these are paired-end we will use `-1` and `-2`

#### SRR585570

```bash
salmon quant -i index/Homo_sapiens/long_index -l A \
        -1 data/fastq/gastric_cancer/SRR585570_1.fastq.gz \
        -2 data/fastq/gastric_cancer/SRR585570_2.fastq.gz \
        -o data/quant/gastric_cancer/SRR585570 \
        --gcBias --seqBias --biasSpeedSamp 5
```

#### SRR585574

```bash
salmon quant -i index/Homo_sapiens/long_index -l A \
        -1 data/fastq/gastric_cancer/SRR585574_1.fastq.gz \
        -2 data/fastq/gastric_cancer/SRR585574_2.fastq.gz \
        -o data/quant/gastric_cancer/SRR585574 \
        --gcBias --seqBias --biasSpeedSamp 5
```

Navigate to `data` > `quant` > `gastric_cancer` > `SRR585570` > `aux_info` and open `meta_info.json`.
Look for a field called `percent_mapped` -- what value does this sample have?

---

### Note on mapping validation

Newer versions of Salmon, including the one that we are using (`0.12.0`), can be run in mapping validation mode that is used when we provide the [`--validateMappings`](https://salmon.readthedocs.io/en/latest/salmon.html#validatemappings) flag to `salmon quant`.
This improves quantification estimates.
When using this option, it is recommended that you trim low quality bases and adapter content prior to quantification with Salmon.
The [authors of Salmon recommend](https://github.com/COMBINE-lab/salmon/releases/tag/v0.13.1) all users move to mapping validation as it is likely to be the default behavior in future releases.
In near- future versions of the workshop, we will very likely perform trimming and use `--validateMappings`. Stay tuned!