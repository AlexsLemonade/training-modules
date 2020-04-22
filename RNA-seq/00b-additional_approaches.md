# Additional bulk RNA-seq tools and approaches

Our modules are designed to take up approximately half a day.
As a result, we do not cover every tool or method that would complement what is included in our main bulk RNA-seq module.
For your reference, we've put together what we might teach if we had more time or if they required fewer computational resources.
This is intended to be an introduction to these approaches and a "jumping off point" for further reading or experimentation.
We present them in order of increasing difficulty and/or departure from what is presented in training.
 
####  Contents

* [MultiQC](#multiqc)
* [tximeta](#tximeta)
* [Decoy sequence-aware selective alignment with Salmon](#decoy-sequence-aware-selective-alignment-with-salmon)

## MultiQC

[MultiQC](https://multiqc.info/) is a tool that aggregates output from many tools (80 are currently supported; see the list [here](https://multiqc.info/docs/#multiqc-modules)) into a single HTML report ([Ewels et al. _Bioinformatics._ 2016.](http://dx.doi.org/10.1093/bioinformatics/btw354)).
It can be used to combine information, such as FastQC reports, for multiple RNA-seq samples in a project.
This can be helpful to get an overall picture of the samples in your experiment.
For example, fastp reports statistics before and after it processes samples.
If you were to perform _quality trimming_ with fastp and look at the before and after information across samples, it may tell you that, after trimming, a large portion of your reads were too short and were filtered out. 
That's likely information you'd want to know when performing downstream analyses.
MultiQC supports the three tools we present in the bulk RNA-seq module: FastQC, fastp, and Salmon. 
In addition, MultiQC is not limited to RNA-seq data; the website has example reports for Hi-C data and whole genome sequencing.

MultiQC is installed on the RStudio Server we use for training.
If you want to run it on the gastric cancer samples, we provide instructions below. 
Note that you will not have FastQC or fastp output for most samples, so there will be some information missing from the report.

#### Running MultiQC on the gastric cancer samples

Run the following steps in the `Terminal` tab of RStudio.

Set your working directory to the top-level of the RNA-seq folder:

```bash
cd training-modules/RNA-seq
```

**Now you're ready to run MultiQC.**
Because we're only interested in a report for the gastric cancer samples, we'll specify what directories QC and Salmon output can be found in.

Run the following command:

```bash
multiqc \
  QC/gastric_cancer \
  data/quant/gastric_cancer \
  --outdir QC/gastric_cancer/
```

Once this completes, you should have a report at `QC/gastric_cancer/multiqc_report.html`.

## tximeta

We use Salmon and tximport to quantify genes in bulk RNA-seq data in this module.
The authors of these tools have developed [`tximeta`](https://bioconductor.org/packages/release/bioc/html/tximeta.html), a Bioconductor package that has functionality to support reproducible research by automatically adding metadata (e.g., Salmon version used) to a specialized object in R called a [`SummarizedExperiment`](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html) and by linking the transcriptomes used for quantification to public sources (e.g., Ensembl). 
You can import quantifications from Salmon with `tximeta` and follow a similar path to what we present including the summarization to the gene level and the creation of a DESeq2 dataset.

**Read more about `tximeta` in [the package vignette]((https://bioconductor.org/packages/release/bioc/vignettes/tximeta/inst/doc/tximeta.html)) and [the README on GitHub](https://github.com/mikelove/tximeta/blob/master/README.md).**

If you have time to use `tximeta` with your own data, we encourage you to try it out, follow along the package documentation, and let us know what you think!

`tximeta` is not installed on the RStudio Server we use for training, but it can be installed by running the following command in the R console:

```R
BiocManager::install("tximeta")
```

## Decoy sequence-aware selective alignment with Salmon

In mid-2019, the folks that develop Salmon posted a preprint called _Alignment and mapping methodology influence transcript abundance estimation_ ([Srivastava et al. _bioRxiv._ 2019.](https://doi.org/10.1101/657874)).
We're summarizing our understanding of the work here, but we encourage you to take a look at the preprint yourself.

Srivastava et al. demonstrate that transcript abundance estimates from lightweight mapping approaches like Salmon can differ from abundance estimates derived from traditional alignment approaches in experimental data.
This differs from the takeaway of most prior work comparing lightweight mapping and traditional alignment; this is very likely due to the typical focus on _simulated data_ rather than experimental data.

In the preprint, the authors introduced a new approach termed "selective alignment" that is less computationally costly than traditional alignment while still offering improvements over lightweight mapping. 

The current version of Salmon (as of writing this) `v1.2.0` allows users to input sequences from unannotated genomic loci that are similar to annotated transcripts, termed decoy sequences, to avoid falsely mapping fragments that arise from these unannotated regions to transcripts.
This is termed a `salmon_partial_sa_index` [here](https://github.com/COMBINE-lab/salmon/tree/91091fc3650a3220f657a9f31616916513f0ad02#pre-computed-decoy-transcriptomes). 
As of `v1.0.0`, you can use the _full genome_ as decoy ([ref](https://github.com/COMBINE-lab/salmon/tree/91091fc3650a3220f657a9f31616916513f0ad02#pre-computed-decoy-transcriptomes)).

### Why don't we use decoy-aware selective alignment in training?

We select the methods we cover in training in part because it is feasible to run them on a standard laptop used for scientific computing if necessary. 
Some laptops (including ones your instructors own!) are not well-equipped to run `salmon quant` with indices that include the decoy sequence information; they are considerably larger than the indices we use in training, which are cDNA-only indices.

### How can I try it out if I'm interested?

We have not extensively explored the memory and runtime requirements for selective alignment-mode, but we have successfully run human samples using Salmon `v0.14.0` with an index that includes partial decoy sequences on a Linux Desktop with 64 GB of RAM.

We recommend checking out [this section of the Salmon README](https://github.com/COMBINE-lab/salmon/tree/91091fc3650a3220f657a9f31616916513f0ad02#pre-computed-decoy-transcriptomes) and following the [selective alignment tutorial](https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/) from the Salmon authors.
