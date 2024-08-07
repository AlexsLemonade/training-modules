# Additional bulk RNA-seq tools and approaches

Our modules are designed to take up approximately half a day.
As a result, we do not cover every tool or method that would complement what is included in our main bulk RNA-seq module.
For your reference, we've put together what we might teach if we had more time or if they required fewer computational resources.
This is intended to be an introduction to these approaches and a "jumping off point" for further reading or experimentation.
We present them in order of increasing difficulty and/or departure from what is presented in training.

####  Contents

* [MultiQC](#multiqc)
* [Decoy sequence-aware selective alignment with Salmon](#decoy-sequence-aware-selective-alignment-with-salmon)

## MultiQC

[MultiQC](https://multiqc.info/) is a tool that aggregates output from many tools (nearly 150 are currently supported; see the list [here](https://multiqc.info/modules/)) into a single HTML report ([Ewels et al. _Bioinformatics._ 2016.](http://dx.doi.org/10.1093/bioinformatics/btw354)).
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

There are datasets already uploaded to the RStudio Server that you may want use for practicing these approaches or other techniques we've discussed.
For running MultiQC, take a look in the `~/shared-data/training-modules/RNA-seq/QC/` directory to see QC results from all of the datasets we have processed for training.
Note that we have to point you to `~/shared-data/` above because we didn't link all of this to your personal `training-modules` directory.

Run the following steps in the `Terminal` tab of RStudio.

Set your working directory to the top-level of the RNA-seq folder:

```bash
cd ~/training-modules/RNA-seq
```

**Now you're ready to run MultiQC.**
Because we're only interested in a report for the gastric cancer samples, we'll specify what directories QC and Salmon output can be found in.
We are using data and QC that was prepared previously by your instructors.

Run the following command:

```bash
multiqc \
  ~/shared-data/training-modules/RNA-seq/QC/gastric-cancer \
  ~/shared-data/training-modules/RNA-seq/data/gastric-cancer/salmon_quant/ \
  --outdir QC/gastric-cancer
```

Once this completes, you should have a report at `QC/gastric-cancer/multiqc_report.html`.

## Decoy sequence-aware selective alignment with Salmon

Transcript abundance estimates from lightweight mapping approaches like Salmon can differ from abundance estimates derived from traditional alignment approaches in experimental data ([Srivastava _et al._ 2020](https://doi.org/10.1186/s13059-020-02151-8)).
This differs from the takeaway of most prior work comparing lightweight mapping and traditional alignment; this is very likely due to the typical focus on _simulated data_ rather than experimental data.

To provide better results than lightweight mapping alone, Salmon includes a "selective alignment" method that is less computationally costly than traditional alignment while still offering improvements over lightweight mapping.

To apply the selective alignment method, you will need to first have an index that includes not only the transcripts of interest, but also a set of other potentially mapped sequences.
Instructions on creating such an index can be found in [the Salmon documentation](https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode).
One option is to include "decoy sequences" in Salmon index from genomic loci that are similar to annotated transcripts, to avoid falsely mapping fragments that arise from these unannotated regions to the transcripts of interest.
Alternatively (recommended), you can use the [_full genome_ as a decoy](https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/).

### Why don't we use decoy-aware selective alignment in training?

We select the methods we cover in training in part because it is feasible to run them on a standard laptop used for scientific computing if necessary.
Some laptops (including ones your instructors own!) are not well-equipped to run `salmon quant` with indices that include the decoy sequence information; they are considerably larger than the indices we use in training, which are cDNA-only indices.

### How can I try it out if I'm interested?

We have not extensively explored the memory and runtime requirements for the selective alignment mode, but we have successfully run human samples using Salmon with an index that includes partial decoy sequences on a Linux Desktop with 64 GB of RAM.

We recommend checking out [the Salmon documentation on indices](https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode) for more information, and following the [selective alignment tutorial](https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/) from the Salmon authors.
