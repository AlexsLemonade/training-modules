# Bulk RNA-seq Cheatsheet

#### The tables below consist of valuable functions or commands that will help you through this module.
##### Each table represents a different library/tool and its corresponding commands.
> Please note that these tables are not intended to tell you all the information you need to know about each command.
>
> The hyperlinks found in each piece of code will take you to the documentation for further information on the usage of each command.
Please be aware that the documentation will generally provide information about the given function's most current version (or a recent version, depending on how often the documentation site is updated).
This will usually (but not always!) match what you have installed on your machine.
If you have a different version of R or other R packages, the documentation may differ from what you have installed.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [Base `R`](#base-r)
- [`DESeq2`](#deseq2)
- [`FastQC` and `fastp`](#fastqc-and-fastp)
- [`ggplot2`](#ggplot2)
- [`tximeta` and `SummarizedExperiment`](#tximeta-and-summarizedexperiment)
- [`stringr`, `readr`, `dplyr`](#stringr-readr-dplyr)
- [`matrixStats`](#matrixstats)
- [`ComplexHeatmap` and `EnhancedVolcano`](#complexheatmap-and-enhancedvolcano)
- [Salmon](#salmon)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

<div style="page-break-after: always;"></div>

### Base `R`

Read the [Base `R` package documentation](https://rdrr.io/r/).

| Library/Package               | Piece of Code                                                 | What it's called      | What it does                                                             |
|-------------------------------|--------------------------------------------------------------|--------------------------------|--------------------------------------------------------------------------|
| Base `R`                | [`list.files()`](https://rdrr.io/r/base/list.files.html)             | List files                | Produces a character vector of files or directories in the specified directory |
| Base `R`                | [`names()`](https://rdrr.io/r/base/names.html)                       | Names                     | Gets or sets the names of an object                                      |
| Base `R`                | [`colnames()`](https://rdrr.io/r/base/colnames.html)           | Column names              | Gets or sets the column names of a matrix or data frame                  |
| Base `R`                | [`all.equal()`](https://rdrr.io/r/base/all.equal.html)               | All equal                 | Checks if two R objects are nearly equal                                 |
| Base `R`                | [`attr()`](https://rdrr.io/r/base/attr.html)                         | Object Attributes         | Gets or sets the attributes of an object                                 |
| Base `R`                | [`rowSums()`](https://rdrr.io/r/base/colSums.html)                   | Row Sums                  | Returns the sum of the rows in a numeric matrix-like object (i.e.. a matrix, data.frame, etc.) |
| Base `R`                | [`rowMeans()`](https://rdrr.io/r/base/colSums.html)                  | Row Means                 | Returns the mean of the rows in a number matrix-like object (i.e.. a matrix, data.frame, etc.) |
| Base `R`                | [`relevel()`](https://rdrr.io/r/stats/relevel.html)                  | Relevel                   | Reorders the levels of a factor as specified                             |
| Base `R`                | [`summary()`](https://rdrr.io/r/base/summary.html)                   | Object Summary            | Returns a result summary of an object                                    |
| Base `R`                | [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)       | Data Frame                | Checks if an object is a data.frame, and transforms the object into one, if possible |

<div style="page-break-after: always;"></div>

### `DESeq2`

Read the [`DESeq2` package documentation (PDF)](https://www.bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf), and the [package vignette by Love, Anders, and Huber](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).

| Library/Package               | Piece of Code                                                 | What it's called      | What it does                                                             |
|-------------------------------|--------------------------------------------------------------|--------------------------------|--------------------------------------------------------------------------|
| `DESeq2`                | [`vst()`](https://rdrr.io/bioc/DESeq2/man/vst.html)                   | Variance Stabilizing Transformation                  | Applies variance stabilizing transformation to data (log2-like scale)                                      |
| `DESeq2`                | [`DESeqDataSet()`](https://rdrr.io/bioc/DESeq2/man/DESeqDataSet.html)| DESeqDataSet constructor that can take a `SummarizedExperiment`| Creates a DESeqDataSet object                                            |
| `DESeq2`                | [`DESeqDataSetFromMatrix()`](https://rdrr.io/bioc/DESeq2/man/DESeqDataSet.html)| DESeqDataSet constructor        | Creates a DESeqDataSet object from a matrix of count data                                           |
| `DESeq2`                | [`DESeq()`](https://rdrr.io/bioc/DESeq2/man/DESeq.html)                 | Differential Expression Analysis Based on the Negative Binomial Distribution | Estimates size factors, estimates dispersion, and performs negative binomial fitting and Wald statistics as steps in the default DESeq2 differential expression analysis |
| `DESeq2`                | [`plotPCA()`](https://rdrr.io/bioc/DESeq2/man/plotPCA.html)             | PCA plot                                        | Produces a principal component analysis plot for transformed data. It can be used to visually inspect the data, which might allow an analyst to identify batch effects.   |
| `DESeq2`                | [`counts()`](https://rdrr.io/bioc/DESeq2/man/counts.html)              | Counts                                          | Returns count matrix from `DESeqDataSet` object                         |
| `DESeq2`                | [`results()`](https://rdrr.io/bioc/DESeq2/man/results.html)             | Results                                         | Returns the results table from a DESeq2 analysis                         |
| `DESeq2`                | [`assay()`](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#extracting-transformed-values)             | Assay                                           | Returns matrix from the `assay` slot of a `DESeqDataSet` object                                        |

<div style="page-break-after: always;"></div>

### `FastQC` and `fastp`

Read the [`FastQC` documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and the [`fastp` documentation](https://github.com/OpenGene/fastp).

| Library/Package               | Piece of Code                                                 | What it's called      | What it does                                                             |
|-------------------------------|--------------------------------------------------------------|--------------------------------|--------------------------------------------------------------------------|
| `fastp`                 | [`fastp`](https://github.com/OpenGene/fastp)                    | FASTQ preprocessor                              | Preprocesses FASTQ files through adapter trimming, quality filtering, length filtering, and a number of additional options  |
| `FastQC`                | [`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)                   | FASTQC (Quality Control)                        | Performs quality control checks on raw sequence data and outputs a QC (quality control) report   |

### `ggplot2`

Read the [`ggplot2` package documentation](https://ggplot2.tidyverse.org/), an [overall reference for `ggplot2` functions](https://ggplot2.tidyverse.org/reference/index.html), and a [vignette on the usage of the `ggplot2` aesthetics](https://ggplot2.tidyverse.org/articles/ggplot2-specs.html).
Additional vignettes are available from the "Articles" dropdown menu on this webpage.

| Library/Package               | Piece of Code                                                 | What it's called      | What it does                                                             |
|-------------------------------|--------------------------------------------------------------|--------------------------------|--------------------------------------------------------------------------|
| `ggplot2`               | [`ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)                 | GG Save                                         | Saves the last plot in working directory                                 |
| `ggplot2`               | [`last_plot()`](https://ggplot2.tidyverse.org/reference/last_plot.html)           | Last plot                                       | Returns the last plot produced                                           |
| `ggplot2`               | [`geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html)         | Geom point                                      | Creates a scatterplot (when added to the `ggplot()` function)            |
| `ggplot2`               | [`labs()`; `xlab()`; `ylab()`](https://ggplot2.tidyverse.org/reference/labs.html)           | X Axis Labels; Y Axis Labels                    | Modifies the plot labels, with specific functions for the x axis and y axis, respectively        |
| `ggplot2`               | [`coord_fixed()`](https://ggplot2.tidyverse.org/reference/coord_fixed.html)       | Cartesian Coordinates with Fixed Aspect Ratio   | Coerces the coordinates on the plot to represent a fixed specified ratio |

<div style="page-break-after: always;"></div>

### `tximeta` and `SummarizedExperiment`

Read the [`tximeta` package documentation (PDF)](https://bioconductor.org/packages/release/bioc/manuals/tximeta/man/tximeta.pdf), and the [package vignette by Love _et al._](https://bioconductor.org/packages/devel/bioc/vignettes/tximeta/inst/doc/tximeta.html).
Read the [`SummarizedExperiment` package documentation (PDF)](http://bioconductor.org/packages/release/bioc/manuals/SummarizedExperiment/man/SummarizedExperiment.pdf), and the [package vignette by Morgan _et al._](https://www.bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html).

| Library/Package               | Piece of Code                                                 | What it's called      | What it does                                                             |
|-------------------------------|--------------------------------------------------------------|--------------------------------|--------------------------------------------------------------------------|
| `tximeta`| [`tximeta()`](https://rdrr.io/bioc/tximeta/man/tximeta.html)           | tximeta           | Imports transcript-level estimates, attaches transcriptome annotation, and returns a `SummarizedExperiment` object  |
| `tximeta`| [`makeLinkedTxome()`](https://rdrr.io/bioc/tximeta/man/linkedTxome.html)| Make Linked Transcriptome | Sets up transcriptome annotation to be used by the `tximeta()` function (Only necessary if `tximeta()` fails to find annotation, like for non-human, non-mouse species data) |
| `tximeta`| [`summarizeToGene()`](https://rdrr.io/bioc/tximeta/man/summarizeToGene.html)| Summarize to Gene | Takes a `SummarizedExperiment` that was set up by `tximeta` and summarizes transcript data to the gene-level |
| `SummarizedExperiment`| [`rowData()`; `colData()`](https://rdrr.io/bioc/SummarizedExperiment/man/SummarizedExperiment-class.html)| Row/Col Data| Accesses the row or column data from a `SummarizedExperiment` object|
| `SummarizedExperiment`| [`assay()`; `assayNames()`](https://rdrr.io/bioc/SummarizedExperiment/man/SummarizedExperiment-class.html)| Assay or Assay Names| Accesses the assay data or the names of the assays from `SummarizedExperiment` object|


### `stringr`, `readr`, `dplyr`

Documentation for each of these packages can be accessed by clicking the package name in the table below.

| Library/Package               | Piece of Code                                                 | What it's called      | What it does                                                             |
|-------------------------------|--------------------------------------------------------------|--------------------------------|--------------------------------------------------------------------------|
| [`stringr`](https://stringr.tidyverse.org/index.html)            |[`word()`](https://stringr.tidyverse.org/reference/word.html)       | Word                            | Extracts words from a character vector                         |
| [`readr`](https://readr.tidyverse.org/index.html)                 |[`write_rds()`](https://readr.tidyverse.org/reference/read_rds.html) | Write RDS                      | Writes data to a .RDS output file                                 |
| [`dplyr`](https://dplyr.tidyverse.org/)                 | [`pull()`](https://dplyr.tidyverse.org/reference/pull.html)                      | Pull                              | Extracts a variable (column) as a vector                                          |


### `matrixStats`

Read the [`matrixStats` package documentation (PDF)](https://cran.r-project.org/web/packages/matrixStats/matrixStats.pdf) and [summary of functions](https://cran.r-project.org/web/packages/matrixStats/vignettes/matrixStats-methods.html).

| Library/Package               | Piece of Code                                                 | What it's called      | What it does                                                             |
|-------------------------------|--------------------------------------------------------------|--------------------------------|--------------------------------------------------------------------------|
| `matrixStats`                 | [`rowVars()`](https://www.rdocumentation.org/packages/matrixStats/versions/1.3.0/topics/rowVars) | Row variance                   | Estimates the variance for each row in a matrix |
| `matrixStats`                 | [`rowSds()`](https://search.r-project.org/CRAN/refmans/matrixStats/html/rowSds.html) | Row standard deviations | Estimates the standard deviations for each row in a matrix |

<div style="page-break-after: always;"></div>

### `ComplexHeatmap` and `EnhancedVolcano`

Read the [`ComplexHeatmap` Complete Reference](https://jokergoo.github.io/ComplexHeatmap-reference/book/).
Read the [`EnhancedVolcano` manual](https://bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html).

| Library/Package               | Piece of Code                                                 | What it's called      | What it does                                                             |
|-------------------------------|--------------------------------------------------------------|--------------------------------|--------------------------------------------------------------------------|
| `ComplexHeatmap`              | [`Heatmap()`](https://rdrr.io/bioc/ComplexHeatmap/man/Heatmap.html) | Heatmap constructor | Constructs a `Heatmap` class object that can then be used to plot a heatmap    |
| `ComplexHeatmap`              | [`HeatmapAnnotation()`](https://rdrr.io/bioc/ComplexHeatmap/man/HeatmapAnnotation.html) | Heatmap annotation constructor | Constructs a `HeatmapAnnotation` class object that can be used to annotate a heatmap |
| `EnhancedVolcano`             | [`EnhancedVolcano()`](https://rdrr.io/bioc/EnhancedVolcano/man/EnhancedVolcano.html) | Enhanced volcano plot | Constructs a highly customizable volcano plot with colored and labeled points |


### Salmon

Read the [Salmon tool documentation](https://salmon.readthedocs.io/en/latest/).

| Tool               | Piece of Code                                                 | What it's called      | What it does                                                             |
|-------------------------------|--------------------------------------------------------------|--------------------------------|--------------------------------------------------------------------------|
| Salmon                | [`salmon index`](https://salmon.readthedocs.io/en/latest/salmon.html)             | Salmon index                                    | Builds a transcriptome index which is required for Salmon quantification (from the command line) |
| Salmon                | [`salmon quant`](https://salmon.readthedocs.io/en/latest/salmon.html)             | Salmon quantification                           | Runs Salmon’s quantification of transcript expression (from the command line)                    |

<div style="page-break-after: always;"></div>

