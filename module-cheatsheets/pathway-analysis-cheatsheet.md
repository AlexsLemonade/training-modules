# Pathway Analysis Cheatsheet

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

- [`msigdbr`](#msigdbr)
- [AnnotationDbi](#annotationdbi)
- [Base `R`](#base-r)
- [`enrichplot`](#enrichplot)
- [`clusterProfiler`](#clusterprofiler)
- [`GSVA`](#gsva)
- [`qusage`](#qusage)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

<div style="page-break-after: always;"></div>

### `msigdbr`

Read the [`msigdbr` documentation](https://rdrr.io/cran/msigdbr/f/README.md).

| Library/Package| Piece of Code| What it's called| What it does |
|----------------|--------------|-----------------|--------------|
| `msigdbr`| [`msigdbr()`](https://rdrr.io/cran/msigdbr/man/msigdbr.html)| Retrieve the MSigDB gene sets data frame | Retrieves the specified MSigDB dataset |

### AnnotationDbi

Read the [`AnnotationDbi` package vignette (PDF)](http://bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf).

| Library/Package                                  | Piece of Code                               | What it's called                               | What it does                                                                                                                                   |
|-----------------------------------------|---------------------------------------------|-------------------------|------------------------------------------------------------------------|
| `AnnotationDbi`        | [`keytypes()`](https://rdrr.io/bioc/AnnotationDbi/man/AnnotationDb-class.html)      | Keytypes     | Returns a character vector of column names/types of gene identifiers (e.g. `ENSEMBL` ) available in an `AnnotationDbi` package.                                              |
| `AnnotationDbi`       | [`mapIds()`](https://rdrr.io/bioc/AnnotationDbi/man/AnnotationDb-class.html) | Mapped IDs       | Extracts the mapped IDs for a set of gene identifiers. The types of gene identifiers (e.g. `ENSEMBL` or `ENTREZ`) are supplied to arguments: `keytype` (type of gene identifiers we are providing in the `keys` argument) and  `column` (type of gene identifiers we want returned).                                 |

### Base `R`

Read the [Base `R` documentation](https://rdrr.io/r/).

|Library/Package|Piece of Code|What it's called| What it does|
|---------------|-------------|----------------|-------------|
| Base `R`| [`fisher.test()`](https://rdrr.io/r/stats/fisher.test.html)| Fisher's Exact Test | Performs the Fisher's exact test for testing the null of independence of rows and columns for a given matrix or data.frame with count data |
| Base `R` | [`setdiff()`](https://rdrr.io/cran/probs/man/setdiff.html) | Set difference | Returns the difference of two sets (e.g., vectors) |

<div style="page-break-after: always;"></div>

<!--
### `DESeq2`

Read the [`DESeq2` package documentation (PDF)](https://www.bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf), and the [package vignette by Love, Anders, and Huber](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).

| Library/Package               | Piece of Code                                                 | What it's called      | What it does                                                             |
|-------------------------------|--------------------------------------------------------------|--------------------------------|--------------------------------------------------------------------------|
| `DESeq2`                | [`lfcShrink()`](https://rdrr.io/bioc/DESeq2/man/lfcShrink.html)            | Shrink Log Fold Changes                                        | Adds shrunken log2 fold changes to the results of a `DESeqDataSet` object                         |
| `DESeq2`                | [`results()`](https://rdrr.io/bioc/DESeq2/man/results.html)             | Results                                         | Returns the results table from a DESeq2 analysis                         |
| `DESeq2`                | [`resultsNames()`](https://rdrr.io/bioc/DESeq2/man/results.html)            | Results Names                                        | Returns the names of the estimated effects or coefficients of the `DESeq` model                         |

-->

### `enrichplot`

Read the [`enrichplot` package documentation (PDF)](https://bioconductor.org/packages/devel/bioc/manuals/enrichplot/man/enrichplot.pdf).

| Library/Package               | Piece of Code                                                 | What it's called      | What it does                                                             |
|-------------------------------|--------------------------------------------------------------|--------------------------------|--------------------------------------------------------------------------|
| `enrichplot`                | [`dotplot()`](https://yulab-smu.github.io/clusterProfiler-book/chapter12.html#dot-plot)             | Dot plot                                         | Produces a dot plot for given enrichment results                         |
| `enrichplot`                | [`upsetplot()`](https://yulab-smu.github.io/clusterProfiler-book/chapter12.html#upset-plot)             | Upset plot                                       | Produces an upset plot, which shows the overlapping genes between gene sets, for given enrichment results                      |
| `enrichplot`                | [`gseaplot()`](https://yulab-smu.github.io/clusterProfiler-book/chapter12.html#running-score-and-preranked-list-of-gsea-result)             | GSEA plot                                       | Produces a plot visualization displaying the distribution of gene set and enrichment score                      |

### `clusterProfiler`

Read the [`clusterProfiler` package documentation (PDF)](https://www.bioconductor.org/packages/release/bioc/manuals/clusterProfiler/man/clusterProfiler.pdf).

| Library/Package               | Piece of Code                                                 | What it's called      | What it does                                                             |
|-------------------------------|--------------------------------------------------------------|--------------------------------|--------------------------------------------------------------------------|
| `clusterProfiler`                | [`enricher()`](https://rdrr.io/bioc/clusterProfiler/man/enricher.html)             | Enricher                                         | Performs a universal over-representation analysis for a given list of genes and gene sets or pathways             |
| `clusterProfiler`                | [`GSEA()`](https://rdrr.io/bioc/clusterProfiler/man/GSEA.html)             | Gene Set Enrichment Analysis (GSEA)                                         | Performs a universal gene set enrichment analysis on given preranked (sorted) named vector of statistics, where the names in the vector are gene identifiers of gene sets                     |

<div style="page-break-after: always;"></div>

### `GSVA`

Read the [`GSVA` package documentation](https://www.bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.html).

| Library/Package               | Piece of Code                                                 | What it's called      | What it does                                                             |
|-------------------------------|--------------------------------------------------------------|--------------------------------|--------------------------------------------------------------------------|
| `GSVA`                | [`gsva()`](https://rdrr.io/github/rcastelo/GSVA/man/gsva.html)             | Gene Set Variation Analysis (GSVA)                                         | Estimates gene set variation analysis enrichment scores on given gene expression matrix                       |
| `GSVA`                | [`gsvaParam()`](https://rdrr.io/github/rcastelo/GSVA/man/gsvaParam-class.html)             | Gene Set Variation Analysis (GSVA)  Parameters                                       | Specify parameters to use with `gsva()`               |


### `qusage`

Read the [`qusage` package documentation](https://rdrr.io/bioc/qusage/man/qusage.html).

| Library/Package|Piece of Code| What it's called| What it does  |
|----------------|-------------|-----------------|---------------|
| `qusage` | [`read.gmt()`](https://rdrr.io/bioc/qusage/man/read.gmt.html) | Read in `.gmt` files | Reads in gene set information from `.gmt` files |
