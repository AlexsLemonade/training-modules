# Pathway Analysis Cheatsheet

#### The tables below consist of valuable functions or commands that will help you through this module.
##### Each table represents a different library/tool and its corresponding commands.
> Please note that these tables are not intended to tell you all the information you need to know about each command.
>
> The hyperlinks found in each piece of code will take you to the documentation for further information on the usage of each command.

<div style="page-break-after: always;"></div>

<br>
<br>

### `msigdbr`

Read the `msigdbr` documentation [**here**](https://rdrr.io/cran/msigdbr/).

| Library/Package| Piece of Code| What it's called| What it does |
|----------------|--------------|-----------------|--------------|
| `msigdbr`| [`msigdbr()`](https://rdrr.io/cran/msigdbr/man/msigdbr.html)| **there was a typo here ("geom vertical line" aint it), but how do we PRONOUNCE this?**| Retrieves the specified MSigDB dataset |

<div style="page-break-after: always;"></div>



### AnnotationDbi 
Read the `AnnotationDbi` package vignette [**here**](http://bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf).

| Library/Package                                  | Piece of Code                               | What it's called                               | What it does                                                                                                                                   |
|-----------------------------------------|---------------------------------------------|-------------------------|------------------------------------------------------------------------|
| `AnnotationDbi`        | [`keytypes()`](https://rdrr.io/bioc/AnnotationDbi/man/AnnotationDb-class.html)      | Keytypes     | Returns a character vector of column names/keytypes (e.g., type of gene identifiers) available in an `AnnotationDbi` package.                                              |
| `AnnotationDbi`       | [`mapIDs()`](https://rdrr.io/bioc/AnnotationDbi/man/AnnotationDb-class.html) | Mapped IDs       | Extracts the mapped ids for a set of keys (e.g., gene identifiers) of a specific keytype       

<div style="page-break-after: always;"></div>


### Base `R`
Read the Base `R` documentation [**here**](https://rdrr.io/r/).


|Library/Package|Piece of Code|What it's called| What it does|
|---------------|-------------|----------------|-------------|
| Base `R`| [`fisher.test()`](https://rdrr.io/r/stats/fisher.test.html)| Fisher's Exact Test | Performs the Fisher's exact test for testing the null of independence of rows and columns for a given matrix or data.frame with count data |

<div style="page-break-after: always;"></div>

### `DESeq2`

Read the `DESeq2` package documentation [**here**](https://bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf) and the package vignette by Love, Anders, and Huber [**here**](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).

| Library/Package               | Piece of Code                                                 | What it's called      | What it does                                                             |
|-------------------------------|--------------------------------------------------------------|--------------------------------|--------------------------------------------------------------------------|
| `DESeq2`                | [`results()`](https://rdrr.io/bioc/DESeq2/man/results.html)             | Results                                         | Returns the results table from a DESeq2 analysis                         |
| `DESeq2`                | [`lfcShrink()`](https://rdrr.io/bioc/DESeq2/man/lfcShrink.html)            | Shrink Log Fold Changes                                        | Adds shrunken log2 fold changes to the results of a `DESeqDataSet` object                         |
| `DESeq2`                | [`resultsNames()`](https://rdrr.io/bioc/DESeq2/man/results.html)            | Results Names                                        | Returns the names of the estimated effects or coefficients of the `DESeq` model                         |

### `enrichplot`

Read the `enrichplot` package documentation [**here**](https://rdrr.io/bioc/enrichplot/).

| Library/Package               | Piece of Code                                                 | What it's called      | What it does                                                             |
|-------------------------------|--------------------------------------------------------------|--------------------------------|--------------------------------------------------------------------------|
| `enrichplot`                | [`dotplot()`](https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html#dot-plot)             | Dot plot                                         | Produces a dot plot for given enrichment results                         |
| `enrichplot`                | [`upsetplot()`](https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html#upset-plot             | Upset plot                                       | Produces an upset plot, which shows the overlapping genes between gene sets, for given enrichment results                      |
| `enrichplot`                | [`gseaplot()`](https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html#running-score-and-preranked-list-of-gsea-result)             | GSEA plot                                       | Produces a plot visualization displaying the distribution of gene set and enrichment score                      |

### `clusterProfiler`

Read the `clusterProfiler` package documentation [**here**](https://rdrr.io/bioc/clusterProfiler/).

| Library/Package               | Piece of Code                                                 | What it's called      | What it does                                                             |
|-------------------------------|--------------------------------------------------------------|--------------------------------|--------------------------------------------------------------------------|
| `clusterProfiler`                | [`enricher()`](https://rdrr.io/bioc/clusterProfiler/man/enricher.html)             | Enricher                                         | Performs a universal over-representation analysis for a given list of genes and gene sets or pathways             |
| `clusterProfiler`                | [`GSEA()`](https://rdrr.io/bioc/clusterProfiler/man/GSEA.html)             | Gene Set Enrichment Analysis (GSEA)                                         | Performs a universal gene set enrichment analysis on given preranked (sorted) named vector of statistics, where the names in the vector are gene identifiers of gene sets                     |

<div style="page-break-after: always;"></div>

### `GSVA`

Read the `GSVA` package documentation [**here**](https://rdrr.io/bioc/GSVA/).

| Library/Package               | Piece of Code                                                 | What it's called      | What it does                                                             |
|-------------------------------|--------------------------------------------------------------|--------------------------------|--------------------------------------------------------------------------|
| `GSVA`                | [`gsva()`](https://rdrr.io/bioc/GSVA/man/gsva.html)             | Gene Set Variation Analysis (GSVA)                                         | Estimates gene set variation analysis enrichment scores on given gene expression matrix                       |


### `qusage`

Read the `qusage` package documentation [**here**](https://rdrr.io/bioc/qusage/).

| Library/Package|Piece of Code| What it's called| What it does  |
|----------------|-------------|-----------------|---------------|
| `qusage` | [`read.gmt()`](https://rdrr.io/bioc/qusage/man/read.gmt.html) | Read in `.gmt` files | Reads in gene set information from `.gmt` files |
