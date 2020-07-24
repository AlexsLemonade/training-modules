# Pathway Analysis Cheatsheet

#### The tables below consist of valuable functions or commands that will help you through this module.
##### Each table represents a different library/tool and its corresponding commands.
> Please note that these tables are not intended to tell you all the information you need to know about each command.
>
> The hyperlinks found in each piece of code will take you to the documentation for further information on the usage of each command.

<div style="page-break-after: always;"></div>

### `msigdbr`

Read the `msigdbr` documentation [**here**](https://www.rdocumentation.org/packages/msigdbr/versions/7.1.1).

| Library/Package| Piece of Code| What it's called| What it does |
|----------------|--------------|-----------------|--------------|
| `msigdbr`| [`msigdbr()`](https://www.rdocumentation.org/packages/msigdbr/versions/7.1.1/topics/msigdbr)| geom vertical line| Retrieves the specified MSigDB dataset |

### AnnotationDbi 
Read the `AnnotationDbi` package vignette [**here**](http://bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf).

| Library/Package                                  | Piece of Code                               | What it's called                               | What it does                                                                                                                                   |
|-----------------------------------------|---------------------------------------------|-------------------------|------------------------------------------------------------------------|
| `AnnotationDbi`        | [`keytypes()`](https://bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf)      | Keytypes     | Returns a character vector of column names/types of gene identifiers (e.g. `ENSEMBL` ) available in an `AnnotationDbi` package.                                              |
| `AnnotationDbi`       | [`mapIds()`](https://www.rdocumentation.org/packages/AnnotationDbi/versions/1.44.0/topics/AnnotationDb-objects) | Mapped IDs       | Extracts the mapped IDs for a set of gene identifiers. The types of gene identifiers (e.g. `ENSEMBL` or `ENTREZ`) are supplied to arguments: `keytype` (type of gene identifiers we are providing in the `keys` argument) and  `column` (type of gene identifiers we want returned).                                 |  

### Base `R`

Read the Base `R` documentation [**here**](https://www.rdocumentation.org/packages/base/versions/3.5.1).

|Library/Package|Piece of Code|What it's called| What it does|
|---------------|-------------|----------------|-------------|
| Base `R`| [`fisher.test()`](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/fisher.test)| Fisher's Exact Test | Performs the Fisher's exact test for testing the null of independence of rows and columns for a given matrix or data.frame with count data |

<div style="page-break-after: always;"></div>

### `DESeq2`

Read the `DESeq2` package documentation [**here**](https://bioc.ism.ac.jp/packages/3.8/bioc/manuals/DESeq2/man/DESeq2.pdf) and the package vignette by Love, Anders, and Huber [**here**](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).

| Library/Package               | Piece of Code                                                 | What it's called      | What it does                                                             |
|-------------------------------|--------------------------------------------------------------|--------------------------------|--------------------------------------------------------------------------|
| `DESeq2`                | [`results()`](https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/results)             | Results                                         | Returns the results table from a DESeq2 analysis                         |
| `DESeq2`                | [`lfcShrink()`](https://rdrr.io/bioc/DESeq2/man/lfcShrink.html)            | Shrink Log Fold Changes                                        | Adds shrunken log2 fold changes to the results of a `DESeqDataSet` object                         |
| `DESeq2`                | [`resultsNames()`](https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/results)            | Results Names                                        | Returns the names of the estimated effects or coefficients of the `DESeq` model                         |

### `enrichplot`

Read the `enrichplot` package documentation [**here**](https://bioconductor.org/packages/devel/bioc/manuals/enrichplot/man/enrichplot.pdf).

| Library/Package               | Piece of Code                                                 | What it's called      | What it does                                                             |
|-------------------------------|--------------------------------------------------------------|--------------------------------|--------------------------------------------------------------------------|
| `enrichplot`                | [`dotplot()`](https://yulab-smu.github.io/clusterProfiler-book/chapter12.html#dot-plot)             | Dot plot                                         | Produces a dot plot for given enrichment results                         |
| `enrichplot`                | [`upsetplot()`](https://yulab-smu.github.io/clusterProfiler-book/chapter12.html#upset-plot)             | Upset plot                                       | Produces an upset plot, which shows the overlapping genes between gene sets, for given enrichment results                      |
| `enrichplot`                | [`gseaplot()`](https://yulab-smu.github.io/clusterProfiler-book/chapter12.html#running-score-and-preranked-list-of-gsea-result)             | GSEA plot                                       | Produces a plot visualization displaying the distribution of gene set and enrichment score                      |

### `clusterProfiler`

Read the `clusterProfiler` package documentation [**here**](https://www.rdocumentation.org/packages/clusterProfiler/versions/3.0.4).

| Library/Package               | Piece of Code                                                 | What it's called      | What it does                                                             |
|-------------------------------|--------------------------------------------------------------|--------------------------------|--------------------------------------------------------------------------|
| `clusterProfiler`                | [`enricher()`](https://www.rdocumentation.org/packages/clusterProfiler/versions/3.0.4/topics/enricher)             | Enricher                                         | Performs a universal over-represention/enrichment analysis for a given genes of interest list                  |
| `clusterProfiler`                | [`GSEA()`](https://www.rdocumentation.org/packages/clusterProfiler/versions/3.0.4/topics/GSEA)             | Gene Set Enrichment Analysis (GSEA)                                         | Performs a universal gene set enrichment analysis on given preranked (sorted) named vector of statistics, where the names in the vector are gene identifiers of gene sets                     |

<div style="page-break-after: always;"></div>

### `GSVA`

Read the `GSVA` package documentation [**here**](https://www.rdocumentation.org/packages/GSVA/versions/1.20.0).

| Library/Package               | Piece of Code                                                 | What it's called      | What it does                                                             |
|-------------------------------|--------------------------------------------------------------|--------------------------------|--------------------------------------------------------------------------|
| `GSVA`                | [`gsva()`](https://www.rdocumentation.org/packages/GSVA/versions/1.20.0/topics/gsva)             | Gene Set Variation Analysis (GSVA)                                         | Estimates gene set variation analysis enrichment scores on given gene expression matrix                       |


### `qusage`

Read the `qusage` package documentation [**here**](https://www.rdocumentation.org/packages/qusage/versions/2.4.0).

| Library/Package|Piece of Code| What it's called| What it does  |
|----------------|-------------|-----------------|---------------|
| `qusage` | [`read.gmt()`](https://www.rdocumentation.org/packages/qusage/versions/2.4.0/topics/read.gmt) | Read in `.gmt` files | Reads in gene set information from `.gmt` files |
