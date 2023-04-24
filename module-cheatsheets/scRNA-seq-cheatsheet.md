# scRNA-seq Cheatsheet

#### The tables below consist of valuable functions or commands that will help you through this module.
##### Each table represents a different library/tool and its corresponding commands.
> Please note that these tables are not intended to tell you all the information you need to know about each command.
>
> The hyperlinks found in each piece of code will take you to the documentation for further information on the usage of each command.
> Please be aware that the documentation will generally provide information about the given function's most current version (or a recent version, depending on how often the documentation site is updated).
This will usually (but not always!) match what you have installed on your machine.
If you have a different version of R or other R packages, the documentation may differ from what you have installed.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [Base `R`](#base-r)
- [Salmon and `alevinQC`](#salmon-and-alevinqc)
- [`SingleCellExperiment`, `txmimeta`, and `DropletUtils`](#singlecellexperiment-txmimeta-and-dropletutils)
- [`scran` and `scater`](#scran-and-scater)
- [`purrr`, `stringr`, and `tibble`](#purrr-stringr-and-tibble)
- [`SingleR`](#singler)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

<div style="page-break-after: always;"></div>

### Base `R`

Read the Base [`R` documentation](https://rdrr.io/r/).

|Library/Package|Piece of Code|What it's called| What it does|
|---------------|-------------|----------------|-------------|
| Base `R`| [`rowSums()`](https://rdrr.io/r/base/colSums.html)| Row sums|Calculates sums for each row|
| Base `R`| [`colSums()`](https://rdrr.io/r/base/colSums.html)| Column sums| Calculates sums for each column|
| Base `R`| [`t()`](https://rdrr.io/r/base/t.html)  |Transpose| Returns the transpose of a matrix or data frame|
| Base `R`| [`prcomp()`](https://rdrr.io/r/stats/prcomp.html)| Principal Components Analysis|Executes a principal components analysis on specified matrix or data frame|
| Base `R`| [`<-function(x) { <code> }`](https://adv-r.hadley.nz/functions.html) | Function | Creates a function that would take the defined parameters as input and execute the commands within the curly braces  |


### Salmon and `alevinQC`
Read the command-line tool [Salmon documentation](https://salmon.readthedocs.io/en/latest/salmon.html).

Read the R package [`alevinQC` documentation](https://rdrr.io/bioc/alevinQC/).

Software/package | Piece of Code| What it's called| What it does |
|--------------|-----------------|--------------|-----------------|
| Salmon |[`salmon alevin`](https://salmon.readthedocs.io/en/latest/alevin.html)    | Salmon Alevin     | Runs the Alevin quantification from the command line  |
|[`alevinQC`](http://www.bioconductor.org/packages/devel/bioc//vignettes/alevinQC/inst/doc/alevinqc.html) | [`alevinQCReport()`](http://www.bioconductor.org/packages/devel/bioc//vignettes/alevinQC/inst/doc/alevinqc.html#generate-qc-report) | Alevin QC Report | Produces a QC (quality check) report from the `salmon alevin` output |


<div style="page-break-after: always;"></div>

### `SingleCellExperiment`, `txmimeta`, and `DropletUtils`

Read the [`SingleCellExperiment` package documentation (and e-book)](https://bioconductor.org/books/release/OSCA/), and a [vignette on its usage](https://rdrr.io/bioc/SingleCellExperiment/f/vignettes/intro.Rmd).
Note that some of the `SingleCellExperiment` functions link to documentation from other packages like `SummarizedExperiment` or `ExperimentSubset`.
In fact, `SingleCellExperiment` objects are based around existing Bioconductor functions in those packages, so the function usage is equivalent!

Read the [`tximeta` package documentation](https://rdrr.io/bioc/tximeta/), and a [vignette on its usage](https://rdrr.io/bioc/tximeta/f/README.md).

Read the [`DropletUtils` package documentation](https://rdrr.io/github/MarioniLab/DropletUtils/).



| Library/Package      | Piece of Code      | What it's called    | What it does  |
|----------------------|----------------------------|--------------------------------------------|--------------------------------------------------------------|
| `SingleCellExperiment` | [`SingleCellExperiment()`](https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html) | Single Cell Experiment| Creates a `SingleCellExperiment` object  |
| `SingleCellExperiment`| [`colData()`](https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html#2_Creating_SingleCellExperiment_instances) | Column Data | Extracts and stores cell-level metadata that describes features of the `SingleCellExperiment` object|
| `SingleCellExperiment`| [`rowData()`](https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html#2_Creating_SingleCellExperiment_instances)  | Row Data   | Extracts and stores gene-level metadata that describes features of the `SingleCellExperiment` object|
| `SingleCellExperiment`| [`assay()`](https://rdrr.io/bioc/SummarizedExperiment/man/SummarizedExperiment-class.html)  | Assay | Extracts and stores a given assay from a `SingleCellExperiment` object|
| `SingleCellExperiment`| [`assayNames()`](https://rdrr.io/bioc/SummarizedExperiment/man/SummarizedExperiment-class.html)  | Assay names  | Returns a vector of the names of all assays in a `SingleCellExperiment` object|
| `SingleCellExperiment`| [`logcounts()`](https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html#4_Convenient_access_to_named_assays)| Log counts| Extracts and stores log-transformed single-cell experiment count data as an assay of the `SingleCellExperiment` object|
| `SingleCellExperiment`| [`counts()`](https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html#4_Convenient_access_to_named_assays)| Counts| Extracts and stores raw single-cell experiment count data as an assay of the `SingleCellExperiment` object|
| `SingleCellExperiment`| [`reducedDim()`](https://rdrr.io/bioc/ExperimentSubset/man/reducedDim.html)| Reduced dim| Extracts or stores a given reduced dimension from a `SingleCellExperiment` object|
| `SingleCellExperiment`| [`reducedDimNames()`](https://rdrr.io/bioc/ExperimentSubset/man/reducedDimNames.html)| Reduced dim names| Returns a vector of the names of all reduced dimensions in a `SingleCellExperiment` object|
| `S4Vectors` | [`DataFrame()`](https://rdrr.io/bioc/S4Vectors/man/DataFrame-class.html)| Data frame | Not to be confused with `data.frame()` from Base R. This is a slightly different data frame-like object needed for storing information in `SingleCellExperiment` object's `colData` slot.|
| `tximeta` | [`tximeta()`](https://rdrr.io/bioc/tximeta/man/tximeta.html) | Transcript Quantification Import with Automatic Metadata | Load a directory of results produced by Salmon/or alevin output, including the associated metadata |
| `DropletUtils` | [`read10xCounts()`](https://rdrr.io/github/MarioniLab/DropletUtils/man/read10xCounts.html) | Read 10x counts | Load data from a 10x Genomics experiment into R |
| `DropletUtils` | [`emptyDrops()`](https://rdrr.io/github/MarioniLab/DropletUtils/man/emptyDrops.html) | Empty drops | Use the overall gene expression patterns in the sample to identify empty droplets |
| `DropletUtils` | [`emptyDropsCellRanger()`](https://rdrr.io/github/MarioniLab/DropletUtils/man/emptyDropsCellRanger.html) | Empty drops Cell Ranger | Use an  approach analogous to Cell Ranger's algorithm to identify empty droplets |


<div style="page-break-after: always;"></div>


### `scran` and `scater`


Read the [`scran` package documentation](https://rdrr.io/bioc/scran/), and a [vignette on its usage](https://rdrr.io/bioc/scran/f/vignettes/scran.Rmd).

Read the [`scater` package documentation](https://rdrr.io/bioc/scater/), and a [vignette on its usage](http://www.bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/overview.html).

| Library/Package      | Piece of Code      | What it's called    | What it does  |
|----------------------|----------------------------|--------------------------------------------|--------------------------------------------------------------|
| `scran` | [`quickCluster()`](https://rdrr.io/bioc/scran/man/quickCluster.html) | Quick Clustering     | Groups similar cells into clusters which are stored in the `SingleCellExperiment` object and are used for the calculation of size factors by `scran::computeSumFactors`|
| `scran` | [`computeSumFactors()`](https://rdrr.io/bioc/scran/man/computeSumFactors.html)  | Compute Sum Factors| Returns a numeric vector of computed sum factors for each cell cluster stored in the `SingleCellExperiment` object. The cluster-based size factors are deconvolved into cell-based size factors that are stored in the `SingleCellExperiment` object and used by the `scran::normalize` function for the normalization of each cell's gene expression profile|
| `scran`| [`getTopHVGs()`](https://rdrr.io/bioc/scran/man/getTopHVGs.html)| Get top highly variable genes | Identify variable genes in a `SingleCellExperiment` object, based on variance |
| `scran`| [`modelGeneVar()`](https://rdrr.io/bioc/scran/man/modelGeneVar.html)| model per gene variance | Model the per gene variance of a `SingleCellExperiment` object |
| `scran`| [`findMarkers()`](https://rdrr.io/bioc/scran/man/findMarkers.html)| Find marker genes | Find candidate marker genes for clusters of cells |
| `scater`| [`logNormCounts()`](https://rdrr.io/bioc/scuttle/man/logNormCounts.html) | Normalize log counts| Returns the `SingleCellExperiment` object with normalized expression values for each cell, using the size factors stored in the object
| `scater`| [`addPerCellQC()`](https://rdrr.io/bioc/scuttle/man/addPerCellQC.html)| Add per cell quality control | For a `SingleCellExperiment` object, calculate and add quality control per cell and store in `colData`  |
| `scater`| [`addPerFeatureQC()`](https://rdrr.io/bioc/scuttle/man/addPerCellQC.html)| Add per feature quality control | For a `SingleCellExperiment` object, calculate and add quality control per feature (genes usually) and store in `rowData`|
| `scater`| [`calculatePCA()`](https://rdrr.io/bioc/scater/man/runPCA.html)| Calculate PCA | Calculates principal components analysis on a `SingleCellExperiment` object, returning a PCA matrix |
| `scater`| [`runPCA()`](https://rdrr.io/bioc/scater/man/runPCA.html)| Run PCA | Calculates principal components analysis on a `SingleCellExperiment` object, returning an SCE object with a PCA reduced dimension |
| `scater`| [`calculateUMAP()`](https://rdrr.io/bioc/scater/man/runUMAP.html)| Calculate UMAP| Calculates uniform manifold approximate projection on a `SingleCellExperiment` object, returning a UMAP matrix |
| `scater`| [`runUMAP()`](https://rdrr.io/bioc/scater/man/runUMAP.html)| Run UMAP | Calculates uniform manifold approximate projection on a `SingleCellExperiment` object, returning an SCE object with a UMAP reduced dimension |
| `scater`| [`calculateTSNE()`](https://rdrr.io/bioc/scater/man/runTSNE.html)| Calculate t-SNE | Calculates t-stochastic neighbor embedding on a `SingleCellExperiment` object, returning an SCE object with a TSNE reduced dimension |
| `scater`| [`runTSNE()`](https://rdrr.io/bioc/scater/man/runTSNE.html)| Calculate UMAP | Calculates t-stochastic neighbor embedding on a `SingleCellExperiment` object, returning a t-SNE matrix  |
| `scater`| [`plotReducedDim()`](https://rdrr.io/bioc/scater/man/plotReducedDim.html)| Plot reduced dimensions | Plot a given reduced dimension slot from a `SingleCellExperiment` object by its name |
| `scater`| [`plotPCA()`](https://rdrr.io/bioc/scater/man/plot_reddim.html)| Plot PCA | Plot the "PCA"-named reduced dimension slot from a `SingleCellExperiment` object |
| `scater`| [`plotUMAP()`](https://rdrr.io/bioc/scater/man/plot_reddim.html)| Plot UMAP | Plot the "UMAP"-named reduced dimension slot from a `SingleCellExperiment` object |


<div style="page-break-after: always;"></div>


### `purrr`, `stringr`, and `tibble`

Read the [`purrr` package documentation](https://purrr.tidyverse.org/).

Read the [`stringr` package documentation](https://stringr.tidyverse.org/).

Read the [`tibble` package documentation](https://tibble.tidyverse.org/).


| Library/Package| Piece of Code| What it's called| What it does |
|----------------|--------------|-----------------|--------------|
| `purrr`| [`map()`](https://purrr.tidyverse.org/reference/map.html)| map | Apply a function across each element of list; return a list |
| `purrr`| [`map_df()`](https://purrr.tidyverse.org/reference/map.html)| map df |  Apply a function across each element of list; return a data frame |
| `purrr`| [`imap()`](https://purrr.tidyverse.org/reference/imap.html)| imap |  Apply a function across each element of list and its index/names |
| `stringr`| [`str_remove()`](https://stringr.tidyverse.org/reference/str_remove.html)| String remove | Remove matched string patterns |
| `tibble` |[`as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html) | As tibble | Coerce `data.frame` or matrix to a tibble |


Note that `purrr::map()` functions can take advantage of R's new (as of version 4.1.0) [anonymous function syntax](https://rdrr.io/r/base/function.html):

```r
# One-line syntax:
\(x) # function code goes here #

# Multi-line syntax:
\(x) {
  # function code goes      #
  # inside the curly braces #
}

# Example: Use an anonymous function with `purrr::map()`
# to get the colData's rownames for each SCE in `list_of_sce_objects`
purrr::map(
  list_of_sce_objects,
  \(x) rownames(colData(x))
)
```

<div style="page-break-after: always;"></div>

### `SingleR`

Read the [`SingleR` package documentation](https://rdrr.io/bioc/SingleR/), and an [e-book on its usage](http://bioconductor.org/books/release/SingleRBook/).

| Library/Package      | Piece of Code      | What it's called    | What it does  |
|----------------------|--------------------|---------------------|---------------|
| `SingleR` | [`trainSingleR()`](https://rdrr.io/bioc/SingleR/man/trainSingleR.html) | Train the SingleR classifier | Build a `SingleR` classifier model object from an annotated reference dataset |
| `SingleR` | [`classifySingleR()`](https://rdrr.io/bioc/SingleR/man/classifySingleR.html) | Classify cells with SingleR | Use a `SingleR` model object to assign cell types to the cells in an `SCE` object |
| `SingleR` | [`SingleR()`](https://rdrr.io/bioc/SingleR/man/SingleR.html) | Annotate scRNA-seq data | Combines `trainSingleR()` and `classifySingleR()` to assign cell types to an `SCE` object from an annotated reference dataset |
