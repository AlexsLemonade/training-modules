# Advanced scRNA-seq Cheatsheet

### The tables below consist of valuable functions or commands that will help you through this module.
### Each table represents a different library/tool and its corresponding commands.

#### You may also be interested in the following additional cheatsheets:

- Download the PDF for the [Introduction to R and Tidyverse cheatsheet](https://github.com/AlexsLemonade/training-modules/raw/master/module-cheatsheets/intro-to-R-tidyverse-cheatsheet.pdf)
- Download the PDF for the [Introduction to Single-Cell RNA sequencing cheatsheet](https://github.com/AlexsLemonade/training-modules/raw/master/module-cheatsheets/scRNA-seq-cheatsheet.pdf)
- Download the PDF for the [`SCE` and `Seurat` object comparison cheatsheet](https://github.com/AlexsLemonade/training-modules/raw/master/module-cheatsheets/SCE-Seurat-cheatsheet.pdf)

> Please note that these tables are not intended to tell you all the information you need to know about each command.
>
> The hyperlinks found in each piece of code will take you to the documentation for further information on the usage of each command.
> Please be aware that the documentation will generally provide information about the given function's most current version (or a recent version, depending on how often the documentation site is updated).
This will usually (but not always!) match what you have installed on your machine.
If you have a different version of R or other R packages, the documentation may differ from what you have installed.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [`scater`](#scater)
- [`miQC`](#miqc)
- [`SingleR`](#singler)
- [`pheatmap` and `EnhancedVolcano`](#pheatmap-and-enhancedvolcano)
- [`DESeq2` and pseudo-bulking functions](#deseq2-and-pseudo-bulking-functions)
- [`tidyverse` functions](#tidyverse-functions)
  - [`purrr` functions](#purrr-functions)
  - [`ggplot2` functions](#ggplot2-functions)
  - [`dplyr`, `tidyr`,`stringr`, and `tibble` functions](#dplyr-tidyrstringr-and-tibble-functions)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


<div style="page-break-after: always;"></div>

## `scater`

Read the [`scater` package documentation](https://rdrr.io/bioc/scater/), and a [vignette on its usage](http://www.bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/overview.html).


| Library/Package      | Piece of Code      | What it's called    | What it does  |
|----------------------|--------------------|---------------------|---------------|
| `scater`| [`plotReducedDim()`](https://rdrr.io/bioc/scater/man/plotReducedDim.html)| Plot reduced dimensions | Plot a given reduced dimension slot from a `SingleCellExperiment` object by its name |
| `scater`| [`plotUMAP()`](https://rdrr.io/bioc/scater/man/plot_reddim.html)| Plot UMAP | Plot the "UMAP"-named reduced dimension slot from a `SingleCellExperiment` object |
| `scater`| [`plotExpression()`](https://rdrr.io/bioc/scater/man/plotExpression.html)| Plot expression | Plot expression values for all cells in a `SingleCellExperiment` object, using the `logcounts` assay by default|


## `miQC`

Read the [`miQC` package documentation](https://rdrr.io/github/greenelab/miQC/), and a [vignette on its usage](https://www.bioconductor.org/packages/release/bioc/vignettes/miQC/inst/doc/miQC.html).


| Library/Package      | Piece of Code      | What it's called    | What it does  |
|----------------------|--------------------|---------------------|---------------|
| `miQC`| [`mixtureModel()`](https://rdrr.io/github/greenelab/miQC/man/mixtureModel.html)| Mixture model | Fit a `miQC` mixture model to a `SingleCellExperiment` object for use in filtering |
| `miQC`| [`filterCells()`](https://rdrr.io/github/greenelab/miQC/man/filterCells.html)| Filter cells | Filter cells from a `SingleCellExperiment` object based on a `miQC` model, returning a filtered `SingleCellExperiment` object |
| `miQC`| [`plotMetrics()`](https://rdrr.io/github/greenelab/miQC/man/plotMetrics.html)| Plot metrics | Plot percent of mitochondrial reads against the number of unique genes found for each cell |
| `miQC`| [`plotModel()`](https://rdrr.io/github/greenelab/miQC/man/plotModel.html)| Plot model | `miQC::plotMetics()` with the `miQC` fitted model overlaid |
| `miQC`| [`plotFiltering()`](https://rdrr.io/github/greenelab/miQC/man/plotFiltering.html)| Plot filtering | Plot percent of mitochondrial reads against the number of unique genes found, coloring points based on whether they will be filtered out or not |

<div style="page-break-after: always;"></div>
## `batchelor` and `harmony`

Read the [`batchelor` package documentation](https://rdrr.io/cran/batchelor/), and a [vignette on its usage](https://rdrr.io/bioc/batchelor/f/vignettes/correction.Rmd).

Read the [`harmony` package documentation](https://rdrr.io/github/immunogenomics/harmony/man/harmony.html), and a [vignette on its usage](https://cran.r-project.org/web/packages/harmony/vignettes/quickstart.html).

| Library/Package      | Piece of Code      | What it's called    | What it does  |
|----------------------|--------------------|---------------------|---------------|
| `batchelor`| [`MultiBatchPCA()`](https://rdrr.io/bioc/batchelor/man/multiBatchPCA.html)| Multi-batch PCA | Perform PCA across multiple gene expression matrices, weighted by batch size |
| `batchelor`| [`fastMNN()`](https://rdrr.io/bioc/batchelor/man/fastMNN.html)| Fast mutual nearest neighbors correction | Perform integration on an SCE object with mutual nearest neighbors using the `fastMNN` algorithm, returning an SCE object with batch-corrected principal components |
| `harmony`| [`HarmonyMatrix()`](https://rdrr.io/github/immunogenomics/harmony/man/HarmonyMatrix.html)| Perform `harmony` integration on a matrix | Perform integration with `harmony` on either a matrix of principle components or gene expression, returning a matrix of batch-corrected principal components  |

## `SingleR`

Read the [`SingleR` package documentation](https://rdrr.io/bioc/SingleR/), and an [e-book on its usage](http://bioconductor.org/books/release/SingleRBook/).

| Library/Package      | Piece of Code      | What it's called    | What it does  |
|----------------------|--------------------|---------------------|---------------|
| `SingleR` | [`trainSingleR()`](https://rdrr.io/bioc/SingleR/man/trainSingleR.html) | Train the SingleR classifier | Build a `SingleR` classifier model object from an annotated reference dataset |
| `SingleR` | [`classifySingleR()`](https://rdrr.io/bioc/SingleR/man/classifySingleR.html) | Classify cells with SingleR | Use a `SingleR` model object to assign cell types to the cells in an `SCE` object |
| `SingleR` | [`SingleR()`](https://rdrr.io/bioc/SingleR/man/SingleR.html) | Annotate scRNA-seq data | Combines `trainSingleR()` and `classifySingleR()` to assign cell types to an `SCE` object from an annotated reference dataset |


## `pheatmap` and `EnhancedVolcano`
Read the [`pheatmap` package documentation](https://rdrr.io/cran/pheatmap/).

Read the [`EnhancedVolcano` package documentation](https://rdrr.io/bioc/EnhancedVolcano/), and [vignette on its usage](https://bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html).

| Library/Package      | Piece of Code      | What it's called    | What it does  |
|----------------------|--------------------|---------------------|---------------|
| `pheatmap`| [`pheatmap()`](https://rdrr.io/cran/pheatmap/man/pheatmap.html)| Pretty heatmap | Plot a (pretty!) clustered heatmap |
| `EnhancedVolcano`| [`EnhancedVolcano()`](https://rdrr.io/bioc/EnhancedVolcano/man/EnhancedVolcano.html)| Enhanced volcano | Plot a volcano plot to visualize differential expression analysis results |



## `DESeq2` and pseudo-bulking functions

Read the [`DESeq2` package documentation](https://rdrr.io/bioc/DESeq2/), and a [vignette on its usage](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).


| Library/Package      | Piece of Code      | What it's called    | What it does  |
|----------------------|--------------------|---------------------|---------------|
| `scuttle`| [`aggregateAcrossCells()`](https://rdrr.io/bioc/scuttle/man/aggregateAcrossCells.html)| Aggregate data across groups of cells | Sum counts for each combination of features across groups of cells, commonly used to _pseudo-bulk_ `SCE` counts |
| `DESeq2`| [`DESeqDataSet()`](https://rdrr.io/bioc/DESeq2/man/DESeqDataSet.html)| DESeq Dataset | Establish a `DESeq` object from a pseudo-bulked `SingleCellExperiment` object or a bulk `SummarizedExperiment` object |
| `DESeq2`| [`estimateSizeFactors()`](https://rdrr.io/bioc/DESeq2/man/estimateSizeFactors.html)| Estimate size factors | Estimate size factors which are used to normalize counts for differential expression analysis  |
| `DESeq2`| [`rlog()`](https://rdrr.io/bioc/DESeq2/man/rlog.html)| Apply a regularized log transformation | Log2-transform counts in a `DESeq`  object for differential expression analysis   |
| `DESeq2`| [`plotPCA()`](https://rdrr.io/bioc/DESeq2/man/plotPCA.html)| Sample PCA plot for transformed data | Plot sample PCA from a log-transformed `DESeq` object to check for batch effects  |
| `DESeq2`| [`DESeq()`](https://rdrr.io/bioc/DESeq2/man/DESeq.html)| Perform differential expression analysis | Perform differential expression: Estimate size factors, transform data,  estimate dispersions, and perform testing. |
| `DESeq2`| [`plotDispEsts()`](https://rdrr.io/bioc/DESeq2/man/plotDispEsts.html)|Plot dispersion estimates | Plot dispersion estimates from a fitted `DESeq` object to evaluate model fit |
| `DESeq2`| [`results()`](https://rdrr.io/bioc/DESeq2/man/results.html) |Extract results from a `DESeq` analysis  | Extract results from a fitted `DESeq` object into a data frame |
| `DESeq2`| [`resultsNames()`](https://rdrr.io/bioc/DESeq2/man/results.html) | Extract results names | Return coefficient names from a fitted `DESeq` object |
| `DESeq2`| [`lfcShrink()`](https://rdrr.io/bioc/DESeq2/man/lfcShrink.html) | Shrink log2 fold changes | Add shrunken log2-fold changes to a results table produced by `DESeq2::results()` |


## `tidyverse` functions

### `purrr` functions

Read the [`purrr` package documentation](https://purrr.tidyverse.org/) and a [vignette on its usage](https://purrr.tidyverse.org/articles/base.html), and download the [`purr` package cheatsheet](https://github.com/rstudio/cheatsheets/raw/main/purrr.pdf).

| Library/Package| Piece of Code| What it's called| What it does |
|----------------|--------------|-----------------|--------------|
| `purrr`| [`map()`](https://purrr.tidyverse.org/reference/map.html)| map | Apply a function across each element of list; return a list |
| `purrr`| [`imap()`](https://purrr.tidyverse.org/reference/imap.html)| imap |  Apply a function across each element of list and its index/names; return a list |
| `purrr`| [`map2()`](https://purrr.tidyverse.org/reference/map2.html)| map2 |  Apply a function across each element of two lists at a time; return a list |
| `purrr`| [`reduce()`](https://purrr.tidyverse.org/reference/reduce.html)| Reduce |  Reduce a list to a single value by applying a given function |

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


### `ggplot2` functions

Read the [`ggplot2` package documentation](https://ggplot2.tidyverse.org/) and an [overall reference for `ggplot2` functions](https://ggplot2.tidyverse.org/reference/index.html), and download the [`ggplot2` package cheatsheet](https://github.com/rstudio/cheatsheets/raw/main/data-visualization.pdf).

| Library/Package| Piece of Code| What it's called| What it does |
|----------------|--------------|-----------------|--------------|
| `ggplot2` | [`geom_bar()`](https://ggplot2.tidyverse.org/reference/geom_bar.html) | Barplot  | Creates a barplot of counts for a given categorical variable when added as a layer to a `ggplot()` object |
| `ggplot2` | [`scale_fill_brewer()`](https://ggplot2.tidyverse.org/reference/scale_brewer.html) | Add brewer fill scale  | Apply a Brewer "fill" color palette to a categorical variable in a `ggplot()` object |
| `ggplot2` | [`guides()`](https://ggplot2.tidyverse.org/reference/guides.html) | Guides | Function to customize legend ("guide") appearance |
| `ggplot2` | [`facet_grid()`](https://ggplot2.tidyverse.org/reference/facet_grid.html) | Facet grid | Plot individual panels using specified variables to subset the data across rows and/or columns of a grid |
| `ggplot2` | [`vars()`](https://ggplot2.tidyverse.org/reference/vars.html) | Vars | Helper function to specify variables to `facet_grid()` or `facet_wrap()`|
| `ggplot2` | [`theme_bw()`](https://ggplot2.tidyverse.org/reference/ggtheme.html) | Black and white theme | Display `ggplot` with gridlines but a white background |
| `ggplot2` | [`theme()`](https://ggplot2.tidyverse.org/reference/theme.html) | Theme | Customize elements of a `ggplot` plot theme |
| `ggplot2` | [`element_text()`](https://ggplot2.tidyverse.org/reference/element.html) | Element text | Customize textual elements of a `ggplot` theme |



### `dplyr`, `tidyr`,`stringr`, and `tibble` functions

Read the full documentation and download cheatsheets (where available) for these `tidyverse` packages at the following links:

- [`dplyr` documentation](https://dplyr.tidyverse.org/) and [`dplyr` cheatsheet](https://github.com/rstudio/cheatsheets/raw/main/data-transformation.pdf)
- [`tidyr` documentation](https://tidyr.tidyverse.org/) and [`tidyr` cheatsheet](https://github.com/rstudio/cheatsheets/raw/main/tidyr.pdf)
- [`stringr` documentation](https://stringr.tidyverse.org/) and [`stringr` cheatsheet](https://github.com/rstudio/cheatsheets/raw/main/strings.pdf)
- [`tibble` documentation](https://tibble.tidyverse.org/)

| Library/Package| Piece of Code| What it's called| What it does |
|----------------|--------------|-----------------|--------------|
| `dplyr`| [`pull()`](https://dplyr.tidyverse.org/reference/pull.html) | Pull | Extract a single column from a data frame into a stand-alone vector |
| `dplyr`| [`count()`](https://dplyr.tidyverse.org/reference/count.html) | Count | Count the number of observations in each group of a data frame |
| `dplyr`| [`left_join()`](https://dplyr.tidyverse.org/reference/left_join.html)| Left join | Joins two data frames together, retaining only rows present in the first ("left") argument to the function |
| `dplyr`| [`relocate()`](https://dplyr.tidyverse.org/reference/relocate.html) | Relocate | Change column order in a data frame by relocating one or more columns |
| `dplyr`| [`case_when()`](https://dplyr.tidyverse.org/reference/case_when.html) | Case when | Return a value based on a set of `TRUE`/`FALSE` comparisons; a vectorized `if-else` |
| `tidyr`| [`pivot_longer()`](https://tidyr.tidyverse.org/reference/pivot_longer.html) | Pivot longer | Convert a "wide" format data frame to a "long" format data frame |
| `tibble`| [`as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html) | As tibble | Convert an object to a tibble |
| `stringr`| [`str_detect()`](https://stringr.tidyverse.org/reference/str_detect.html) | String detect | Returns `TRUE`/`FALSE` if a string contains a given substring|
| `stringr`| [`str_starts()`](https://stringr.tidyverse.org/reference/str_starts.html) | String starts | Returns `TRUE`/`FALSE` if a string starts with a given substring|
