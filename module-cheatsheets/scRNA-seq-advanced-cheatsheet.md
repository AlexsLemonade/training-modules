# Advanced scRNA-seq Cheatsheet

#### The tables below consist of valuable functions or commands that will help you through this module.
#### Each table represents a different library/tool and its corresponding commands.

##### You may also be interested in the following additional cheatsheets:

- Download the PDF for the [Introduction to R and Tidyverse cheatsheet](https://github.com/AlexsLemonade/training-modules/raw/master/module-cheatsheets/intro-to-R-tidyverse-cheatsheet.pdf)
- Download the PDF for the [Introduction to Single-Cell RNA sequencing cheatsheet](https://github.com/AlexsLemonade/training-modules/raw/master/module-cheatsheets/scRNA-seq-cheatsheet.pdf)

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
- [`batchelor` and `harmony`](#batchelor-and-harmony)
- [`pheatmap` and `EnhancedVolcano`](#pheatmap-and-enhancedvolcano)
- [`tidyverse` functions](#tidyverse-functions)
  - [`purrr` functions](#purrr-functions)
  - [`ggplot2` functions](#ggplot2-functions)
  - [`dplyr`, `tidyr`,`stringr`, and `tibble` functions](#dplyr-tidyrstringr-and-tibble-functions)
- [`Seurat` and `SCE` object conversion](#seurat-and-sce-object-conversion)
  - [Converting from `Seurat` to `SCE`](#converting-from-seurat-to-sce)
  - [Converting from `SCE` to `Seurat`](#converting-from-sce-to-seurat)
<!-- END doctoc generated TOC please keep comment here to allow auto update -->


<div style="page-break-after: always;"></div>

### `scater`

Read the [`scater` package documentation](https://rdrr.io/bioc/scater/), and a [vignette on its usage](http://www.bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/overview.html).

<br>

| Library/Package      | Piece of Code      | What it's called    | What it does  |
|----------------------|----------------------------|--------------------------------------------|--------------------------------------------------------------|
| `scater`| [`plotReducedDim()`](https://rdrr.io/bioc/scater/man/plotReducedDim.html)| Plot reduced dimensions | Plot a given reduced dimension slot from a `SingleCellExperiment` object by its name |
| `scater`| [`plotUMAP()`](https://rdrr.io/bioc/scater/man/plot_reddim.html)| Plot UMAP | Plot the "UMAP"-named reduced dimension slot from a `SingleCellExperiment` object |
| `scater`| [`plotExpression()`](https://rdrr.io/bioc/scater/man/plotExpression.html)| Plot expression | Plot expression values for all cells in a `SingleCellExperiment` object, using the `logcounts` assay by default|


### `miQC`

Read the [`miQC` package documentation](https://rdrr.io/github/greenelab/miQC/), and a [vignette on its usage](https://www.bioconductor.org/packages/release/bioc/vignettes/miQC/inst/doc/miQC.html).


| Library/Package      | Piece of Code      | What it's called    | What it does  |
|----------------------|----------------------------|--------------------------------------------|--------------------------------------------------------------|
| `miQC`| [`mixtureModel()`](https://rdrr.io/github/greenelab/miQC/man/mixtureModel.html)| Mixture model | Fit a `miQC` mixture model to a `SingleCellExperiment` object for use in filtering |
| `miQC`| [`filterCells()`](https://rdrr.io/github/greenelab/miQC/man/filterCells.html)| Filter cells | Filter cells from a `SingleCellExperiment` object based on a `miQC` model, returning a filtered `SingleCellExperiment` object |
| `miQC`| [`plotMetrics()`](https://rdrr.io/github/greenelab/miQC/man/plotMetrics.html)| Plot metrics | Plot percent of mitochondrial reads against the number of unique genes found for each cell |
| `miQC`| [`plotModel()`](https://rdrr.io/github/greenelab/miQC/man/plotModel.html)| Plot model | `miQC::plotMetics()` with the `miQC` fitted model overlaid |
| `miQC`| [`plotFiltering()`](https://rdrr.io/github/greenelab/miQC/man/plotFiltering.html)| Plot filtering | Plot percent of mitochondrial reads against the number of unique genes found, coloring points based on whether they will be filtered out or not |


### `batchelor` and `harmony`

Read the [`batchelor` package documentation](https://rdrr.io/cran/batchelor/), and a [vignette on its usage](https://rdrr.io/bioc/batchelor/f/vignettes/correction.Rmd).

Read the [`harmony` package documentation](https://rdrr.io/github/immunogenomics/harmony/man/harmony.html), and a [vignette on its usage](https://cran.r-project.org/web/packages/harmony/vignettes/quickstart.html).

| Library/Package      | Piece of Code      | What it's called    | What it does  |
|----------------------|----------------------------|--------------------------------------------|--------------------------------------------------------------|
| `batchelor`| [`MultiBatchPCA()`](https://rdrr.io/bioc/batchelor/man/multiBatchPCA.html)| Multi-batch PCA | Perform PCA across multiple gene expression matrices, weighted by batch size |
| `batchelor`| [`fastMNN()`](https://rdrr.io/bioc/batchelor/man/fastMNN.html)| Fast mutual nearest neighbors correction | Perform integration on an SCE object with mutual nearest neighbors using the `fastMNN` algorithm, returning an SCE object with batch-corrected principal components |
| `harmony`| [`HarmonyMatrix()`](https://rdrr.io/github/immunogenomics/harmony/man/HarmonyMatrix.html)| Perform `harmony` integration on a matrix | Perform integration with `harmony` on either a matrix of principle components or gene expression, returning a matrix of batch-corrected principal components  |


### `pheatmap` and `EnhancedVolcano`
Read the [`pheatmap` package documentation](https://rdrr.io/cran/pheatmap/).

Read the [`EnhancedVolcano` package documentation](https://rdrr.io/bioc/EnhancedVolcano/), and [vignette on its usage](https://bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html).

| Library/Package      | Piece of Code      | What it's called    | What it does  |
|----------------------|----------------------------|--------------------------------------------|--------------------------------------------------------------|
| `pheatmap`| [`pheatmap()`](https://rdrr.io/cran/pheatmap/man/pheatmap.html)| Pretty heatmap | Plot a (pretty!) clustered heatmap |
| `EnhancedVolcano`| [`EnhancedVolcano()`](https://rdrr.io/bioc/EnhancedVolcano/man/EnhancedVolcano.html)| Enhanced volcano | Plot a volcano plot to visualize differential expression analysis results |



### `tidyverse` functions



#### `purrr` functions

Read the [`purrr` package documentation](https://purrr.tidyverse.org/) and a [vignette on its usage](https://purrr.tidyverse.org/articles/base.html), and download the [`purr` package cheatsheet](https://github.com/rstudio/cheatsheets/blob/main/purrr.pdf).

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


#### `ggplot2` functions

Read the [`ggplot2` package documentation](https://ggplot2.tidyverse.org/) and an [overall reference for `ggplot2` functions](https://ggplot2.tidyverse.org/reference/index.html), and download the [`ggplot2` package cheatsheet](https://github.com/rstudio/cheatsheets/blob/main/data-visualization.pdf).

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



#### `dplyr`, `tidyr`,`stringr`, and `tibble` functions

Read the full documentation and download cheatsheets (where available) for these `tidyverse` packages at the following links:

- [`dplyr` documentation](https://dplyr.tidyverse.org/) and [`dplyr` cheatsheet](https://github.com/rstudio/cheatsheets/blob/main/data-transformation.pdf)
- [`tidyr` documentation](https://tidyr.tidyverse.org/) and [`tidyr` cheatsheet](https://github.com/rstudio/cheatsheets/blob/main/tidyr.pdf)
- [`stringr` documentation](https://stringr.tidyverse.org/) and [`stringr` cheatsheet](https://github.com/rstudio/cheatsheets/blob/main/strings.pdf)
- [`tibble` documentation](https://tibble.tidyverse.org/)

<br>

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





### `Seurat` and `SCE` object conversion

The `Seurat` documentation provides a [vignette about converting objects](https://satijalab.org/seurat/articles/conversion_vignette.html) between `SCE` and `Seurat` formats.

In addition, we provide some code examples below for how you can accomplish these conversions.
For all code examples below, it is assumed that the `SingleCellExperiment` library has been loaded into your R environment:

```r
library(SingleCellExperiment)
```


#### Converting from `Seurat` to `SCE`

The following example code assumes you have a `Seurat` object called `seurat_obj`.

```r
# Convert Seurat object to SCE object
sce_object <- Seurat::as.SingleCellExperiment(seurat_obj)
```

By default, all assays present in the `Seurat` object will be ported into the new `SCE` object.
To only specify that certain assays are retained, you can optionally provide the argument `assays`, as in: `assays = c("assays", "to", "keep")`.


#### Converting from `SCE` to `Seurat`

This [documentation from the `ScPCA`](https://scpca.readthedocs.io/en/latest/faq.html#what-if-i-want-to-use-seurat-instead-of-bioconductor) introduces how to convert `SCE` objects to `Seurat` objects.
Although this documentation was written for `ScPCA` datasets, the steps generally apply to any `SCE` object.
Briefly, here is how you can convert a `Seurat` to `SCE` object, focusing on porting over _assays_.
We provide `"RNA"` as the argument to `assay`, as this `Seurat`'s default name for raw count matrices.

The following example code assumes you have an `SCE` object called `sce_object`.

```r
# Create seurat object from the existing `sce_object`'s counts matrix,
seurat_object <- Seueat::CreateSeuratObject(counts = counts(sce_object),
                                            assay = "RNA",
                                            project = "name of your project goes here")

#
```

For further conversion steps, please see the [associated `ScPCA` documentation](https://scpca.readthedocs.io/en/latest/faq.html#what-if-i-want-to-use-seurat-instead-of-bioconductor).

Alternatively, we offer a conversion function `sce_to_seurat()` as part of our [`scpcaTools()` package](https://github.com/AlexsLemonade/scpcaTools/), which holds utilities used in the `ScPCA` workflow.
Again, although this function was written to convert `SCE` objects from `ScPCA`, it should generally work for most `SCE` objects.
Importantly, it will only retain a single assay, the raw counts, in the new `SCE` object, and it will not retain reduced dimension representations (e.g., PCA or UMAP).
Therefore, this function may be moslty useful at the early stages of processing before you have normalized counts and and calculated reduced dimensions.

You can obtain this package using the `remotes` package, which may also need to be installed first:

```r
# Install `remotes`, as needed:
install.packages("remotes")

# Install the current version of `scpcaTools`
remotes::install_github("AlexsLemonade/scpcaTools")

# Now, you can use the function, specifying the argument `assay` for which
seurat_object <- scpcaTools::sce_to_seurat(sce_object)
```