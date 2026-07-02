# Spatial Transcriptomics Cheatsheet

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

- [`SpatialExperiment`](#spatialexperiment)
- [`VisiumIO`](#visiumio)
- [`ggspavis`](#ggspavis)
- [`SpotSweeper`](#spotsweeper)
- [`Banksy`](#banksy)
- [`spacexr`](#spacexr)
- [`scuttle`, `scran` , and `scater`](#scuttle-scran--and-scater)
- [`bluster`](#bluster)
- [`patchwork`](#patchwork)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

<div style="page-break-after: always;"></div>


### `SpatialExperiment`


Read the [`SpatialExperiment` package documentation (and e-book)](https://bioconductor.org/books/release/OSTA/) and a [vignette on its usage](https://www.bioconductor.org/packages/release/bioc/vignettes/SpatialExperiment/inst/doc/SpatialExperiment.html).

`SpatialExperiment` objects are based on [`SingleCellExperiment` objects](https://bioconductor.org/books/release/OSCA/).
There are many functions that are shared between the two, with equivalent usage.
See [this `SingleCellExperiment` vignette](https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html) for example usage of shared functions.
You can also refer to the Data Lab cheatsheets for our "Introduction to scRNA-seq" (`scRNA-seq-cheatsheet.pdf`) and "Advanced scRNA-seq" (`scRNA-seq-advanced-cheatsheet.pdf`) training modules.



| Library/Package      | Piece of Code      | What it's called    | What it does  |
|----------------------|----------------------------|--------------------------------------------|--------------------------------------------------------------|
| `SpatialExperiment` | [`SpatialExperiment()`](https://www.bioconductor.org/packages/release/bioc/vignettes/SpatialExperiment/inst/doc/SpatialExperiment.html) | Spatial Experiment| Creates a `SpatialExperiment` object  |
| `SpatialExperiment` | [`spatialCoords()`](https://www.bioconductor.org/packages/release/bioc/vignettes/SpatialExperiment/inst/doc/SpatialExperiment.html#13_spatialCoords) | Spatial Experiment coordinates | Extracts spatial coordinates of spots in the `SpatialExperiment` object  |
| `SpatialExperiment` | [`imgData()`](https://www.bioconductor.org/packages/release/bioc/vignettes/SpatialExperiment/inst/doc/SpatialExperiment.html#14_imgData) | Spatial Experiment images | View the image data associated with a `SpatialExperiment` object   |




### `VisiumIO`

<!-- TODO: This package has not yet made it to rdrr -->
Read the R package [`VisiumIO` documentation](https://www.bioconductor.org/packages/release/bioc/html/VisiumIO.html) and a [vignette on its usage](https://www.bioconductor.org/packages/release/bioc/vignettes/VisiumIO/inst/doc/VisiumIO.html).

|Library/Package|Piece of Code|What it's called| What it does|
|---------------|-------------|----------------|-------------|
| `VisiumIO` | [`TENxVisium()`](https://www.bioconductor.org/packages/release/bioc/vignettes/VisiumIO/inst/doc/VisiumIO.html#tenxvisium) | 10x Visium | Constructor function for importing 10x Visium data |
| `VisiumIO`| [`import()`](https://www.bioconductor.org/packages/release/bioc/vignettes/VisiumIO/inst/doc/VisiumIO.html#importing-into-spatialexperiment) | Import | Loads 10x Visium data defined with the `TENxVisium()` constructor function |




### `ggspavis`

<!-- TODO: This package has not yet made it to rdrr -->
Read the [`ggspavis` package documentation](https://www.bioconductor.org/packages/release/bioc/html/ggspavis.html) and a [vignette on its usage](https://www.bioconductor.org/packages/release/bioc/vignettes/ggspavis/inst/doc/ggspavis_overview.html).



| Library/Package      | Piece of Code      | What it's called    | What it does  |
|----------------------|----------------------------|--------------------------------------------|--------------------------------------------------------------|
| `ggspavis` | [`plotCoords()`](https://www.bioconductor.org/packages/release/bioc/vignettes/ggspavis/inst/doc/ggspavis_overview.html#spot-shape---slide-seq-and-visium) | Plot coordinates | Creates a spot plot showing spatial locations in x/y with options for spot annotation  |
| `ggspavis` | [`plotVisium()`](https://www.bioconductor.org/packages/release/bioc/vignettes/ggspavis/inst/doc/ggspavis_overview.html#spot-shape---slide-seq-and-visium) | Plot Visium | Creates a spot plot for Visium data specifically showing spatial locations in x/y with options for spot annotation as well as an option to show the spatial image, e.g., the H&E image  |


<div style="page-break-after: always;"></div>



### `SpotSweeper`


<!-- TODO: This package has not yet made it to rdrr -->
Read the [`SpotSweeper` package documentation](https://www.bioconductor.org/packages/release/bioc/html/SpotSweeper.html) and a [vignette on its usage](https://www.bioconductor.org/packages/release/bioc/vignettes/SpotSweeper/inst/doc/getting_started.html).


| Library/Package      | Piece of Code      | What it's called    | What it does  |
|----------------------|----------------------------|--------------------------------------------|--------------------------------------------------------------|
| `SpotSweeper` | [`localOutliers()`](https://www.bioconductor.org/packages/release/bioc/vignettes/SpotSweeper/inst/doc/getting_started.html#identifying-local-outliers-using-spotsweeper)  | Local outliers | Detect local outlier spots in a `SpatialExperiment` based on a given quality control metric |
| `SpotSweeper`| [`plotQCmetrics()`](https://www.bioconductor.org/packages/release/bioc/vignettes/SpotSweeper/inst/doc/getting_started.html#visualizing-local-outliers)| Plot QC metrics | Create a spot plot highlighting outliers detected with `SpotSweeper::localOutliers()` |




### `Banksy`

<!-- TODO: This package has not yet made it to rdrr -->
Read the [`Banksy` package documentation](https://www.bioconductor.org/packages/release/bioc/html/Banksy.html), which contains several vignettes on its usage, including a [vignette on parameter selection](https://www.bioconductor.org/packages/release/bioc/vignettes/Banksy/inst/doc/parameter-selection.html).
Note that many of their vignettes use of image-based spatial transcriptomics data that is more fine-grained compared to 10x Visium.


| Library/Package      | Piece of Code      | What it's called    | What it does  |
|----------------------|----------------------------|--------------------------------------------|--------------------------------------------------------------|
| `Banksy` | [`computeBanksy()`](https://www.bioconductor.org/packages/release/bioc/vignettes/Banksy/inst/doc/domain-segment.html#running-banksy)  | Compute `Banksy` | Compute the component neighborhood matrices for the `Banksy` matrix |
| `Banksy`| [`runBanksyPCA()`](https://www.bioconductor.org/packages/release/bioc/vignettes/Banksy/inst/doc/domain-segment.html#running-banksy)| Run `Banksy` PCA | Run PCA on a `Banksy` matrix using a provided value for the spatial weighting parameter `lambda` |


### `spacexr`

<!-- TODO: This package has not yet made it to rdrr -->
Read the [`spacexr` package documentation](https://www.bioconductor.org/packages/release/bioc/html/spacexr.html), and a [vignette on using their RCTD algorithm](https://www.bioconductor.org/packages/release/bioc/vignettes/spacexr/inst/doc/rctd-tutorial.html).


| Library/Package      | Piece of Code      | What it's called    | What it does  |
|----------------------|----------------------------|--------------------------------------------|--------------------------------------------------------------|
| `spacexr` | [`createRctd()`](https://www.bioconductor.org/packages/release/bioc/vignettes/spacexr/inst/doc/rctd-tutorial.html#step-1-preprocess-data)  | Create an RCTD object | Preprocess data before performing deconvolution with RCTD |
| `spacexr`| [`runRCTD()`](https://www.bioconductor.org/packages/release/bioc/vignettes/spacexr/inst/doc/rctd-tutorial.html#step-2-run-rctd)| Run RCTD | Run the RCTD algorithm to decompose cell type mixtures |
| `spacexr`| [`plotCellTypeWeight()`](https://www.bioconductor.org/packages/release/bioc/vignettes/spacexr/inst/doc/rctd-tutorial.html#single-cell-type)| Plot a single cell type's weights | Plot pixel proportions for a specific cell type on a spatial coordinates plot |

<div style="page-break-after: always;"></div>


### `scuttle`, `scran` , and `scater`


Read the [`scuttle` package documentation](https://rdrr.io/bioc/scuttle/), and a [vignette on its usage](https://rdrr.io/bioc/scuttle/f/vignettes/overview.Rmd).

Read the [`scran` package documentation](https://rdrr.io/bioc/scran/), and a [vignette on its usage](https://rdrr.io/bioc/scran/f/vignettes/scran.Rmd).

Read the [`scater` package documentation](https://rdrr.io/bioc/scater/), and a [vignette on its usage](http://www.bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/overview.html).


| Library/Package      | Piece of Code      | What it's called    | What it does  |
|----------------------|----------------------------|--------------------------------------------|--------------------------------------------------------------|
| `scuttle`| [`addPerCellQC()`](https://rdrr.io/bioc/scuttle/man/addPerCellQC.html)| Add per cell quality control | For a `SingleCellExperiment` object, calculate and add quality control per cell and store in `colData`  |
| `scuttle` | [`computeLibraryFactors()`](https://rdrr.io/bioc/scuttle/man/librarySizeFactors.html)  | Compute Library Factors | Returns a numeric vector of computed size factors for each spot (or cell) stored in a `SpatialExperiment` (or `SingleCellExperiment`) object. The size factor is computed as the library size of each spot/cell after scaling them to have a mean of 1 across all spots/cells |
| `scuttle`| [`logNormCounts()`](https://rdrr.io/bioc/scuttle/man/logNormCounts.html)| Normalize log counts | Returns the `SpatialExperiment` (or `SingleCellExperiment`) object with normalized expression values for each spot (cell), using the size factors stored in the object |
| `scran`| [`getTopHVGs()`](https://rdrr.io/bioc/scran/man/getTopHVGs.html)| Get top highly variable genes | Identify variable genes in a `SingleCellExperiment` object, based on variance |
| `scran`| [`modelGeneVar()`](https://rdrr.io/bioc/scran/man/modelGeneVar.html)| model per gene variance | Model the per gene variance of a `SingleCellExperiment` object |
| `scran`| [`clusterCells()`](https://rdrr.io/github/MarioniLab/scran/man/clusterCells.html)| Cluster cells | Perform clustering on an SCE object using the `bluster` package |
| `scater`| [`runPCA()`](https://rdrr.io/bioc/scater/man/runPCA.html)| Run PCA | Calculates principal components analysis on a `SingleCellExperiment` object, returning an SCE object with a PCA reduced dimension |
| `scater`| [`runUMAP()`](https://rdrr.io/bioc/scater/man/runUMAP.html)| Run UMAP | Calculates uniform manifold approximate projection on a `SingleCellExperiment` object, returning an SCE object with a UMAP reduced dimension |
| `scater`| [`plotUMAP()`](https://rdrr.io/bioc/scater/man/plot_reddim.html)| Plot UMAP | Plot the "UMAP"-named reduced dimension slot from a `SingleCellExperiment` object |


<div style="page-break-after: always;"></div>



### `bluster`

Read the [`bluster` package documentation](https://rdrr.io/bioc/bluster/) and vignettes on its usage:

* [Flexible clustering for Bioconductor](https://rdrr.io/bioc/bluster/f/vignettes/clusterRows.Rmd)
* [Assorted clustering diagnostics](https://rdrr.io/bioc/bluster/f/vignettes/diagnostics.Rmd)


| Library/Package      | Piece of Code      | What it's called    | What it does  |
|----------------------|--------------------|---------------------|---------------|
| `bluster`| [`NNGraphParam()`](https://rdrr.io/bioc/bluster/man/NNGraphParam-class.html)| Graph-based clustering parameters | Set up parameters for nearest-neighbor (NN) graph-based clustering algorithms within `scran::clusterCells()` or `bluster::clusterRows()` |
| `bluster`| [`approxSilhouette()`](https://rdrr.io/bioc/bluster/man/approxSilhouette.html)| Approximate silhouette width | Calculate an approximate silhouette width for each cell given a set of clusters |
| `bluster`| [`neighborPurity()`](https://rdrr.io/bioc/bluster/man/neighborPurity.html)| Compute neighborhood purity | Calculate neighborhood purity for each cell given a set of clusters |
| `bluster`| [`bootstrapStability()`](https://rdrr.io/bioc/bluster/man/bootstrapStability.html)| Assess cluster stability by bootstrapping  | Generate cluster bootstrap replicates to estimate cluster robustness to sampling noise  |


### `patchwork`


Read the [`patchwork` package documentation](https://rdrr.io/cran/patchwork), and a [vignette on its usage](https://patchwork.data-imaginist.com/articles/patchwork.html).


| Library/Package      | Piece of Code      | What it's called    | What it does  |
|----------------------|----------------------------|--------------------------------------------|--------------------------------------------------------------|
| `patchwork` | [`wrap_plots()`](https://rdrr.io/cran/patchwork/man/wrap_plots.html)  | Wrap plots | Wrap multiple `ggplot2` objects into a single multi-panel plot  |
| `patchwork` | [`+`](https://rdrr.io/cran/patchwork/man/wrap_plots.html)  | Plot arithmetic | Place `ggplot2` objects side-by-side by adding them together with a `+`. The `patchwork` package must be loaded to use this symbol with plots |
| `patchwork` | [`/`](https://rdrr.io/cran/patchwork/man/wrap_plots.html)  | Plot arithmetic | Stack `ggplot2` objects on top of one another "dividing" them with a `/`. The `patchwork` package must be loaded to use this symbol with plots |
