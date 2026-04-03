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
- [`scuttle`](#scuttle)
- [`patchwork`](#patchwork)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

<div style="page-break-after: always;"></div>


### `SpatialExperiment`


Read the [`SpatialExperiment` package documentation (and e-book)](https://bioconductor.org/books/release/OSTA/) and a [vignette on its usage](https://www.bioconductor.org/packages/release/bioc/vignettes/SpatialExperiment/inst/doc/SpatialExperiment.html).

`SpatialExperiment` objects are based [`SingleCellExperiment` objects](https://bioconductor.org/books/release/OSCA/) and its existing functions, and the associated function usage is equivalent.
See [this `SingleCellExperiment` vignette](https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html) for example usage of shared functions, as well as cheatsheets for the "Introduction to scRNA-seq" (`scRNA-seq-cheatsheet.pdf`) and "Advanced scRNA-seq" (`scRNA-seq-advanced-cheatsheet.pdf`) training modules.



| Library/Package      | Piece of Code      | What it's called    | What it does  |
|----------------------|----------------------------|--------------------------------------------|--------------------------------------------------------------|
| `SpatialExperiment` | [`SpatialExperiment()`](https://www.bioconductor.org/packages/release/bioc/vignettes/SpatialExperiment/inst/doc/SpatialExperiment.html) | Spatial Experiment| Creates a `SpatialExperiment` object  |
| `SpatialExperiment` | [`spatialCoords()`](https://www.bioconductor.org/packages/release/bioc/vignettes/SpatialExperiment/inst/doc/SpatialExperiment.html#13_spatialCoords) | Spatial Experiment coordinates | Extracts and stores spatial coordinates of spots in the `SpatialExperiment` object  |
| `SpatialExperiment` | [`imgData()`](https://www.bioconductor.org/packages/release/bioc/vignettes/SpatialExperiment/inst/doc/SpatialExperiment.html#14_imgData) | Spatial 
Experiment| Creates a `SpatialExperiment` object  |




### `VisiumIO`

<!-- TODO: This package has not yet made it to rdrr -->
Read the R package [`VisiumIO` documentation](https://www.bioconductor.org/packages/release/bioc/html/VisiumIO.html) and a [vignette on its usage](https://www.bioconductor.org/packages/release/bioc/vignettes/VisiumIO/inst/doc/VisiumIO.html).

|Library/Package|Piece of Code|What it's called| What it does|
|---------------|-------------|----------------|-------------|
| `VisiumIO` | `TENxVisium()` | 10x Visium | Constructor function for importing 10x Visium data |
| `VisiumIO`| `import()`| Import | Loads 10x Visium data defined with the `TENxVisium()` constructor function |




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
| `SpotSweeper` | [`localOutliers()`](https://www.bioconductor.org/packages/release/bioc/vignettes/scuttle/inst/doc/norm.html#2_Computing_size_factors)  | Local outliers | Detect local outlier spots in a `SpatialExperiment` based on a given quality control metric |
| `SpotSweeper`| [`plotQCmetrics()`](https://www.bioconductor.org/packages/release/bioc/vignettes/scuttle/inst/doc/norm.html#2_Computing_size_factors)| Plot QC metrics | Create a spot plot highlighting outliers detected with `SpotSweeper::localOutliers()` |




### `scuttle`


Read the [`scuttle` package documentation](https://rdrr.io/bioc/scuttle/), and a [vignette on its usage](https://rdrr.io/bioc/scuttle/f/vignettes/overview.Rmd).


| Library/Package      | Piece of Code      | What it's called    | What it does  |
|----------------------|----------------------------|--------------------------------------------|--------------------------------------------------------------|
| `scuttle` | [`computeLibraryFactors()`](https://www.bioconductor.org/packages/release/bioc/vignettes/scuttle/inst/doc/norm.html#2_Computing_size_factors)  | Compute Library Factors | Returns a numeric vector of computed size factors for each spot (or cell) stored in a `SpatialExperiment` (or `SingleCellExperiment`) object. The size factor is computed as the library size of each spot/cell after scaling them to have a mean of 1 across all spots/cells |
| `scuttle`| [`logNormCounts()`](https://www.bioconductor.org/packages/release/bioc/vignettes/scuttle/inst/doc/norm.html#2_Computing_size_factors)| Normalize log counts | Returns the `SpatialExperiment` (or `SingleCellExperiment`) object with normalized expression values for each spot (cell), using the size factors stored in the object |


### `patchwork`


Read the [`patchwork` package documentation](https://rdrr.io/cran/patchwork), and a [vignette on its usage](https://patchwork.data-imaginist.com/articles/patchwork.html).


| Library/Package      | Piece of Code      | What it's called    | What it does  |
|----------------------|----------------------------|--------------------------------------------|--------------------------------------------------------------|
| `patchwork` | [`wrap_plots()`](https://rdrr.io/cran/patchwork/man/wrap_plots.html)  | Wrap plots | Wrap multiple `ggplot2` plots into a single multi-panel plot  |
