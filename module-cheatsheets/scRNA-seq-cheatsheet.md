# scRNA-seq Cheatsheet

#### The tables below consist of valuable functions or commands that will help you through this module.
##### Each table represents a different library/tool and its corresponding commands.
> Please note that these tables are not intended to tell you all the information you need to know about each command.
>
> The hyperlinks found in each piece of code will take you to the documentation for further information on the usage of each command.

<div style="page-break-after: always;"></div>

### Base `R`

Read the Base `R` documentation [**here**](https://www.rdocumentation.org/packages/base/versions/3.5.1)

|Library/Package|Piece of Code|What it's called| What it does|
|---------------|-------------|----------------|-------------|
| Base `R`| [`rowSums()`](https://www.rdocumentation.org/packages/raster/versions/3.3-13/topics/rowSums)| Row sums|Calculates sums for each row|
| Base `R`| [`colSums()`](https://www.rdocumentation.org/packages/raster/versions/3.3-13/topics/colSums)| column sums| Calculates sums for each column|
| Base `R`| [`t()`](https://www.rdocumentation.org/packages/base/versions/3.5.1/topics/t)  |Transpose| Returns the transpose of a matrix or data frame|
| Base `R`| [`prcomp()`](https://www.rdocumentation.org/packages/stats/versions/3.5.1/topics/prcomp)| Principal Components Analysis|Executes a principal components analysis on specified matrix or data frame|
| Base `R` [`<-function(x) { <code> }`](http://adv-r.had.co.nz/Functions.html) |Create a function| Creates a function that would take the defined parameters as input and execute the commands within the curly braces  |


### `ggplot2`

Read the `ggplot2` documentation [**here**](https://ggplot2.tidyverse.org/reference/)

| Library/Package| Piece of Code| What it's called| What it does |
| Base `R`| [`<-function(x) { <code> }`](http://adv-r.had.co.nz/Functions.html)| Create a function| Creates a function that would take the defined parameters as input and execute the commands within the curly braces  |

### `scran`, `scater`, `SingleCellExperiment`
Read the `scran` package documentation [**here**](https://www.rdocumentation.org/packages/scran/versions/1.10.2), and a vignette on its usage [**here**](https://bioconductor.org/packages/devel/bioc/vignettes/scran/inst/doc/scran.html).

Read the `scater` package documentation [**here**](https://www.rdocumentation.org/packages/scater/versions/1.10.1), and a vignette on its usage [**here**](https://github.com/davismcc/scater/blob/master/vignettes/vignette-intro.Rmd).

Read the `SingleCellExperiment` package documentation [**here**](https://osca.bioconductor.org/), and a vignette on its usage [**here**](https://bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html).

In addition to the links above, [Amezquita et al.](https://www.biorxiv.org/content/biorxiv/early/2019/03/27/590562.full.pdf) is a useful paper on the single-cell analysis workflow involving the `SingleCellExperiment` object.

<br>

| Library/Package      | Piece of Code      | What it's called    | What it does  |
|----------------------|----------------------------|--------------------------------------------|--------------------------------------------------------------|
| `SingleCellExperiment` | [`SingleCellExperiment()`](https://bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html#2_creating_singlecellexperiment_instances) | Single Cell Experiment      | Creates a `SingleCellExperiment` object  |
| `SingleCellExperiment`| [`logcounts()`](https://bioconductor.org/packages/devel/bioc/vignettes/scater/inst/doc/vignette-intro.html#3_calculating_a_variety_of_expression_values)   | Logcounts      | Stores or extracts log-transformed single-cell experiment count data as an assay of the `SingleCellExperiment` object       |
| `SingleCellExperiment`| [`rowData()`](https://bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html#5_extracting_coldata_and_rowdata)  | Row Data   | Extracts and stores gene-level metadata that describes features of the `SingleCellExperiment` object|
| `scran` | [`quickCluster()`](https://www.rdocumentation.org/packages/scran/versions/1.0.3/topics/Quick%20clustering) | Quick Clustering     | Groups similar cells into clusters which are stored in the `SingleCellExperiment` object and are used for the calculation of size factors by `scran::computeSumFactors`  |
| `scran` | [`computeSumFactors()`](https://www.rdocumentation.org/packages/scran/versions/1.0.3/topics/Deconvolution%20Methods)  | Compute Sum Factors| Returns a numeric vector of computed sum factors for each cell cluster stored in the `SingleCellExperiment` object. The cluster-based size factors are deconvolved into cell-based size factors that are stored in the `SingleCellExperiment` object and used by the `scran::normalize` function for the normalization of each cell's gene expression profile.    |
| `scater`| [`normalize()`](https://www.rdocumentation.org/packages/scater/versions/1.0.4/topics/normalize)| Normalize | Returns the `SingleCellExperiment` object with normalized expression values for each cell, using the size factors stored in the object  |


<div style="page-break-after: always;"></div>

### `alevinQC`, `colorblindr`, `Rtsne`, `tibble`
Documentation for each of these packages can be accessed by clicking the package name in the table below.

| Library/Package|Piece of Code| What it's called| What it does  |
|----------------|-------------|-----------------|---------------|
| [`alevinQC`](http://www.bioconductor.org/packages/devel/bioc//vignettes/alevinQC/inst/doc/alevinqc.html) | [`alevinQCReport()`](http://www.bioconductor.org/packages/devel/bioc//vignettes/alevinQC/inst/doc/alevinqc.html#generate-qc-report) | Alevin QC Report | Produces a QC (quality check) report from the alevin output |
| [`colorblindr`](https://www.rdocumentation.org/packages/colorblindr/versions/0.1.0)|[`scale_color_OkabeIto()`](https://www.rdocumentation.org/packages/colorblindr/versions/0.1.0/topics/scale_colour_OkabeIto) | OkabeIto Color Scale | When added as a layer to a plot, makes the plot colorblind friendly |

<br>
<br>

### Salmon
Read the Salmon documentation [**here**](https://salmon.readthedocs.io/en/latest/salmon.html).

| Piece of Code| What it's called| What it does |  
|--------------|-----------------|--------------|
| [`salmon alevin`](https://salmon.readthedocs.io/en/latest/alevin.html)    | Salmon Alevin     | Runs the Alevin quantification from the command line  |
