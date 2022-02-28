# Machine Learning Cheatsheet

#### The tables below consist of valuable functions or commands that will help you through this module. 
##### Each table represents a different library/tool and its corresponding commands. 
> Please note that these tables are not intended to tell you all the information you need to know about each command. 
> 
> The hyperlinks found in each piece of code will take you to the documentation for further information on the usage of each command. 
 
<br>
<br>

### AnnotationDbi 
Read the `AnnotationDbi` package vignette [**here**](http://bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf).

| Library/Package                                  | Piece of Code                               | What it's called                               | What it does                                                                                                                                   |
|-----------------------------------------|---------------------------------------------|-------------------------|------------------------------------------------------------------------|
| `AnnotationDbi`        | [`keytypes()`](https://www.rdocumentation.org/packages/AnnotationDbi/versions/1.44.0/topics/AnnotationDb-objects)      | Keytypes     | Returns a character vector of column names/keytypes (e.g., type of gene identifiers) available in an `AnnotationDbi` package.                                              |
| `AnnotationDbi`       | [`mapIDs()`](https://www.rdocumentation.org/packages/AnnotationDbi/versions/1.44.0/topics/AnnotationDb-objects) | Mapped IDs       | Extracts the mapped ids for a set of keys (e.g., gene identifiers) of a specific keytype                                    |                                  

 <div style="page-break-after: always;"></div>
 
### Base `R`
Read the Base `R` documentation [**here**](https://rdrr.io/r/).

| Library/Package                                  | Piece of Code                               | What it's called                               | What it does                                                                                                                                   |
|-----------------------------------------|---------------------------------------------|-------------------------|------------------------------------------------------------------------|
| Base `R`               | [`round()`](https://rdrr.io/r/base/Round.html)         | Round    | Rounds the values in the object provided in the first argument to the number of decimal places specified in the second argument         |
| Base `R`               | [`identical()`](https://rdrr.io/r/base/identical.html)     | Identical   | Checks if two objects are exactly equal        |
| Base `R`               | [`prcomp()`](https://rdrr.io/r/stats/prcomp.html)                             | Principal Components Analysis  | Executes a principal components analysis on specified matrix or data frame                                               |   
| Base `R`               | [`rowSums()`](https://rdrr.io/r/base/colSums.html)     | Row Sums  | Returns the sums of the rows in a numeric array, matrix, or data.frame | 
| Base `R`               | [`rowMeans()`](https://rdrr.io/r/base/colSums.html)     | Row Means  | Returns the means of the rows in a numeric array, matrix, or data.frame | 
| Base `R`               | [`quantile()`](https://rdrr.io/r/stats/quantile.html)    | Sample Quantiles  | Returns the sample quantiles for a given numeric vector of data and numeric vector of probabilities  |
| Base `R`               | [`cor()`](https://rdrr.io/r/stats/cor.html)        | Correlation     | Computes correlation between columns using a specified correlation method, and returns a correlation matrix        |
| Base `R`               | [`as.dist()`](https://rdrr.io/r/stats/dist.html)         | Distance matrix computation    | Returns a special object of class `dist`, a distance matrix used by the `hclust()` function         |
| Base `R`               | [`hclust()`](https://rdrr.io/r/stats/hclust.html)    | Hierarchical Clustering  | Performs hierarchical clustering analysis on a set of dissimilarities and methods  |
| Base `R`               | [`table()`](https://rdrr.io/r/base/table.html)    | Create Table       | Creates a contingency table of counts for each combination of factor levels               | 
| Base `R`               | [`duplicated()`](https://rdrr.io/r/base/duplicated.html)      | Duplicated   | Returns a logical vector, where `TRUE` represents elements of the object that are duplicates               |
| Base `R`               | [`any()`](https://rdrr.io/r/base/any.html) | Any  | Checks to see if at least one of the elements are `TRUE` when given a logical vector                  |
| Base `R`               | [`cbind()`](https://rdrr.io/r/base/cbind.html)  | Column bind    | Combines vectors, matrices, or data.frames by columns                                   |
| Base `R`             | [`pairwise.wilcox.test()`](https://rdrr.io/r/stats/pairwise.wilcox.test.html) | Pairwise Wilcoxon Rank Sum Tests | Calculates the pairwise comparisons between group levels    |

 <div style="page-break-after: always;"></div>
 
### `PLIER`
Read the `PLIER` package documentation [**here**](https://rdrr.io/github/wgmao/PLIER/f/README.md), and see its corresponding paper [**here**](https://www.nature.com/articles/s41592-019-0456-1). <br>
A `PLIER` package vignette can be found [**here**](https://github.com/wgmao/PLIER/blob/master/vignettes/vignette.pdf) and can also serve as documentation for the commands in the table below. 

| Library/Package                                  | Piece of Code                               | What it's called                               | What it does                                                                                                                                   |
|-----------------------------------------|---------------------------------------------|-------------------------|------------------------------------------------------------------------|
| `PLIER`               | `combinePaths()`         | Combine Pathways                 | Combines the pathway data obtained from `PLIER` and returns the result as a matrix                                                            |
| `PLIER`                | `commonRows()`           | Common Rows                      | Determines the rows (genes) that are common to the specified data matrices and returns them as a character vector                             |
| `PLIER`                | `rowNorm()`              | Row Normalize                    | Normalizes each row (gene) by z-scoring the expression values                                                                                 |
| `PLIER`                | `num.pc()`               | Number of Principal Components   | Returns the number of significant principal components                                                                                        |
| `PLIER`                | `PLIER()`                | Main PLIER Function              | Main function of the Pathway-Level Information ExtractoR.                           |
| `PLIER`                | `plotU()`                | Plot U Matrix                    | Plots the U matrix obtained from the `PLIER` function results, allowing insight into the pathways or cell types captured by the latent variables |

 <div style="page-break-after: always;"></div>
 
### `ComplexHeatmap`
Read the `ComplexHeatmap` package documentation [**here**](https://rdrr.io/bioc/ComplexHeatmap/).

| Library/Package                                  | Piece of Code                               | What it's called                               | What it does                                                                                                                                   |
|-----------------------------------------|---------------------------------------------|-------------------------|------------------------------------------------------------------------|
| `ComplexHeatmap`       | [`Heatmap()`](https://rdrr.io/bioc/ComplexHeatmap/man/Heatmap.html)   | Complex Heatmap     | Constructs a heatmap whose graphics and features can be defined                                           |
| `ComplexHeatmap`       | [`HeatmapAnnotation()`](https://rdrr.io/bioc/ComplexHeatmap/man/HeatmapAnnotation.html)    | Heatmap Annotation Constructor   | Creates an annotation object to be used in conjunction with a Heatmap    |

<br>
<br>
 
### `ggplot2`
Read the `ggplot2` package documentation [**here**](https://ggplot2.tidyverse.org/). <br>
A vignette on the usage of the `ggplot2` aesthetics can be found [**here**](https://ggplot2.tidyverse.org/articles/ggplot2-specs.html). Additional vignettes are available from the "Articles" dropdown menu on this webpage.

| Library/Package                                  | Piece of Code                               | What it's called                               | What it does                                                                                                                                   |
|-----------------------------------------|---------------------------------------------|-------------------------|------------------------------------------------------------------------|
| `ggplot2`              | [`geom_jitter()`](https://ggplot2.tidyverse.org/reference/geom_jitter.html)          | Jittered Points  | Adds a small amount of random variation at each point’s location on a plot               |
| `ggplot2`              | [`labs()`](https://ggplot2.tidyverse.org/reference/labs.html)        | Labels       | Sets the axis, legend, and plot labels if specified                     |
| `ggplot2`              | [`theme()`](https://ggplot2.tidyverse.org/reference/theme.html)      | Theme        | Modify the theme to sets the specified non-data elements of a plot (i.e. plot title, legend spacing, text size, etc.)       |

<br>
<br>

### `tidyr`
Read the `tidyr` package documentation [**here**](https://tidyr.tidyverse.org/index.html). <br>
A vignette on the usage of the `tidyr` package can be found [**here**](https://tidyr.tidyverse.org/articles/tidy-data.html).


| Library/Package                                  | Piece of Code                               | What it's called                               | What it does                                                                                                                                   |
|-----------------------------------------|---------------------------------------------|-------------------------|------------------------------------------------------------------------|
| `tidyr`            | [`separate()`](https://tidyr.tidyverse.org/reference/separate.html)    | Separate  | Separates a character column into multiple columns with a given regular expression or numeric locations     |
| `tidyr`            | [`pivot_longer()`](https://tidyr.tidyverse.org/reference/pivot_longer.html)    | Pivot Longer  | Pivots data in a data.frame from wide to long format     |

<br>
<br> 

<div style="page-break-after: always;"></div>

### Other packages and functions
Documentation for each of these packages can be accessed by clicking the package name in the table below.

| Library/Package                                  | Piece of Code                               | What it's called                               | What it does             |
|-----------------------------------------|---------------------------------------------|-------------------------|-----------------------------------------------|
| [`data.table`](https://rdrr.io/cran/data.table/)       | [`fread()`](https://rdrr.io/cran/data.table/man/fread.html)    | F read           | Reads in data faster than base R                                                                         |
| [`purrr`](https://purrr.tidyverse.org/)            | [`discard()`](https://purrr.tidyverse.org/reference/keep.html)    | Discard  | Discards the given elements        |
| [`dplyr`](https://dplyr.tidyverse.org/)            | [`pull()`](https://dplyr.tidyverse.org/reference/pull.html)    | Pull  | Pulls a single variable out of a given table of data        |
| [`matrixStats`](https://rdrr.io/rforge/matrixStats/man/matrixStats-package.html)            | [`rowSds()`](https://rdrr.io/rforge/matrixStats/man/rowSds.html)    | Row Standard Deviations  | Returns the standard deviation estimates for each row in a matrix       |
| [`matrixStats`](https://rdrr.io/rforge/matrixStats/man/matrixStats-package.html)            | [`rowVars()`](https://rdrr.io/rforge/matrixStats/man/rowVars.html)    | Row Variances  | Returns the variance estimates for each row in a matrix       |
| [`umap`](https://rdrr.io/cran/umap/)            | [`umap()`](https://rdrr.io/cran/umap/man/umap.html)    | Uniform Manifold Approximation and Projection (UMAP)  | Computes a manifold approximation and projection on a given matrix or data.frame     |
| [`ConsensusClusterPlus`](https://rdrr.io/bioc/ConsensusClusterPlus/) | [`ConsensusClusterPlus()`](https://rdrr.io/bioc/ConsensusClusterPlus/man/ConsensusClusterPlus.html) | Consensus Clustering             | Finds the consensus across multiple runs of the clustering algorithm                              |
| [`plotly`](https://rdrr.io/cran/plotly/)            | [`plot_ly()`](https://rdrr.io/cran/plotly/man/plot_ly.html)    | Plotly Visualization  | Initiates a plotly visualization with given R objects     |
| [`ggsignif`](https://rdrr.io/cran/ggsignif/)             | [`geom_signif()`](https://rdrr.io/cran/ggsignif/man/stat_signif.html)    | Create Significance Layer   | Adds significance information to the plot. It can be used to run statistical tests and display the significance information from those tests. We use it differently, in a way that gives us more control, in the notebook.        |
