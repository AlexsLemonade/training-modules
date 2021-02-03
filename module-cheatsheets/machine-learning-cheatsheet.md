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
Read the Base `R` package documentation [**here**](https://www.rdocumentation.org/packages/base/versions/3.5.1). 

| Library/Package                                  | Piece of Code                               | What it's called                               | What it does                                                                                                                                   |
|-----------------------------------------|---------------------------------------------|-------------------------|------------------------------------------------------------------------|
| Base `R`               | [`round()`](https://www.rdocumentation.org/packages/base/versions/3.5.1/topics/Round)         | Round    | Rounds the values in the object provided in the first argument to the number of decimal places specified in the second argument         |
| Base `R`               | [`identical()`](https://www.rdocumentation.org/packages/base/versions/3.5.1/topics/identical)     | Identical   | Checks if two objects are exactly equal        |
| Base `R`               | [`prcomp()`](https://www.rdocumentation.org/packages/stats/versions/3.5.1/topics/prcomp)                             | Principal Components Analysis  | Executes a principal components analysis on specified matrix or data frame                                               |   
| Base `R`               | [`rowSums()`](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/colSums)     | Row Sums  | Returns the sums of the rows in a numeric array, matrix, or data.frame | 
| Base `R`               | [`rowMeans()`](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/colSums)     | Row Means  | Returns the means of the rows in a numeric array, matrix, or data.frame | 
| Base `R`               | [`quantile()`](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/quantile)    | Sample Quantiles  | Returns the sample quantiles for a given numeric vector of data and numeric vector of probabilities  |
| Base `R`               | [`cor()`](https://www.rdocumentation.org/packages/stats/versions/3.6.1/topics/cor)        | Correlation     | Computes correlation between columns using a specified correlation method, and returns a correlation matrix        |
| Base `R`               | [`as.dist()`](https://www.rdocumentation.org/packages/base/versions/3.5.1/topics/Round)         | Distance matrix computation    | Returns a special object of class `dist`, a distance matrix used by the `hclust()` function         |
| Base `R`               | [`hclust()`](https://www.rdocumentation.org/packages/stats/versions/3.5.1/topics/hclust)    | Hierarchical Clustering  | Performs hierarchical clustering analysis on a set of dissimilarities and methods  |
| Base `R`               | [`table()`](https://www.rdocumentation.org/packages/base/versions/3.5.1/topics/table)    | Create Table       | Creates a contingency table of counts for each combination of factor levels               | 
| Base `R`               | [`duplicated()`](https://www.rdocumentation.org/packages/base/versions/3.5.1/topics/duplicated)           | Duplicated              | Returns a logical vector, where `TRUE` represents elements of the object that are duplicates               |
| Base `R`               | [`any()`](https://www.rdocumentation.org/packages/base/versions/3.5.1/topics/any)             | Any      | Checks to see if at least one of the elements are `TRUE` when given a logical vector                                      |
| Base `R`               | [`cbind()`](https://www.rdocumentation.org/packages/base/versions/3.5.1/topics/cbind)                | Column bind                          | Combines vectors, matrices, or data.frames by columns                                         |
| Base `R`             | [`pairwise.wilcox.test()`](https://www.rdocumentation.org/packages/stats/versions/3.6.1/topics/pairwise.wilcox.test) | Pairwise Wilcoxon Rank Sum Tests | Calculates the pairwise comparisons between group levels                   |

 <div style="page-break-after: always;"></div>
 
### `PLIER`
Read the `PLIER` package documentation [**here**](https://www.biorxiv.org/content/10.1101/116061v2). <br>
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
Read the `ComplexHeatmap` package documentation [**here**](https://www.rdocumentation.org/packages/ComplexHeatmap/versions/1.20.0).

| Library/Package                                  | Piece of Code                               | What it's called                               | What it does                                                                                                                                   |
|-----------------------------------------|---------------------------------------------|-------------------------|------------------------------------------------------------------------|
| `ComplexHeatmap`       | [`Heatmap()`](https://www.rdocumentation.org/packages/ComplexHeatmap/versions/1.10.2/topics/Heatmap)   | Complex Heatmap     | Constructs a heatmap whose graphics and features can be defined                                           |
| `ComplexHeatmap`       | [`HeatmapAnnotation()`](https://www.rdocumentation.org/packages/ComplexHeatmap/versions/1.10.2/topics/HeatmapAnnotation)    | Heatmap Annotation Constructor   | Creates an annotation object to be used in conjunction with a Heatmap    |

<br>
<br>
 
### `ggplot2`
Read the `ggplot2` package documentation [**here**](https://www.rdocumentation.org/packages/ggplot2/versions/3.1.0). <br>
A vignette on the usage of the `ggplot2` package can be found [**here**](https://cran.r-project.org/web/packages/ggplot2/vignettes/ggplot2-specs.html).

| Library/Package                                  | Piece of Code                               | What it's called                               | What it does                                                                                                                                   |
|-----------------------------------------|---------------------------------------------|-------------------------|------------------------------------------------------------------------|
| `ggplot2`              | [`geom_jitter()`](https://www.rdocumentation.org/packages/ggplot2/versions/3.1.0/topics/geom_jitter)          | Jittered Points  | Adds a small amount of random variation at each point’s location on a plot               |
| `ggplot2`              | [`labs()`](https://www.rdocumentation.org/packages/ggplot2/versions/3.1.0/topics/labs)        | Labels       | Sets the axis, legend, and plot labels if specified                     |
| `ggplot2`              | [`theme()`](https://www.rdocumentation.org/packages/ggplot2/versions/3.1.0/topics/theme)      | Theme        | Sets the specified non-data elements of a plot (i.e. plot title, legend spacing, text size, etc.)       |

<br>
<br>

### `tidyr`
Read the `tidyr` package documentation [**here**](https://tidyr.tidyverse.org/index.html). <br>
A vignette on the usage of the `tidyr` package can be found [**here**](https://tidyr.tidyverse.org/articles/tidy-data.html).
| Library/Package                                  | Piece of Code                               | What it's called                               | What it does                                                                                                                                   |
|-----------------------------------------|---------------------------------------------|-------------------------|------------------------------------------------------------------------|
| `tidyr`            | [`separate()`](https://tidyr.tidyverse.org/reference/separate.html)    | Separate  | Separates a character column into multiple columns with a given regular expression or numeric locations     |
| `tidyr`            | [`pivot_longer()`](https://tidyr.tidyverse.org/reference/pivot_longer.html)    | Pivot Longer  | Pivots data in a data.frame from wide to long format     |


<div style="page-break-after: always;"></div>

### Other packages and functions
Documentation for each of these packages can be accessed by clicking the package name in the table below.

| Library/Package                                  | Piece of Code                               | What it's called                               | What it does                                                                                                                                   |
|-----------------------------------------|---------------------------------------------|-------------------------|------------------------------------------------------------------------|
| [`data.table`](https://www.rdocumentation.org/packages/data.table/versions/1.11.8)       | [`fread()`](https://www.rdocumentation.org/packages/data.table/versions/1.11.8/topics/fread)    | F read           | Reads in data faster than base R                                                                         |
| [`purrr`](https://purrr.tidyverse.org/)            | [`discard()`](https://purrr.tidyverse.org/reference/keep.html)    | Discard  | Discards the given elements        |
| [`dplyr`](https://www.rdocumentation.org/packages/dplyr/versions/0.7.8)            | [`pull()`](https://www.rdocumentation.org/packages/dplyr/versions/0.7.8/topics/pull)    | Pull  | Pulls a single variable out of a given table of data        |
| [`matrixStats`](https://purrr.tidyverse.org/)            | [`rowSds()`](https://www.rdocumentation.org/packages/matrixStats/versions/0.5.0/topics/rowSds)    | Row Standard Deviations  | Returns the standard deviation estimates for each row in a matrix       |
| [`matrixStats`](https://purrr.tidyverse.org/)            | [`rowVars()`](https://www.rdocumentation.org/packages/matrixStats/versions/0.5.0/topics/rowVars)    | Row Variances  | Returns the variance estimates for each row in a matrix       |
| [`umap`](https://www.rdocumentation.org/packages/umap/versions/0.2.6.0)            | [`umap()`](https://www.rdocumentation.org/packages/umap/versions/0.2.6.0/topics/umap)    | Uniform Manifold Approximation and Projection (UMAP)  | Computes a manifold approximation and projection on a given matrix or data.frame     |
| [`ConsensusClusterPlus`](https://www.rdocumentation.org/packages/ConsensusClusterPlus/versions/1.46.0) | [`ConsensusClusterPlus()`](https://www.rdocumentation.org/packages/ConsensusClusterPlus/versions/1.36.0/topics/ConsensusClusterPlus) | Consensus Clustering             | Finds the consensus across multiple runs of the clustering algorithm                              |
| [`plotly`](https://www.rdocumentation.org/packages/plotly/versions/4.9.2.1)            | [`plot_ly()`](https://www.rdocumentation.org/packages/plotly/versions/4.9.2.1/topics/plot_ly)    | Plotly Visualization  | Initiates a plotly visualization with given R objects     |
| [`ggsignif`](https://www.rdocumentation.org/packages/ggsignif/versions/0.4.0)             | [`geom_signif()`](https://www.rdocumentation.org/packages/ggsignif/versions/0.4.0/topics/stat_signif)    | Create Significance Layer   | Adds significance information to the plot. It can be used to run statistical tests and display the significance information from those tests. We use it differently, in a way that gives us more control, in the notebook.        |
