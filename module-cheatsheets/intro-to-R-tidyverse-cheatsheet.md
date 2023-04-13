# Introduction to R and Tidyverse Cheatsheet

#### The tables below consist of valuable functions and commands that will help you through this module.
##### Each table represents a different library/tool and the corresponding commands.
> Please note that these tables are not intended to tell you all the information you need to know about each command.
>
> The hyperlinks found in each piece of code will take you to the documentation for further information on the usage of each command.
Please be aware that the documentation will generally provide information about the given function's most current version (or a recent version, depending on how often the documentation site is updated).
This will usually (but not always!) match what you have installed on your machine.
If you have a different version of R or other R packages, the documentation may differ from what you have installed.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [Base `R`](#base-r)
- [`tidyverse`](#tidyverse)
  - [`dplyr`](#dplyr)
  - [`ggplot2`](#ggplot2)
  - [`readr`, `fs`, `tibble` `tidyr`](#readr-fs-tibble-tidyr)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


<div style="page-break-after: always;"></div>

### Base `R`
Read the [Base `R` documentation](https://rdrr.io/r/).

| Library/Package                    | Piece of code        | What it's called                | What it does                                                                                                      |
|------------------------------------|----------------------|---------------------------------|---------------------------------------------------------------------------------------------------------------------------------|
| Base `R`          | [`library()`](https://rdrr.io/r/base/library.html)                    | Library   | Loads and attaches additional packages to the R environment.                                                                                                                                                  |
| Base `R`          | [`<-`](https://rdrr.io/r/base/assignOps.html)                    | Assignment operator   | Assigns a name to something in the R environment.                                                                                                                                                  |
| Base `R`           | [`|>`](https://rdrr.io/r/base/pipeOp.html)        | Pipe operator         | Funnels an object from the output of one function to input of the next function |
| Base `R`          | [`c()`](https://rdrr.io/r/base/c.html)                     | Combine           | Combines values into a vector or list.                                                                                                                                                             |
| Base `R`          | [`%in%`](https://rdrr.io/r/base/match.html)                                  | "in" logical operator           | Checks if the given value(s) on the left side of the operator are in the vector or other R object defined on the right side of the operator. It returns a logical `TRUE` or `FALSE` statement. [This resource](http://www.datasciencemadesimple.com/in-operator-in-r/) also provides a helpful explanation about its usage.                                                                                                                                                         |
| Base `R`          | [`rm(x)`](https://rdrr.io/r/base/rm.html)              | Remove                | Removes object(s) `x` from your environment.                                                                                                                                                         |
| Base `R`          | [`==, <=, >=, !=`](https://rdrr.io/r/base/Comparison.html)       | Relational Operators  | These are binary operators which allow for the comparison of values in an object.                                                                                                                 |
| Base `R`          | [`str(x)`](https://rdrr.io/r/utils/str.html)               | Object Structure      | Gets a summary of the object `x` structure.                                                                                                                                                         |
| Base `R`          | [`class(x)`](https://rdrr.io/r/base/class.html)                            | Object Class      | Returns the type of the values in object `x`.                                                                |
| Base `R`          | [`nrow(x)`; `ncol(x)`](https://rdrr.io/r/base/nrow.html)                            | Number of Rows; Number of Columns    | Get the number of rows and the number of columns in an object `x`, respectively.                                                          |
| Base `R`          | [`length(x)`](https://rdrr.io/r/base/length.html)                            | Length   | Returns how long the object `x` is.                                                               |
| Base `R`          | [`min(x)`](https://rdrr.io/r/base/Extremes.html)                            | Minimum           | Returns the minimum value of all values in an object `x`.                                                           |
| Base `R`          | [`sum(x)`](https://rdrr.io/r/base/sum.html)                            | Sum               | Returns the sum of all values (values must be integer, numeric, or logical) in object `x`.                                                           |
| Base `R`          | [`mean(x)`](https://rdrr.io/r/base/mean.html)                            | Mean               | Returns the arithmetic mean of all values (values must be integer or numeric) in object `x` or logical vector `x`.                                                           |
| Base `R`          | [`log(x)`](https://rdrr.io/r/base/Log.html)                            | Logarithm        | Gives the natural logarithm of object `x`. `log2(x)` can be used to give the logarithm of the object in base 2. Or the base can be specified as an argument.                                                            |
| Base `R`          | [`head()`; `tail()`](https://rdrr.io/r/utils/head.html)             | Head; Tail               | Returns the top 6 (`head()`) or bottom 6 (`tail()`) rows of an object in the environment by default. You can specify how many rows you want by including the `n = ` argument.                                                        |
| Base `R`          | [`factor(x)` or `as.factor(x)`](https://rdrr.io/r/base/factor.html)       | Factor             | Coerces object x into a factor (which is used to represent categorical data). This function can be used to coerce object `x` into other data types, i.e., `as.character`, `as.numeric`, `as.data.frame`, `as.matrix`, etc. |
| Base `R`          | [`levels(x)`](https://rdrr.io/r/base/levels.html)                | Levels attributes    | Returns or sets the value of the levels in an object `x`.                                                                                                                                                               |
| Base `R`          | [`summary(x)`](https://rdrr.io/r/base/summary.html)                | Object summary    | Returns a summary of the values in object `x`.                                                                                                                                                               |
| Base `R`          | [`data.frame()`](https://rdrr.io/r/base/data.frame.html)   | Data Frame            | Creates a data frame where the named arguments will be the same length.                                                                                             |
| Base `R`          | [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html) | Session Information   | Returns the R version information, the OS, and the attached packages in the current R session.                                                                                                 |
| Base `R`          | [`file.path()`](https://rdrr.io/r/base/file.path.html)             | File path | Constructs the path to a desired file.                                                                                                                                                              |
| Base `R`          | [`dir()`](https://rdrr.io/r/base/list.files.html)             | Directory | Lists the names of the files and/or directories in the named directory.                                                                                                                                                              |
| Base `R`          | [`getwd()`](https://rdrr.io/r/base/getwd.html)             | Get working directory | Finds the current working directory.                                                                                                                                                              |
| Base `R`          | [`setwd()`](https://rdrr.io/r/base/getwd.html)             | Set working directory | Changes the current working directory.                                                                                                                                                            |
| Base `R`          | [`dir.exists()`](https://rdrr.io/r/base/files2.html)         | Directory exists      | Checks the file path to see if the directory exists there.                                                                                                                                         |
| Base `R`          | [`dir.create()`](https://rdrr.io/r/base/files2.html)             | Create directory      | Creates a directory at the specified path.                                                                                                                                                   |
| Base `R`          | [`apply()`](https://rdrr.io/r/base/apply.html)             | Apply                 | Returns a vector or list of values after applying a specified function to values in each row/column of an object.                                                                                  |
| Base `R`          | [`round()`](https://rdrr.io/r/base/Round.html)             | Round                 | Rounds the values of an object to the specified number of decimal places (default is 0).                                                                                                          |
| Base `R`          | [`names()`](https://rdrr.io/r/base/names.html)                       | Names                     | Gets or sets the names of an object.                                      |
| Base `R`          | [`colnames()`](https://rdrr.io/r/base/colnames.html)           | Column names              | Gets or sets the column names of a matrix or data frame.                  |
| Base `R`                | [`all.equal()`](https://rdrr.io/r/base/all.equal.html)               | All equal                 | Checks if two R objects are nearly equal.                                 |
| Base `R`                | [`all()`](https://rdrr.io/r/base/all.html)               | All               | Checks if all of the values are `TRUE` in a logical vector.                                 |
| Base `R`               | [`t()`](https://rdrr.io/r/base/t.html)                                        | Transpose       | Returns the transpose of a matrix or data frame. If given a data frame, returns a matrix.                                                          |

<div style="page-break-after: always;"></div>

### `tidyverse`
Read the [`tidyverse` package documentation](https://tidyverse.tidyverse.org/), as well as the [philosophy behind the `tidyverse`](https://tidyverse.tidyverse.org/articles/manifesto.html).

#### `dplyr`
Read the [`dplyr` package documentation](https://dplyr.tidyverse.org/), and a [vignette on its usage](https://dplyr.tidyverse.org/articles/dplyr.html).

| Library/Package                                 | Piece of code                        | What it's called                       | What it does                                                                                                                                                                                      |
|-------------------------------------------------|--------------------------------------|----------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `dplyr`/`magrittr`           | [`%>%`](https://rdrr.io/cran/magrittr/man/pipe.html)        | Pipe operator         | Funnels an object from the output of one function to input of the next function (used like the base pipe `|>` found in `R` 4.1 and later, but can be used in earlier versions of `R`)   |
| `dplyr`           | [`filter()`](https://dplyr.tidyverse.org/reference/filter.html)            | Filter                | Returns a subset of rows matching the conditions of the specified logical argument                                                                                                                |
| `dplyr`           | [`arrange()`](https://dplyr.tidyverse.org/reference/arrange.html)          | Arrange               | Reorders rows in ascending order. `arrange(desc())` would reorder rows in descending order.                                                                                                         |
| `dplyr`           | [`select()`](https://dplyr.tidyverse.org/reference/select.html)            | Select                | Selects columns that match the specified argument                                                                                                                                                |
| `dplyr`           | [`mutate()`](https://dplyr.tidyverse.org/reference/mutate.html)            | Mutate                | Adds a new column that is a function of existing columns                                                                                                                                          |
| `dplyr`           | [`summarize()`](https://dplyr.tidyverse.org/reference/summarise.html)      | Summarize             | Summarizes multiple values in an object into a single value. This function can be used with other functions to retrieve a single output value for the grouped values. `summarize` and `summarise` are synonyms in this package. However, note that this function does not work in the same manner as the base R `summary` function.                                                                                                   |
| `dplyr`           | [`rename()`](https://dplyr.tidyverse.org/reference/rename.html)            | Rename                | Renames designated columns while keeping all variables of the data.frame                                                                                                                                                                        |
| `dplyr`           | [`group_by()`](https://dplyr.tidyverse.org/reference/group_by.html)        | Group By              | Groups data into rows that contain the same specified value(s)                                                                                                                                    |
| `dplyr`           | [`inner_join()`](https://dplyr.tidyverse.org/reference/mutate-joins.html)          | Inner Join            | Joins data from two data frames, retaining only the rows that are in both datasets.    |

<div style="page-break-after: always;"></div>


#### `ggplot2`
Read the [`ggplot2` package documentation](https://ggplot2.tidyverse.org/), an [overall reference for `ggplot2` functions](https://ggplot2.tidyverse.org/reference/index.html), and a [vignette on the usage of the `ggplot2` aesthetics](https://ggplot2.tidyverse.org/articles/ggplot2-specs.html).
Additional vignettes are available from the "Articles" dropdown menu on this webpage.

| Library/Package                      | Piece of code                                   | What it's called       | What it does                                                                                                                                                                                      |
|--------------------------------------|-----------------------------------------------------|--------------------------|--------------------------------------------------------------------------------------------------------------------------------------|
| `ggplot2`         | [`ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html)                   | GG Plot               | Begins a plot that is finished by adding layers.        |
| `ggplot2`         | [`aes()`](https://ggplot2.tidyverse.org/reference/aes.html)                   | Aesthetic Mappings              | Designates how variables in the data object are mapped to the visual properties of the ggplot.        |
| `ggplot2`         | [`geom_boxplot()`](https://ggplot2.tidyverse.org/reference/geom_boxplot.html)       | Boxplot               | Creates a boxplot when added as a layer to a `ggplot()` object.            |
| `ggplot2`         | [`geom_density()`](https://ggplot2.tidyverse.org/reference/geom_density.html)       | Density Plot               | Creates a smoothed plot when added as a layer to a `ggplot()` object based on the computed density estimate.            |
| `ggplot2`         | [`geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html)       | Scatterplot               | Creates a scatterplot when added as a layer to a `ggplot()` object.            |
| `ggplot2`         | [`geom_line()`](https://ggplot2.tidyverse.org/reference/geom_path.html)       | Line plot               | Creates a line plot when added as a layer to a `ggplot()` object by connecting the points in order of the x axis variable.            |
| `ggplot2`         | [`geom_hline()`](https://ggplot2.tidyverse.org/reference/geom_abline.html)       | Horizontal line              | Annotates a plot with a horizontal line when added as a layer to a `ggplot()` object
| `ggplot2`         | [`geom_vline()`](https://ggplot2.tidyverse.org/reference/geom_abline.html)       | Vertical line              | Annotates a plot with a vertical line when added as a layer to a `ggplot()` object
| `ggplot2`         | [`theme_classic()`](https://ggplot2.tidyverse.org/reference/ggtheme.html)               | Classic Theme   | Displays `ggplot` without gridlines. The [`ggtheme` documentation](https://ggplot2.tidyverse.org/reference/ggtheme.html) has descriptions on additional themes that can be used.                                                                         |
| `ggplot2`               | [`labs()`](https://ggplot2.tidyverse.org/reference/labs.html)           | Labels   | Modify labels (axis, title, legends) on a `ggplot()` object.   |
| `ggplot2`               | [`xlab()`; `ylab()`; `ggtitle()`](https://ggplot2.tidyverse.org/reference/labs.html)           |X Axis Labels; Y Axis Labels; GG Title                    | Alternative individual functions to add individual plot labels: x-axis, y-axis, and title, respectively.|
| `ggplot2`            | [`facet_wrap()`](https://ggplot2.tidyverse.org/reference/facet_wrap.html)                     | Facet Wrap      | Plots individual graphs using specified variables to subset the data.                                                                   |
| `ggplot2`               | [`ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)                 | GG Save                                         | Saves the last plot in working directory.                                 |
| `ggplot2`               | [`last_plot()`](https://ggplot2.tidyverse.org/reference/last_plot.html)           | Last plot                                       | Returns the last plot produced.                                           |


<div style="page-break-after: always;"></div>

#### `readr`, `fs`, `tibble` `tidyr`
Read the [`readr` package documentation](https://readr.tidyverse.org/) and a [vignette on its usage](https://readr.tidyverse.org/articles/readr.html).
Read the [`fs` package documentation](https://fs.r-lib.org/).
Read the [`tibble` package documentation](https://tibble.tidyverse.org/) and a [vignette on its usage](https://tibble.tidyverse.org/articles/tibble.html).
Read the [`tidyr` package documentation](https://tidyr.tidyverse.org/) and a [vignette on its usage](https://tidyr.tidyverse.org/articles/tidy-data.html).

| Library/Package                      | Piece of code                                   | What it's called       | What it does                                                                                                                                                                                      |
|--------------------------------------|-----------------------------------------------------|---------------------------|-----------------------------------------------------------------------------------------------------------------------------------|
| `readr`           | [`read_tsv()`](https://readr.tidyverse.org/reference/read_delim.html)           | Read TSV             | Reads in a TSV file from a specified file path. This function can be tailored to read in other common types of files, e.g. `read_csv()`, `read_rds()`, etc.                                          |
| `fs`           | [`dir_create()`](https://fs.r-lib.org/reference/create.html)           | Create directory             | Create a directory, unless the directory already exists.
| `tibble`          | [`column_to_rownames()`](https://tibble.tidyverse.org/reference/rownames.html) | Column to Rownames    | Transforms an existing column called by a string into the rownames.                            |
| `tibble`          | [`rownames_to_column()`](https://tibble.tidyverse.org/reference/rownames.html)      | Rownames to Column    | Transforms the rownames of a data frame into a column (which is added to the start of the data frame).  The string supplied as an argument will be the name of the new column.            |
| `tidyr`          | [`pivot_longer()`](https://tidyr.tidyverse.org/reference/pivot_longer.html)                  | Pivot Longer    | Lengthens a data frame by increasing the number of rows and decreasing the number of columns.              |
