# Intro to R and Tidyverse Cheatsheet

#### The tables below consist of valuable functions and commands that will help you through this module. 
##### Each table represents a different library/tool and the corresponding commands. 
> Please note that these tables are not intended to tell you all the information you need to know about each command. 
> 
> The hyperlinks found in each piece of code will take you to the documentation for further information on the usage of each command. 

<div style="page-break-after: always;"></div>

### Base `R`
Read the Base `R` documentation [**here**](https://www.rdocumentation.org/packages/base/versions/3.5.1).
                                                                                                                 
| Library/Package                    | Piece of code        | What it's called                | What it does                                                                                                      |                                                                                 
|------------------------------------|----------------------|---------------------------------|---------------------------------------------------------------------------------------------------------------------------------|
| Base `R`          | [`<-`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/assignOps.html)                    | Assignment operator   | Assigns a name to something in the R environment                                                                                                                                                  |
| Base `R`          | [`c()`](https://www.rdocumentation.org/packages/base/versions/3.5.1/topics/c)                     | Concatenate           | Combines values into a vector or list                                                                                                                                                             |
| Base `R`          | [`rm(x)`](https://www.rdocumentation.org/packages/base/versions/3.5.1/topics/remove)              | Remove                | Removes object(s) `x` from your environment                                                                                                                                                         |
| Base `R`          | [`==, <=, >=, !=`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/Comparison.html)       | Relational Operators  | These are binary operators which allow for the comparison of values in an object.                                                                                                                 |
| Base `R`          | [`str(x)`](https://www.rdocumentation.org/packages/utils/versions/3.5.1/topics/str)               | Object Structure      | Gets a summary of the object `x` structure.                                                                                                                                                         |
| Base `R`          | [`head()`](https://www.rdocumentation.org/packages/utils/versions/3.5.1/topics/head)             | Head                  | Returns the top 6 rows of an object in the environment by default. You can specify how many rows you want by including the `n = `argument.                                                        |
| Base `R`          | [`tail()`](https://www.rdocumentation.org/packages/utils/versions/3.5.1/topics/head)             | Tail                  | Returns the bottom 6 rows of an object in the environment by default. You can specify how many rows you want by including the `n =` argument.                                                                  |
| Base `R`          | [`as.factor(x)`](https://www.rdocumentation.org/packages/base/versions/3.5.1/topics/factor)       | As Factor             | Coerces object x into a factor (which is used to represent categorical data). This function can be used to coerce object `x` into other data types, i.e., `as.character`, `as.data.frame`, `as.matrix`, etc. |
| Base `R`          | [`dim(x)`](https://www.rdocumentation.org/packages/base/versions/3.5.1/topics/dim)                | Dimensions summary    | Returns the dimensions of object `x`.                                                                                                                                                               |
| Base `R`          | [`data.frame()`](https://www.rdocumentation.org/packages/base/versions/3.5.1/topics/data.frame)   | Data Frame            | Creates a data.frame where the named arguments will be the same length                                                                                             |
| Base `R`          | <a href = "https://www.rdocumentation.org/packages/base/versions/3.5.1/topics/Extract">`[[x]]`</a>   |  Double bracket       | Extracts `x` from a nested list                                |
| Base `R`          | [`sessionInfo()`](https://www.rdocumentation.org/packages/utils/versions/3.5.1/topics/sessionInfo) | Session Information   | Returns the R version information, the OS, and the attacked packages in the current R session                                                                                                 |
| Base `R`          | [`getwd()`](https://www.rdocumentation.org/packages/base/versions/3.5.1/topics/getwd)             | Get working directory | Finds the current working directory.                                                                                                                                                              |
| Base `R`          | [`setwd()`](https://www.rdocumentation.org/packages/base/versions/3.5.1/topics/getwd)             | Set working directory | Changes the current working directory.                                                                                                                                                            |
| Base `R`          | [`dir.exists()`](https://www.rdocumentation.org/packages/logcondens.mode/versions/1.0.1/topics/dir.exists)         | Directory exists      | Checks the file path to see if the directory exists there                                                                                                                                         |
| Base `R`          | [`dir.create()`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/files2.html)             | Create directory      | Creates a directory at the specified file path.                                                                                                                                                   |
| Base `R`          | [`apply()`](https://www.rdocumentation.org/packages/base/versions/3.5.1/topics/apply)             | Apply                 | Returns a vector or list of values after applying a specified function to values in each row/column of an object                                                                                  |
| Base `R`          | [`round()`](https://www.rdocumentation.org/packages/base/versions/3.5.1/topics/Round)             | Round                 | Rounds the values of an object to the specified number of decimal places (default is 0).                                                                                                          |
| Base `R`          | [`rnorm()`](https://www.rdocumentation.org/packages/compositions/versions/1.40-2/topics/rnorm)    | R Norm                | Generates a vector of random numbers from a standard normal distribution                      |

<div style="page-break-after: always;"></div>
                                                                                                
### `dplyr`
Read the `dplyr` package documentation [**here**](https://www.rdocumentation.org/packages/dplyr/versions/0.7.8). <br>
A vignette on the usage of the `dplyr` package can be found [**here**](https://www.rdocumentation.org/packages/dplyr/versions/0.7.8/vignettes/dplyr.Rmd).

| Library/Package                                 | Piece of code                        | What it's called                       | What it does                                                                                                                                                                                      |
|-------------------------------------------------|--------------------------------------|----------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `dplyr`           | [`%>%`](https://www.rdocumentation.org/packages/magrittr/versions/1.5/topics/%25%26gt%3B%25)        | Pipe operator         | Funnels a data.frame through tidyverse operations                                                                                                                                                 |
| `dplyr`           | [`filter()`](https://www.rdocumentation.org/packages/dplyr/versions/0.7.8/topics/filter)            | Filter                | Returns a subset of rows matching the conditions of the specified logical argument                                                                                                                |
| `dplyr`           | [`arrange()`](https://www.rdocumentation.org/packages/dplyr/versions/0.7.8/topics/arrange)          | Arrange               | Reorders rows in ascending order. arrange(desc()) would reorder rows in descending order.                                                                                                         |
| `dplyr`           | [`select()`](https://www.rdocumentation.org/packages/dplyr/versions/0.7.8/topics/select)            | Select                | Selects columns that match the specified argument                                                                                                                                                |
| `dplyr`           | [`mutate()`](https://www.rdocumentation.org/packages/dplyr/versions/0.7.8/topics/mutate)            | Mutate                | Adds a new column that is a function of existing columns                                                                                                                                          |
| `dplyr`           | [`summarise()`](https://www.rdocumentation.org/packages/dplyr/versions/0.7.8/topics/summarise)      | Summarise             | Summarises multiple values in an object into a single value. This function can be used with other functions to retrieve a single output value for the grouped values. `summarize` and `summarise` are synonyms in this package.                                                                                                    |
| `dplyr`           | [`rename()`](https://www.rdocumentation.org/packages/dplyr/versions/0.7.8/topics/select)            | Rename                | Renames designated columns while keeping all variables of the data.frame                                                                                                                                                                        |
| `dplyr`           | [`group_by()`](https://www.rdocumentation.org/packages/dplyr/versions/0.7.8/topics/group_by)        | Group By              | Groups data into rows that contain the same specified value(s)                                                                                                                                    |
| `dplyr`           | [`inner_join()`](https://www.rdocumentation.org/packages/dplyr/versions/0.7.8/topics/join)          | Inner Join            | Joins data from two data frames, retaining only the rows that are in both datasets.    |

<div style="page-break-after: always;"></div>
                                                                                                                                                                                                                 
### Docker 
Read documentation on the Docker container [**here**](https://www.docker.com/resources/what-container). 

Piece of code                                   | What it's called                                | What it does                           |
|-------------------------------------------------|------------------------|--------------------------------------------------------------------------------------|
| [`docker pull`](https://docs.docker.com/engine/reference/commandline/pull/) | Docker Pull      | Pulls an image from a docker container |
| [`docker run`](https://docs.docker.com/v17.12/engine/reference/run/)  | Docker Run             | Runs processes in a docker container   |


### `ggplot2`
Read the `ggplot2` package documentation [**here**](https://www.rdocumentation.org/packages/ggplot2/versions/3.1.0). <br>
A vignette on the usage of the `ggplot2` package can be found [**here**](https://cran.r-project.org/web/packages/ggplot2/vignettes/ggplot2-specs.html).

| Library/Package                      | Piece of code                                   | What it's called       | What it does                                                                                                                                                                                      |
|--------------------------------------|-----------------------------------------------------|--------------------------|--------------------------------------------------------------------------------------------------------------------------------------|
| `ggplot2`         | [`ggplot()`](https://www.rdocumentation.org/packages/ggplot2/versions/3.1.0/topics/ggplot)                   | GG Plot               | Begins a plot that is finished by adding layers.        |
| `ggplot2`         | [`geom_boxplot()`](https://www.rdocumentation.org/packages/ggplot2/versions/3.1.0/topics/geom_boxplot)       | Boxplot               | Creates a boxplot when combined with ggplot()            |


### `readr` and `tibble`
Read the `readr` package documentation [**here**](https://www.rdocumentation.org/packages/readr/versions/1.3.0) and the package vignette [**here**](https://www.rdocumentation.org/packages/readr/versions/1.3.0/vignettes/readr.Rmd). <br>
Read the `tibble` package documentation [**here**](https://www.rdocumentation.org/packages/tibble/versions/1.4.2) and the package vignette [**here**](https://www.rdocumentation.org/packages/tibble/versions/1.4.2/vignettes/tibble.Rmd).

| Library/Package                      | Piece of code                                   | What it's called       | What it does                                                                                                                                                                                      |
|--------------------------------------|-----------------------------------------------------|---------------------------|-----------------------------------------------------------------------------------------------------------------------------------|
| `readr`           | [`read_tsv()`](https://www.rdocumentation.org/packages/readr/versions/1.3.0/topics/read_delim)           | Read TSV             | Reads in a TSV file from a specified file path. This function can be tailored to read in other common types of files. i.e. read_csv(), read_rds(), etc.                                          |
| `tibble`          | [`column_to_rownames()`](https://www.rdocumentation.org/packages/tibble/versions/1.4.2/topics/rownames) | Column to Rownames    | Transforms an existing column into the rownames.                            |                                                                                         