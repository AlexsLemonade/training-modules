#!/usr/bin/env Rscript
#
# Make live versions of .Rmd files in training modules
#
# Replaces code in chunks with a chunk option of `live = TRUE`
# Comments are preserved

# Install exrcise package if needed.
if (!"exrcise" %in% installed.packages()){
  remotes::install_github("AlexsLemonade/exrcise", dependencies = TRUE, upgrade = "never")
}

# find the project root
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# list of files to transform
infiles <- c(file.path(root_dir, "intro-to-R-tidyverse",
                       c("01-intro_to_base_R.Rmd",
                         "02-intro_to_ggplot2.Rmd",
                         "03-intro_to_tidyverse.Rmd")),
             file.path(root_dir, "RNA-seq",
                       c("02-gastric_cancer_tximport.Rmd",
                         "03-gastric_cancer_exploratory.Rmd",
                         "05-nb_cell_line_DESeq2.Rmd")),
             file.path(root_dir, "scRNA-seq",
                       c("01-filtering_scRNA-seq.Rmd",
                         "02-normalizing_scRNA-seq.Rmd",
                         "04-tag-based_scRNA-seq_processing.Rmd",
                         "05-dimension_reduction_scRNA-seq.Rmd"))
             )


# Rerender notebooks
purrr::map(infiles, rmarkdown::render, envir = new.env(), quiet = TRUE)


# new files will be made with -live.Rmd suffix
outfiles <- stringr::str_replace(infiles, "(.*)\\.Rmd$", "\\1-live.Rmd")

# Generate live versions
# capture to avoid printing to stdout
out <- purrr::map2(infiles, outfiles, exrcise::exrcise, replace_flags = "live")
