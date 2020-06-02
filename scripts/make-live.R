#!/usr/bin/env Rscript
#
# Make live versions of .Rmd files in training modules

# Install exrcise package if needed.
if (!"exrcise" %in% installed.packages()){
  remotes::install_github("AlexsLemonade/exrcise", dependencies = TRUE)
}

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

infiles <- c(file.path(root_dir, "intro-to-R-tidyverse", 
                       c("01-intro_to_base_R.Rmd", 
                         "02-intro_to_ggplot2.Rmd",
                         "03-intro_to_tidyverse.Rmd"))
)

outfiles <- stringr::str_replace(infiles, "(.*)\\.Rmd$", "\\1-live.Rmd")

purrr::map2(infiles, outfiles, exrcise::exrcise, replace_flags = "live")