#!/usr/bin/env Rscript
#
# Make live versions of .Rmd files in training modules
#
# Replaces code in chunks with a chunk option of `live = TRUE`
# Comments are preserved
#
# If --rendering option is FALSE, the rendering is skipped.
# Default is TRUE -- rendering is not skipped. 

# Load library:
library(optparse)

# Set up optparse options
option_list <- list(
  make_option(
    opt_str = "--rendering", action = "character",
    default = "TRUE", help = "Needs a 'TRUE/FALSE' to determine whether the markdown::render() steps will be run for all notebooks."
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Turn character into logical
opt$rendering <- as.logical(opt$rendering)

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
                         "05-dimension_reduction_scRNA-seq.Rmd")),
             file.path(root_dir,  "machine-learning",
                       c("01-openpbta_heatmap.Rmd",
                         "02-openpbta_consensus_clustering.Rmd",
                         "03-openpbta_PLIER.Rmd",
                         "04-openpbta_plot_LV.Rmd")),
             file.path(root_dir, "pathway-analysis",
                       c("01-overrepresentation_analysis.Rmd",
                         "02-gene_set_enrichment_analysis.Rmd",
                         "03-gene_set_variation_analysis.Rmd"))
             )


# Rerender notebooks if --skiprendering is FALSE
if (!opt$skiprendering) {
  purrr::map(infiles, rmarkdown::render, envir = new.env(), quiet = TRUE)
}

# new files will be made with -live.Rmd suffix
outfiles <- stringr::str_replace(infiles, "(.*)\\.Rmd$", "\\1-live.Rmd")

# Generate live versions
# capture to avoid printing to stdout
out <- purrr::map2(infiles, outfiles, exrcise::exrcise, replace_flags = "live")
