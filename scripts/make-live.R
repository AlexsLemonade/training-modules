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
    opt_str = "--notebook", type = "character",
    help = "The notebook file to process."),
  make_option(
    opt_str = "--render", type = "character",
    default = "TRUE", help = "Needs a 'TRUE/FALSE' to determine whether the markdown::render() step will be run."
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

if (! tolower(opt$render) %in% c("true", "false", "t", "f")){
  stop("--render option must be TRUE or FALSE")
}

# Turn character into logical
render <- as.logical(opt$render)

# Install exrcise package if needed.
if (!"exrcise" %in% installed.packages()){
  remotes::install_github("AlexsLemonade/exrcise", dependencies = TRUE, upgrade = "never")
}


message("Processing ", opt$notebook)
# Rerender notebooks if --rendering is TRUE
if (render) {
  rmarkdown::render(opt$notebook, env = new.env(), quiet = TRUE)
}

# new files will be made with -live.Rmd suffix
outfile <- stringr::str_replace(opt$notebook, "(.*)\\.Rmd$", "\\1-live.Rmd")

# Generate live versions
# capture to avoid printing to stdout
exrcise::exrcise(opt$notebook, outfile, replace_flags = "live")
