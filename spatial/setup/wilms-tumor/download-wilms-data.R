#!/usr/bin/env Rscript

library(ScPCAr)
library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--email"),
    type = "character",
    help = "Email to use for authenticating with ScPCA Portal"
  ),
  make_option(
    opt_str = c("--outdir"),
    type = "character",
    default = "../../data/wilms-tumor/SCPCS000190",
    help = "Path where spaceranger directory should live"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))


stopifnot(
  "Must provide an email with --email" = !is.null(opt$email),
  "Must provide an outdir with --outdir" = !is.null(opt$outdir)
)
fs::dir_create(opt$outdir)


# Get an authentication token for use with the ScPCA Portal
auth_token <- get_auth(email = opt$email, agree = TRUE)


# Download a data for a sample
file_paths <- download_sample(
  sample_id = "SCPCS000190",
  auth_token = auth_token,
  destination = "scpca_data",
  format = "spatial",
  overwrite = TRUE
)

# capture output directory name
downloaded_dir <- file.path(
  "scpca_data",
  dir("scpca_data")
)


# Move & rename what we're keeping
fs::file_move(
  file.path(downloaded_dir, "SCPCL000429_spatial"),
  file.path(opt$outdir, "spaceranger")
)

# Delete what we don't want - keep the filtered output and summary report
fs::file_delete(
  file.path(opt$outdir, "spaceranger", "SCPCL000429_metadata.json")
)

# Clean up
fs::dir_delete("scpca_data")
