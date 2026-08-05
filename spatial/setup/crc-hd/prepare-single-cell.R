#!/usr/bin/env Rscript
#
# Create a normalized SingleCellExperiment (SCE) object to use as reference for
# annotating the Visium HD colon cancer dataset

# Steps:
#   1. Import Cell Ranger output as an SCE
#   2. Import QC and cell type annotations
#   3. Remove low quality cells
#   4. Normalize cells
#   5. Downsample to 1000 cells per cell type
#   6. Write the resulting SCE to an RDS file

# Libraries --------------------------------------------------------------------

library(optparse)
library(SingleCellExperiment)

set.seed(2026)

# Command line options ---------------------------------------------------------

option_list <- list(
  make_option(
    "--input_h5_file",
    type = "character",
    help = "Path to the filtered feature bc h5ad output file from Cell Ranger"
  ),
  make_option(
    "--output_file",
    type = "character",
    help = "Path to write the normalized SCE"
  )
)

opts <- parse_args(OptionParser(option_list = option_list))

# Read in data -----------------------------------------------------------------
# url to file on Github with reference cell type annotations
annotations_url <- "https://raw.githubusercontent.com/10XGenomics/HumanColonCancer_VisiumHD/main/MetaData/SingleCell_MetaData.csv.gz"

# read in sce and annotations
sce <- DropletUtils::read10xCounts(opts$input_dir)
annotations_df <- readr::read_csv(annotations_url)

# Prep SCE ---------------------------------------------------------------------
# add in cell type information and QCFilter column to colData
colData(sce) <- colData(sce) |>
  as.data.frame() |>
  dplyr::left_join(annotations_df, by = c("Barcode")) |>
  DataFrame(row.names = colnames(sce))


# remove bad cells based on QCFilter column
sce <- sce[, sce$QCFilter == "Keep"]

# normalize on the full dataset
qclust <- scran::quickCluster(sce)
sce <- scran::computeSumFactors(sce, clusters = qclust)
sce <- scater::logNormCounts(sce)

# Downsample -------------------------------------------------------------------

# this object has 260K cells in it so we will downsample to 1000 cells per cell type
# this makes running SingleR much faster
downsampled_cells <- colData(sce) |>
  as.data.frame() |>
  dplyr::group_by(Level1) |>
  dplyr::slice_sample(n = 1000) |>
  dplyr::pull(Barcode)

sce <- sce[ ,downsampled_cells]

# Export -----------------------------------------------------------------------

# export the downsampled and normalized object
readr::write_rds(sce, opts$output_file)
