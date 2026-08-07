#!/usr/bin/env Rscript
#
# Create a normalized SpatialExperiment (SPE) object from Visium HD data
# input data is colon cancer data from 10x
#
# Steps:
#   1. Import Space Ranger segmented outputs as an SPE
#   2. Crop to a central region of the tissue
#   3. Calculate per-cell QC metrics
#   4. Identify and remove spatially local outliers
#   5. Library-size normalize and log transform
#   6. Write the resulting SPE to an RDS file

# Libraries --------------------------------------------------------------------

library(optparse)
library(SpatialExperiment)

# Constants --------------------------------------------------------------------

# Fraction of the tissue's full x and y extent to trim from each side when
# cropping, i.e. 1/3 keeps the central third of the tissue
crop_fraction <- 1 / 3

# Number of workers to use for SpotSweeper outlier detection
n_workers <- 4


# Command line options ---------------------------------------------------------

option_list <- list(
  make_option(
    "--input_dir",
    type = "character",
    help = "Path to the segmented outputs directory from Space Ranger"
  ),
  make_option(
    "--output_file",
    type = "character",
    help = "Path to write the normalized SPE"
  )
)

opts <- parse_args(OptionParser(option_list = option_list))

# Import data ------------------------------------------------------------------

spe <- VisiumIO::TENxVisiumHD(
  format = "mtx", # use mtx to ensure counts assay is an actual matrix and not just reference to h5 file
  images = "lowres",
  segmented_outputs = opts$input_dir
) |>
  VisiumIO::import()

# Pull image(s) fully into memory so the SPE no longer depends on the temporary
# directory used during import
imgs <- imgData(spe)
imgs$data <- imgs$data |>
  purrr::map(\(img) as(img, "LoadedSpatialImage"))
imgData(spe) <- imgs

spe$in_tissue <- TRUE

# Crop region of interest ------------------------------------------------------

# Zoom in on the center of the tissue to keep downstream steps tractable
xy <- spatialCoords(spe)
xs <- range(xy[, 1])
ys <- range(xy[, 2])

dx <- diff(xs) * crop_fraction
dy <- diff(ys) * crop_fraction

box <- list(
  xmin = xs[1] + dx,
  xmax = xs[2] - dx,
  ymin = ys[1] + dy,
  ymax = ys[2] - dy
)

sub <- spe[, xy[, 1] > box$xmin & xy[, 1] < box$xmax &
      xy[, 2] > box$ymin & xy[, 2] < box$ymax]

# Calculate QC metrics ---------------------------------------------------------

# Identify mitochondrial genes using gene symbols
mito_genes <- grepl("^MT-", rowData(sub)$Symbol)

# add qc metrics to coldata
sub <- scuttle::addPerCellQC(
  sub,
  subsets = list(mito = mito_genes)
)

# Detect spatial outliers ------------------------------------------------------

# Flag cells whose total UMI count is low relative to their local neighborhood
sub <- SpotSweeper::localOutliers(
  sub,
  metric = "sum",
  direction = "lower",
  log = TRUE,
  workers = n_workers
)

# Flag cells with few detected genes relative to their local neighborhood
sub <- SpotSweeper::localOutliers(
  sub,
  metric = "detected",
  direction = "lower",
  log = TRUE,
  workers = n_workers
)

# Flag cells with high mitochondrial percentage; not logged, since percentages
# are already on a bounded scale
sub <- SpotSweeper::localOutliers(
  sub,
  metric = "subsets_mito_percent",
  direction = "higher",
  log = FALSE,
  workers = n_workers
)

# Combine all outlier calls into a single column to filter on
sub$local_outliers <- sub$sum_outliers |
  sub$detected_outliers |
  sub$subsets_mito_percent_outliers

# Filter cells -----------------------------------------------------------------

# Remove local outliers, then any remaining cells with no counts
filtered_spe <- sub[, !sub$local_outliers]
filtered_spe <- filtered_spe[, filtered_spe$sum > 0]

# Normalize counts -------------------------------------------------------------

# Per-cell size factors from library size, then log transform
filtered_spe <- scuttle::computeLibraryFactors(filtered_spe)
filtered_spe <- scuttle::logNormCounts(filtered_spe)

# Export -----------------------------------------------------------------------

readr::write_rds(filtered_spe, opts$output_file)
