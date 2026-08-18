#!/usr/bin/env Rscript

# This script calculates Moran's I with Monte Carlo P values
# output is saved to a tsv file with the gene ids and moran results

suppressPackageStartupMessages({
  library(optparse)
  library(SpatialExperiment)
})

option_list <- list(
  make_option(
    opt = c("-i", "--input"),
    type = "character",
    help = "Path to the input osteo rds SPE file"
  ),
  make_option(
    opt = c("-o", "--output"),
    type = "character",
    help = "Path output table"
  ),
  make_option(
    opt = c("-r", "--reps"),
    type = "integer",
    default = 100,
    help = "Number of Monte Carlo replicates"
  ),
  make_option(
    opt = c("-t", "--threads"),
    type = "integer",
    default = 8,
    help = "Number of threads to use for parallelization"
  )
)
# parse arguments
opt <- parse_args(OptionParser(option_list = option_list))

# read file
spe <- readr::read_rds(opt$input)


gene_variance <- scran::modelGeneVar(spe)

# calculate the top 2000 HVGs
hvg_vector <- scran::getTopHVGs(gene_variance, n = 2000)

hvg_spe <- spe[hvg_vector, ]

# Convert to SFE

sfe <- SpatialFeatureExperiment::toSpatialFeatureExperiment(hvg_spe)

# use the Visium barcodes to define the graph
SpatialFeatureExperiment::colGraph(
  sfe
) <- SpatialFeatureExperiment::findVisiumGraph(sfe)

# Calculate Moran's I with MC stats

sfe <- Voyager::runUnivariate(
  sfe,
  type = "moran.mc",
  nsim = opt$reps,
  BPPARAM = BiocParallel::MulticoreParam(opt$threads)
)

moran_mc_result <- rowData(sfe) |>
  as.data.frame() |>
  dplyr::select(
    ID,
    dplyr::starts_with("moran.mc_")
  )

readr::write_tsv(moran_mc_result, opt$output)
