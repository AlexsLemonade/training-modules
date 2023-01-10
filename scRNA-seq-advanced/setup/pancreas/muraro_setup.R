#!/usr/bin/env RScript

# Script to prepare sample files for integration exercise


# Setup -----------------

suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
})

set.seed(2023)

# Parse command line options
options <- list(
  make_option(
    "--outdir",
    help = "Path to output directory for split data files."
  )
)
opts <- parse_args(OptionParser(option_list = options))

# Make output directory if it doesn't exist
if (!(dir.exists(opts$outdir))) {
  dir.create(opts$outdir)
}


# Download and process pancreas data -----------

# install scRNAseq package if needed
if(!requireNamespace("scRNAseq", quietly = TRUE)){
  BiocManager::install("scRNAseq", update = FALSE)
}

# get muraro data
muraro_sce <- scRNAseq::MuraroPancreasData(ensembl = TRUE)

# add ensembl column to rowData & remove originalName
rowData(muraro_sce)$ensembl <- rownames(muraro_sce)
rowData(muraro_sce)$originalName <- NULL

# filtering based on http://bioconductor.org/books/3.16/OSCA.workflows/muraro-human-pancreas-cel-seq.html

stats <- scuttle::perCellQCMetrics(muraro_sce)
qc <- scuttle::quickPerCellQC(stats, 
                              # use ERCC spike-in as QC measure
                              sub.fields = "altexps_ERCC_percent",
                              batch = muraro_sce$donor, 
                              # calculate cutoff stats without outlier sample
                              subset = muraro_sce$donor!="D28")
muraro_sce <- muraro_sce[, !qc$discard]


# Split merged data by donor ------

# get donor list
muraro_sce_list <- unique(muraro_sce$donor) |>
  purrr::set_names() |>
  # select cells for each donor
  purrr::map( \(id) muraro_sce[, muraro_sce$donor == id] )

# Remove the donor-specific cell labels
muraro_sce_list <- muraro_sce_list |>
  purrr::map( \(sce){
    colnames(sce) <- stringr::str_sub(colnames(sce), 5)
    sce
  })



# Normalization and dimension reduction ----

muraro_sce_list <- muraro_sce_list |>
  purrr::map(\(sce){
    
    clusters <- scran::quickCluster(sce)
    sce <- sce |>
      scran::computeSumFactors(clusters = clusters) |>
      scuttle::logNormCounts()
    
    # use spike-in with plate blocking to model gene variance
    var_model <- scran::modelGeneVarWithSpikes(sce, "ERCC", block=sce$plate)
    hvgs <- scran::getTopHVGs(var_model, prop=0.1)
    
    # dimensionality reduction
    sce <- sce |>
      scater::runPCA(subset_row = hvgs) |>
      scater::runUMAP(dimred = "PCA")
    sce
  })

# Remove donor column and spike-in experiment. 

muraro_sce_list <- muraro_sce_list |>
  purrr::map(\(sce){
    sce$donor <- NULL
    altExp(sce) <- NULL
    # add QC fields
    sce <- scuttle::addPerCellQC(sce)
    sce
  })

# Write output files ------
muraro_sce_list |> 
  purrr::iwalk(\(sce, id){
    readr::write_rds(sce, 
                     file.path(opts$outdir, glue::glue("{id}_sce.rds")))
  })

# log session info
sessionInfo()

