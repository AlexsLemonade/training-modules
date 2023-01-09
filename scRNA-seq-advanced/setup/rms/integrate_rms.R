# This script reads in a list of processed SCE files, merges into one SCE object,
# and then performs integration using fastMNN.
# The returned object contains the merged object with the corrected expression and
# corrected PCA and UMAP results.

# Load libraries
suppressPackageStartupMessages({
  library(optparse)
  library(magrittr)
  library(SingleCellExperiment)
})

# Set the seed for reproducibility
set.seed(12345)

# Parse command line options
options <- list(
  make_option(
    "--input_file_list",
    help = "Quoted space separated list of processed SCE objects to integrate."
  ),
  make_option(
    "--integrated_sce_file",
    help = "Output file for integrated SCE object"
  )
)
opt <- parse_args(OptionParser(option_list = options))

# Setup ------------------------------------------------------------------------

# create list of SCE objects 
if(is.null(opt$input_file_list)){
  stop("List of input files containing individual SCE objects to merge is missing.")
} else {
  sce_paths <- unlist(stringr::str_split(opt$input_file_list, ' '))
} 

# check that at least 2 SCE files have been provided 
if(length(sce_paths) < 2){
  stop("Must provide at least 2 SCE files to integrate.")
}

# extract sample ids from file names 
sample_ids <- stringr::str_remove(basename(sce_paths), ".rds")

# read in list of SCE objects 
sce_list <- purrr::map(
  sce_paths,
  readr::read_rds
)
# name list with sample ids 
names(sce_list) <- sample_ids

# remove SCPCL000482
sce_list <- sce_list[which(names(sce_list) != "SCPCL000482")]


# Make output directory if it doesn't exist
output_dir <- dirname(opt$integrated_sce_file)
if (!(dir.exists(output_dir))) {
  dir.create(output_dir)
}

# Merge SCE objects ------------------------------------------------------------

# set up function to re-format individual SCEs prior to merging
format_sce_list <- function(sce, sce_name, shared_genes) {
  
  # Reduce to shared genes
  sce <- sce[shared_genes,]
  
  # Update the colData row names so we can know where cells came from back
  rownames(colData(sce)) <- glue::glue("{sce_name}-{rownames(colData(sce))}")
  
  # Update the rowData names, except `gene_symbol`
  names(rowData(sce)) <- c("gene_symbol", 
                           glue::glue("{sce_name}-mean"),
                           glue::glue("{sce_name}-detected"))
                                
  # Add in sample-level information as stand-alone column
  colData(sce)$sample <- sce_name
  
  # Return
  return(sce)
}

# grab shared genes 
shared_genes <- sce_list %>%
  # get rownames for each entry in sce_list
  purrr::map(rownames) %>%
  # find their intersection
  purrr::reduce(intersect)


# Prepare the sce_list for combining
sce_list <- purrr::imap(sce_list,
                        format_sce_list,
                        shared_genes = shared_genes)

# combine into one sce
combined_sce <- do.call(cbind, sce_list)

# Dimension Reduction ----------------------------------------------------------

# select combined highly variable genes
gene_variance <- scran::modelGeneVar(combined_sce,
                                     # For combined SCEs, specify the group column groups:
                                     block = combined_sce$sample)
hvg_list <- scran::getTopHVGs(gene_variance,
                              n = 2000)

# store HVG in combined object for future use
metadata(combined_sce)$combined_hvg <- hvg_list

# Add PCA and UMAP into the SCE
combined_sce <- combined_sce %>%
  scater::runPCA(subset_row = hvg_list) %>%
  scater::runUMAP(dimred = "PCA")

# Integration ------------------------------------------------------------------

# perform fastMNN integration 
integrated_sce <- batchelor::fastMNN(combined_sce,
                                     batch = combined_sce$sample,
                                     subset.row = hvg_list,
                                     correct.all = TRUE)

# add reconstructed object back to combined sce
assay(combined_sce, "fastmnn_corrected") <- assay(integrated_sce, "reconstructed")

# add integrated PCs to combined sce 
reducedDim(combined_sce, "fastmnn_PCA") <- reducedDim(integrated_sce, "corrected")

# add integrated UMAP to combined SCE 
combined_sce <- scater::runUMAP(combined_sce,
                                dimred = "fastmnn_PCA",
                                name = "fastmnn_UMAP")

# write out combined SCE file with added integration PCs
readr::write_rds(combined_sce, opt$integrated_sce_file, compress = "gz")
