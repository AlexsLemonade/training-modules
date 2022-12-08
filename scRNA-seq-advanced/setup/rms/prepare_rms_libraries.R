# Script to prepare integration data for use in the `scRNA-seq-advanced` workshop
# Libraries are from SCPCP000005, and all `_filtered.rds` files are expected to
#  be present in the `input_sce_path` argument.

## Setup -----------------

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
    "--input_sce_file",
    help = "Path to input _filtered.rds file."
  ),
  make_option(
    "--output_sce_file",
    help = "Output file for processed SCE with celltypes"
  ),
  make_option(
    "--celltypes_file",
    help = "Path to `celltypes.tsv`"
  )
)
opts <- parse_args(OptionParser(option_list = options))

# Make output directory if it doesn't exist
output_dir <- dirname(opts$output_sce_file)
if (!(dir.exists(output_dir))) {
  dir.create(output_dir)
}

# Read SCE
raw_sce <- readr::read_rds(opts$input_sce_file)

# Define library_name for later use parsing cell types
library_name <- stringr::str_extract(opts$input_sce_file, "SCPCL\\d+")

# Read celltypes TSV and create new columns for "broad" and "fine" cell types
celltypes_df <- readr::read_tsv(opts$celltypes_file) %>%
  dplyr::mutate(celltype_broad = stringr::str_remove(celltype, "-[ABCD]$")) %>%
  # We'll also rename `celltype` to `celltype_fine`
  dplyr::rename(celltype_fine = celltype)


### Filter and normalize the SCE ------------------

# Filter based on miQC:
filtered_sce <- raw_sce[, which(raw_sce$miQC_pass)]

# Filter to detected >= 500
# Note that `sum` filtered was handled in scpca-nf
filtered_sce <- filtered_sce[, filtered_sce$detected >= 500]

# Normalize:
qclust <- scran::quickCluster(filtered_sce)
filtered_sce <- scran::computeSumFactors(filtered_sce, clusters = qclust)
norm_sce <- scater::logNormCounts(filtered_sce)


### Dimension reductions ----------------------------------

# calculate high-variance genes
gene_variance <- scran::modelGeneVar(norm_sce)
hvg_list <- scran::getTopHVGs(gene_variance,
                              n = 2000)

# Add PCA and UMAP into the SCE
norm_sce <- norm_sce %>%
  scater::runPCA(subset_row = hvg_list) %>%
  scater::runUMAP(dimred = "PCA")


### Add cell type information to SCE ----------------------


# Add a matching barcode variable into the colData
colData(norm_sce)$barcode <- colnames(norm_sce)

# filter to cells in this library
celltypes_df <- celltypes_df %>%
  dplyr::filter(library == library_name) %>%
  dplyr::select(-library)

# join colDdata to celltype
celltype_coldata <-  colData(norm_sce) %>%
  as.data.frame() %>%
  # Use a left_join since there may be more cells in the SCE than there are in
  #  the celltypes
  dplyr::left_join(celltypes_df, by = "barcode") %>%
  #convert back to a DataFrame with row names
  S4Vectors::DataFrame(row.names = .$barcode)

# check that the names of the cells are still consistent
if (all(rownames(celltype_coldata) == colnames(norm_sce))){
  colData(norm_sce) <- celltype_coldata
} else {
  stop("Error joining celltype_df to colData")
}


### Save -------------------

# Write the SCE file
readr::write_rds(norm_sce, opts$output_sce_file)



### sessionInfo to log -------------
sessionInfo()










