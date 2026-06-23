#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(optparse)
})


option_list <- list(
  make_option(
    opt_str = c("--counts_file"),
    type = "character",
    default = "GSE189955_RNA_Count_matrix.csv.gz",
    help = "Counts CSV file downloaded from GEO"
  ),
  make_option(
    opt_str = c("--metadata_file"),
    type = "character",
    default = "GSE189955_RNA_meta_data.csv.gz",
    help = "Metadata (with cell type annotations) CSV file downloaded from GEO"
  ),
  make_option(
    opt_str = c("--features_file"),
    type = "character",
    default = "../../data/ovarian-carcinoma/spaceranger/filtered_feature_bc_matrix/features.tsv.gz",
    help = "Path to the features.tsv output associated with the eventual query object. Used to support gene symbol -> ensembl conversion"
  ),
  make_option(
    opt_str = c("--output_file"),
    type = "character",
    default = "../../data/reference/GSE189955_ovarian_ref.rds",
    help = "Path to output RDS with prepared reference"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(
  "The provided counts_file does not exist" = file.exists(opt$counts_file),
  "The provided metadata_file does not exist" = file.exists(opt$metadata_file),
  "The provided features_file does not exist" = file.exists(opt$features_file),
  "output_file was not provided" = !is.null(opt$output_file)
)

# This setting is needed to read in the super wide CSV from GEO
# 131072 is default, and 10 times that seems to be enough
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10)


# First, we'll prepare the counts matrix -------------
# Convert from gene symbol -> ensembl, using the features tsv for ID reference
# worth noting the colnames aren't actually barcodes but they are unique ids!
raw_matrix <- readr::read_csv(opt$counts_file)
full_matrix <- as.matrix(raw_matrix[, -1])
rownames(full_matrix) <- raw_matrix$...1

features_df <- readr::read_tsv(
  opt$features_file,
  col_names = c("ensembl_id", "gene_symbol", "type")
)

# match gene symbols in GEO matrix to features_df, get their ensembl IDs
keep_features <- intersect(rownames(full_matrix), features_df$gene_symbol)

# subset, and convert to ensembl
subset_matrix <- full_matrix[keep_features, ]

# convert rownames from symbol to ensembl
symbol_to_ensembl <- setNames(features_df$ensembl_id, features_df$gene_symbol)
rownames(subset_matrix) <- symbol_to_ensembl[rownames(subset_matrix)]

# Now, we can create the SCE
# the output dimensions are 17201 5329
sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = subset_matrix)
)

# Second, add annotations to SCE  ------------------------
annotations_df <- readr::read_csv("GSE189955_RNA_meta_data.csv.gz") |>
  tibble::column_to_rownames("...1")

# reorder to match SCE column order
annotations_df <- annotations_df[colnames(sce), ]

# Add it as colData; celltypes are in `Cell_Type`
colData(sce) <- DataFrame(annotations_df)

# Finally, remove cells where the annotation is NA
sce <- sce[, !is.na(sce$Cell_Type)]

# Export --------
readr::write_rds(sce, opt$output_file, compress = "gz")
