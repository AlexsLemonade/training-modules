#!/usr/bin/env Rscript

# This script reformats the osteo reference:
# - downsamples to <= 10000 cells per L1 cell type
# - converts rownames to Ensembl, using a provided SPE
# - converts to SCE

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SpatialExperiment)
})


option_list <- list(
  make_option(
    opt = c("-i", "--input"),
    type = "character",
    help = "Path to the input osteo reference RDS file."
  ),
  make_option(
    opt = c("-o", "--output"),
    type = "character",
    help = "Path to the output formatted osteo reference RDS file."
  ),
  make_option(
    opt = c("--spe"),
    type = "character",
    help = "Path to SPE file to use for mapping IDs"
  )
)
# parse arguments
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(
  "Must provide the input osteo ref file with -i" = file.exists(opt$input),
  "Must provide the output osteo ref file with -o" = !is.null(opt$output),
  "Must provide an SPE to use for mapping IDs with --spe" = file.exists(opt$spe)
)

# Read in SPE
spe <- readr::read_rds(opt$spe)

## Read in the original reference object with qs
mm_mets <- qs::qread(opt$input)

# Downsample for reproducibility
cells_keep <- mm_mets@meta.data |>
  tibble::rownames_to_column("barcode") |>
  dplyr::group_by(Ann_Level1) |>
  dplyr::slice_sample(n = 10000) |>
  dplyr::pull(barcode)

mm_mets <- mm_mets[, cells_keep]

# Convert to SCE
ref_sce <- Seurat::as.SingleCellExperiment(mm_mets)

# Convert row names to Ensembl as close as we can
map_df <- rowData(spe) |>
  as.data.frame() |>
  dplyr::select(gene_symbol = Symbol, ensembl  = ID)

# Match gene symbols in reference to the map_df
match_idx <- match(rownames(ref_sce), map_df$gene_symbol)

# Track unmatched genes for posterity
unmatched <- rownames(ref_sce)[is.na(match_idx)]
cat(sprintf("Genes lost (no Ensembl match): %d / %d (%.1f%%)\n",
            length(unmatched), nrow(ref_sce),
            100 * length(unmatched) / nrow(ref_sce)))
# ---> Genes lost (no Ensembl match): 3077 / 16807 (18.3%)

# Subset reference to matched genes only and toggle in ensembl rownames
matched_mask <- !is.na(match_idx)
ref_sce_ensembl <- ref_sce[matched_mask, ]
rownames(ref_sce_ensembl) <- map_df$ensembl[match_idx[matched_mask]]

# Remove extra items in the object to save a bit of space
reducedDims(ref_sce_ensembl) <- NULL
logcounts(ref_sce_ensembl) <- NULL

# Export
readr::write_rds(ref_sce_ensembl, opt$output)
