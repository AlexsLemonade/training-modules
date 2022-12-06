## ----setup---------------------------------------------------------------------------------
# Load libraries
library(optparse)
library(magrittr)
library(SingleCellExperiment)


# Setting the seed for reproducibility
set.seed(12345)

# parse command line options
options <- list(
  make_option(
    "--raw_matrix_dir",
    help = "Cell Ranger raw data matrix directory"
  ),
  make_option(
    "--mito_file",
    help = "tsv file with mitochondrial genes"
  ),
  make_option(
    "--outfile",
    help = "Output file path for normalized SCE file"
  )
)
opts <- parse_args(OptionParser(option_list = options))

# Read SCE
raw_sce <- DropletUtils::read10xCounts(opts$raw_matrix_dir) 

# Split out the ADT into an altExp
raw_sce <- splitAltExps(raw_sce, 
                        # "Gene Expression" or "Antibody Capture"
                        f = rowData(raw_sce)$Type, 
                        # The Type which should be the main experiment
                        ref = "Gene Expression")


# name the alt experiment
altExpNames(raw_sce) <- "ADT"

# run emptyDropsCellRanger
droplet_df <- DropletUtils::emptyDropsCellRanger(
  counts(raw_sce), 
  BPPARAM = BiocParallel::MulticoreParam(4)
)

cells_to_retain <- which(droplet_df$FDR <= 0.01)
filtered_sce <- raw_sce[, cells_to_retain]


# read in mitochondrial gene files
mito_genes <- readr::read_tsv(opts$mito_file) %>%
  dplyr::filter(gene_id %in% rownames(filtered_sce)) %>%
  dplyr::pull(gene_id)


# calculate cell stats
filtered_sce <- scuttle::addPerCellQC(filtered_sce,
                                      subsets = list(mito = mito_genes))

message("filtering with miQC")

# calculate miQC model and filter 
miqc_model <- miQC::mixtureModel(filtered_sce) 
filtered_sce <- miQC::filterCells(filtered_sce, model = miqc_model)
message("filtering on ADTs")
# generate ambient model from empty droplets
adt_ambient <- DropletUtils::ambientProfileEmpty(altExp(raw_sce))

# Run QC on the tags
adt_qc_df <- DropletUtils::cleanTagCounts(
  altExp(filtered_sce),
  ambient = adt_ambient
)

message(paste("removing", sum(adt_qc_df$discard), "cells due to ADT filtering"))
# remove discard-flagged cells
filtered_sce <- filtered_sce[, which(!adt_qc_df$discard)]

# Normalize RNA
message("normalizing RNA")
qclust <- scran::quickCluster(filtered_sce)
filtered_sce <- scran::computeSumFactors(filtered_sce, clusters = qclust)
normalized_sce <- scuttle::logNormCounts(filtered_sce)


# Normalize ADT
message("Normalizing ADT counts")
adt_sf <- scuttle::medianSizeFactors(
  altExp(normalized_sce), 
  reference = adt_ambient )

altExp(normalized_sce) <- scater::logNormCounts(
  altExp(normalized_sce), 
  size.factors = adt_sf
)


# Identify to 2000 HVGs by modelling
message("Identifying HVGs")
num_genes <- 2000

gene_variance <- scran::modelGeneVar(normalized_sce)
hv_genes <- scran::getTopHVGs(gene_variance,
                              n = num_genes)

# add PCA & UMAP
message("running dimensionality reduction")
normalized_sce <- normalized_sce %>% 
  scater::runPCA(subset_row = hv_genes) %>% 
  scater::runUMAP(dimred = "PCA")


# save output file
saveRDS(normalized_sce, file = opts$outfile)


## ----session info--------------------------------------------------------------------------
sessionInfo()

