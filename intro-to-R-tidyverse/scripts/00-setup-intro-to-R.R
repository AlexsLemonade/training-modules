# Set up gene level statistics table for intro to R module

# C. Savonen for ALSF - CCDL

# 2020

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Establish R dir
intro_r_dir <- file.path(root_dir, "intro-to-R-tidyverse")

# Establish input and output directories
data_dir <- file.path(intro_r_dir, "data")
plots_dir <- file.path(output_dir, "plots")

# Create these output directories if they don't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Install limma if not installed
if (!("limma" %in% installed.packages())) {
  # Install limma
  BiocManager::install("limma", update = FALSE)
}

# Install the Human annotation package
if (!("org.Hs.eg.db" %in% installed.packages())) {
  # Install org.Hs.eg.db
  BiocManager::install("org.Hs.eg.db", update = FALSE)
}

# Read in the dataset
gene_df <- readr::read_tsv(file.path(input_dir, 
                                     "GSE44971.tsv")) %>% 
  tibble::column_to_rownames("Gene")

# Read in metadata
metadata <- readr::read_tsv(file.path(input_dir, 
                                      "metadata_GSE44971.tsv")) %>%
  dplyr::select(-which(apply(is.na(.), 2, all))) %>% 
  dplyr::mutate(tissue = as.factor(tissue))

## Put samples' data and metadata in the same order
# Make the data in the order of the metadata
gene_df <- gene_df %>% 
  dplyr::select(metadata$geo_accession)

# Create the design matrix from the genotype information
des_mat <- model.matrix(~metadata$tissue)

# Apply linear model to data
fit <- limma::lmFit(gene_df, design = des_mat)

# Apply empirical Bayes to smooth standard errors
fit <- limma::eBayes(fit)

# Apply multiple testing correction and obtain stats
stats_df <- limma::topTable(fit, number = nrow(gene_df)) %>% 
  # Move ensembl IDs to their own column
  tibble::rownames_to_column("ensembl_ids") %>%
  # Map ensembl IDs to their associated gene symbols
  dplyr::mutate("gene_symbols" = AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db, 
    keys = ensembl_ids, 
    column = "SYMBOL", 
    keytype = "ENSEMBL")) %>%  
  # Clean up column names and reorder 
  dplyr::select(
    ensembl_ids, 
    gene_symbols,
    avg_expression = AveExpr,
    log_fold_change = logFC, 
    t_value = t, 
    p_value = P.Value, 
    adj_p_value = adj.P.Val, 
    log_odds = B
    ) %>% 
  # Write this to TSV
  readr::write_tsv(file.path(data_dir, "gene_results_GSE44971.tsv"))

