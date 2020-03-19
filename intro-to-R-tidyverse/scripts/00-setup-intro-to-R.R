# Set up statistics table for intro to R module

# C. Savonen 

# 2020

`%>%` <- dplyr::`%>%`

if (!("limma" %in% installed.packages())) {
  # Install limma
  BiocManager::install("limma", update = FALSE)
}

input_dir <- "data"
output_dir <- "results"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

plots_dir <- file.path(output_dir, "plots")

if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
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
stats_df <- limma::topTable(fit, number = nrow(gene_df)) 
%>% 
  readr::write_tsv(file.path(output_dir, "gene_results_GSE44971.tsv"))

