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
  dplyr::select(-which(apply(is.na(.), 2, all)), 
                sex = refinebio_sex) %>% 
  dplyr::mutate(
    # Clean up character names
    tissue = gsub(" ", "_", tolower(tissue)),
    # Turn our testing variables in factors
    #tissue = as.factor(tissue),
    #sex = as.factor(sex)
    )

## Put samples' data and metadata in the same order
# Make the data in the order of the metadata
gene_df <- gene_df %>% 
  dplyr::select(metadata$geo_accession)

#This set up is based on the example 11.4 Estrogen Data: A 2x2 Factorial Experiment in 
# https://bioc.ism.ac.jp/packages/2.9/bioc/vignettes/limma/inst/doc/usersguide.pdf
 
# Set up combined sex and tissue variable for model
sex_tissue <- as.factor(paste0(metadata$tissue, "_", metadata$sex))

# Set variable levels
sex_tissue = factor(sex_tissue,
                    labels = c("control_female", 
                               "control_male", 
                               "tumor_female", 
                               "tumor_male")
                    )

# Set up contrasts
contrasts(sex_tissue) <- cbind(sex = c(0, 1, 0, 1), 
                               tumor_male = c(0, 0, 0, 1),
                               tumor_female = c(0, 0, 1, 0))

# Create the design matrix from the genotype information
des_mat <- model.matrix(~ sex_tissue)

# Neaten up column names
colnames(des_mat) <- c("intercept", "sex", "tumor_male", "tumor_female")

# Apply linear model to data
fit <- limma::lmFit(gene_df, design = des_mat)

# Set up a contrast matrix that focuses on tumors
contrast_mat <- cbind(tumor_male = c(0, 0, 0, 1),
                      tumor_female = c(0, 0, 1, 0))

# Do final fitting
fit2 <- contrasts.fit(fit, contrast_mat)

# Apply empirical Bayes to smooth standard errors
fit2 <- eBayes(fit2)

# Apply multiple testing correction and obtain stats
stats_df <- limma::topTable(fit2, number = nrow(gene_df)) %>% 
  # Move ensembl IDs to their own column
  tibble::rownames_to_column("ensembl_id") %>%
  # Map ensembl IDs to their associated gene symbols
  dplyr::mutate("gene_symbol" = AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db, 
    keys = ensembl_id, 
    column = "SYMBOL", 
    keytype = "ENSEMBL")) %>%  
  # Clean up column names and reorder 
  dplyr::select(
    ensembl_id, 
    gene_symbol,
    tumor_female, 
    tumor_male, 
    avg_expression = AveExpr,
    #log_fold_change = logFC, 
    #t_value = t, 
    p_value = P.Value, 
    adj_p_value = adj.P.Val, 
    #log_odds = B
    ) %>% 
  # Write this to TSV
  readr::write_tsv(file.path(data_dir, "gene_results_GSE44971.tsv"))

colnames(stats_df)
