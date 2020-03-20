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
gene_df <- readr::read_tsv(file.path(data_dir, 
                                     "GSE44971.tsv")) %>% 
  tibble::column_to_rownames("Gene")

# Read in metadata
metadata <- readr::read_tsv(file.path(data_dir, 
                                      "metadata_GSE44971.tsv")) %>%
  dplyr::select(-which(apply(is.na(.), 2, all)), 
                sex = refinebio_sex) %>% 
  dplyr::mutate(
    # Clean up character names
    tissue = gsub(" ", "_", tolower(tissue))
    )

## Put samples' data and metadata in the same order
# Make the data in the order of the metadata
gene_df <- gene_df %>% 
  dplyr::select(metadata$geo_accession)

# This set up is based on the sum to zero parametrization example from 9.5.4 Classic Interaction Models in 
# https://bioconductor.org/packages/3.1/bioc/vignettes/limma/inst/doc/usersguide.pdf

# Set up combined sex and tissue variable for model
sex <- factor(metadata$sex)
tissue <- factor(metadata$tissue)
design <- model.matrix(~sex*tissue)

# Neaten up column names
colnames(design) <- c("intercept", "male_female", "astrocytoma_normal", "interaction")

# Apply linear model to data
fit <- limma::lmFit(gene_df, design)

# Apply empirical Bayes to smooth standard errors
fit2 <- limma::eBayes(fit)

# collect results for each contrast
stats_df <- dplyr::bind_rows(
  male_female = limma::topTable(fit2, coef = "male_female", number = Inf, sort = "none" ) %>% 
    tibble::rownames_to_column("ensembl_id"),
  astrocytoma_normal = limma::topTable(fit2, coef = "astrocytoma_normal", number = Inf, sort = "none")%>% 
    tibble::rownames_to_column("ensembl_id"),
  interaction = limma::topTable(fit2, coef = "interaction", number = Inf, sort = "none")%>% 
    tibble::rownames_to_column("ensembl_id"),
  .id = "contrast"
)

# Apply multiple testing correction and obtain stats
stats_df <- stats_df %>% 
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
    contrast,
    avg_expression = AveExpr,
    t,
    p_value = P.Value, 
    adj_p_value = adj.P.Val, 
    ) %>% 
  # Write this to TSV
  readr::write_tsv(file.path(data_dir, "gene_results_GSE44971.tsv"))

