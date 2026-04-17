# Set up gene level statistics table for intro to R module

# C. Savonen for ALSF - CCDL

# 2020

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
                                     "GSE44971.tsv")) |>
  tibble::column_to_rownames("Gene")

# Read in metadata
metadata <- readr::read_tsv(file.path(data_dir,
                                      "metadata_GSE44971.tsv")) |>
  dplyr::select(sample_id = refinebio_accession_code,
                -which(apply(is.na(.), 2, all)),
                sex = refinebio_sex,
                tissue,
                brain_location = "brain location",
                title) |>
  dplyr::mutate(
    # Clean up character names
    tissue = gsub(" ", "_", tolower(tissue))
    )

# Write to a TSV file
readr::write_tsv(metadata,
                 file.path(data_dir, "cleaned_metadata_GSE44971.tsv"))

## Put samples' data and metadata in the same order
# Make the data in the order of the metadata
gene_df <- gene_df |>
  dplyr::select(metadata$sample_id)

############# Set up combined sex and tissue variable for model
# This set up is based on the sum to zero parametrization example from
# 9.5.4 Classic Interaction Models in this manual:
# https://bioconductor.org/packages/3.1/bioc/vignettes/limma/inst/doc/usersguide.pdf

# We will test for sex and tissue variables.
# We need these to be factor variables
sex <- factor(metadata$sex)
tissue <- factor(metadata$tissue)

# Here we are declaring our interaction model, this will give us a matrix
# of zeroes and ones that corresponds to which samples belong to which group
design <- model.matrix(~sex*tissue)

# Neaten up column names
colnames(design) <- c("intercept",
                      "male_female",
                      "astrocytoma_normal",
                      "interaction")

# Apply linear model to data
fit <- limma::lmFit(gene_df, design)

# Apply empirical Bayes to smooth standard errors
fit2 <- limma::eBayes(fit)

########## Extract a results table for each contrast we care about
# limma::topTable function will apply multiple testing correction and obtain
# summary statistics on each

# This first table extracts the results for male_female contrast
male_female <- limma::topTable(fit2,
                               coef = "male_female", # This is what we named this contrast
                                                     # in the test
                               number = Inf, # We can tell limma how many results we want back.
                                             # we want them all, so we said `Inf``
                               sort = "none" # We don't want the data sorted, but default is to sort
                               ) |>
  tibble::rownames_to_column("ensembl_id") # We want the ensembl_id to be its own column
                                           # for each of these tables, its easier for storage

# Extract test results for astrocytoma tumor and normal
astrocytoma_normal <- limma::topTable(fit2,
                                      coef = "astrocytoma_normal",
                                      number = Inf,
                                      sort = "none") |>
  tibble::rownames_to_column("ensembl_id")

# Extract test results for the interaction of sex and tissue
interaction <- limma::topTable(fit2,
                               coef = "interaction",
                               number = Inf,
                               sort = "none") |>
  tibble::rownames_to_column("ensembl_id")

# Bind together the 3 results tables for each contrast into one big table
# dplyr::bind_rows will match by column names and bind the rows.
stats_df <- dplyr::bind_rows(
  "male_female" = male_female,
  "astrocytoma_normal" = astrocytoma_normal,
  "interaction" = interaction,
  .id = "contrast" # This argument will create a column that labels for us which results
                   # data.frame it is originally from
)

# Neaten up this data.frame in preparation for saving to a TSV file
stats_df <- stats_df |>
  # Map ensembl IDs to their associated gene symbols by default we are only
  # taking the first mapped value. Can change this with the multiVals argument.
  dplyr::mutate("gene_symbol" = AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = ensembl_id,
    column = "SYMBOL",
    keytype = "ENSEMBL")) |>
  # Clean up column names and reorder
  dplyr::select(
    ensembl_id,
    gene_symbol,
    contrast,
    log_fold_change = logFC,
    avg_expression = AveExpr, # We want our column names to be consistent format
    t_statistic = t, # There is a function called `t` so for disambiguation purposes, we will name this t_statistic
    p_value = P.Value,
    adj_p_value = adj.P.Val,
    ) |>
  # Add filter
  dplyr::filter(avg_expression > 5) |>
  # Write this to TSV
  readr::write_tsv(file.path(data_dir, "gene_results_GSE44971.tsv"))
