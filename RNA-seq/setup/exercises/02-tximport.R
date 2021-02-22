# J. Taroni for ALSF CCDL 2019
# Process SRP049821 RNA-seq data with tximport for bulk RNA-seq exercise
# notebook
#
# Command line usage: Rscript 02-tximport.R
# Input (paths hardcoded): quant files from SRP049821, tx2gene mapping
# Output (path hardcoded): output of tximport

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# create a directory to hold the tximport output if it does not yet exist
output_directory <- file.path(root_dir, "bulk-rnaseq", "data", "tximport",
                              "leukemia")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

# salmon output -- need sample names to easier downstream use
quant_files <- list.files(file.path(root_dir, "bulk-rnaseq", "data", "quants"),
                          pattern = "quant.sf",
                          recursive = TRUE, full.names = TRUE)
sample_names <- stringr::word(quant_files, 3, sep = "/")
names(quant_files) <- sample_names

# transcript to gene mapping
# from https://github.com/AlexsLemonade/training-txome-prep/blob/9b1c326f4a0431d080e1ad19c0718b038ba32a69/tx2gene/Mus_musculus.GRCm38.95_tx2gene.tsv
tx2gene_file <- file.path("..", "index", "Mus_musculus.GRCm38.95_tx2gene.tsv")
tx2gene_df <- readr::read_tsv(tx2gene_file)

# tximport + save to file
txi <- tximport::tximport(quant_files, type = "salmon", tx2gene = tx2gene_df,
                          countsFromAbundance = "no", ignoreTxVersion = TRUE)
saveRDS(txi, file = file.path(output_directory, "leukemia_stem_cell_txi.RDS"))
