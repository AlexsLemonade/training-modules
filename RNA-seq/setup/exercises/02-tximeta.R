# J. Shapiro for ALSF CCDL 2021 (modified from J. Taroni)
# Process SRP049821 RNA-seq data with tximeta for bulk RNA-seq exercise
# notebook
#
# Command line usage: Rscript 02-tximeta.R
# Input (paths hardcoded): quant files from SRP049821
# Output (path hardcoded): output of tximeta (gene summarized)

data_dir <- "/shared/data/training-modules/RNA-seq/data/leukemia"

meta_file <- file.path(data_dir, "SRP049821_metadata.tsv")

# create a directory to hold the tximport output if it does not yet exist
output_directory <- file.path(data_dir, "tximport")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

# salmon output -- need sample names to easier downstream use
quant_files <- list.files(file.path(data_dir, "salmon_quant"),
                          pattern = "quant.sf",
                          recursive = TRUE, full.names = TRUE)
sample_names <- stringr::word(quant_files, -2, sep = "/")
names(quant_files) <- sample_names

coldata <- data.frame(files = quant_files, 
                      names = sample_names)

## read in metadata 
metadata <- readr::read_tsv(meta_file)

coldata <- dplyr::inner_join(coldata, metadata, 
                             by = c("names" = "Run"))

# tximport + save to file
txi <- tximeta::tximeta(coldata)
gene_txi <- tximeta::summarizeToGene(txi_data) 
saveRDS(gene_txi, file = file.path(output_directory, "leukemia_stem_cell_txi.RDS"))
