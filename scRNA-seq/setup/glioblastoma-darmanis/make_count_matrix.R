# CCDL for ALSF 2018
# C. Savonen (updated by J. Shapiro)
#
# Purpose: After salmon has been run successfully on your samples, this script
# will use tximport to quantify by transcripts and put all your samples together
# in one gene matrix tsv file.

# Options:
# "-d" - Directory of where individual samples' salmon folders are located.(default is current dir)
# "-s" - sample information file, a csv of sample data from ena (Optional)
# "-o" - Outfile for count matrix
# "-r" - Outfile for tximport RDS matrix
# "-p" - Outfile for plot
# "-m" - Percent mapped reads (reported as a decimal) cutoff for filtering
#        samples. Default is 0.5.
#
# Command line example:
#
# Rscript scripts/make_count_matrix.R \
# -d data/salmon_quants \
# -s data/sample_list.csv \
# -o data/counts.tsv \
# -r data/txi.RDS \
# -p data/mapping_plot.png \
# -m 0.5

#-------------------------- Get necessary packages-----------------------------#
# Attach needed libraries
library(optparse)


#--------------------------------Set up options--------------------------------#
option_list <- list(
    make_option(c("-d", "--dir"), type = "character", default = getwd(),
                help = "Directory where salmon quantification folders are located"),
    make_option(c("-s", "--sample-file"), type = "character", default = NULL,
                dest = "sample_file",
                help = "sample info table (csv)"),
    make_option(c("--sample-id"), type = "character", default = NULL,
                dest = "sample_id",
                help = "sample id field from sample info table to to use in counts matrix"),
    make_option(c("-t", "--tx2gene"), type = "character",
                help = "tx2gene table (tsv) for tximport"),
    make_option(c("-o", "--outfile"), type = "character",
                help = "Output file counts matrix (tsv)"),
    make_option(c("-r", "--rds"), type = "character",
                help = "Output file for tximport RDS object",
                default = NULL),
    make_option(c("-p", "--plot"), type = "character",
                help = "Output file for mapping rate plot (png)",
                default = NULL),
    make_option(c("-m", "--mapped"), type = "numeric",
                default = "0.5",
                help = "Cutoff for what percent mapped_reads samples should have.
                Any samples with less than the cutoff will be removed.")
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

#------------------------------Import Salmon reads-----------------------------#
# Get quant files
quant_files <- list.files(opt$dir, recursive = TRUE, full.names = TRUE,
                          pattern = "quant.sf")

# Get sample names from files
sample_names <- stringr::word(quant_files, -2, sep = "/")

# Get transcript IDs
transcripts <- read.table(quant_files[[1]], header = TRUE,
                          colClasses = c("character", rep("NULL", 4)))

# Read tx2gene file
tx2gene_df <- readr::read_tsv(opt$tx2gene)

# Do the thing
txi <- tximport::tximport(quant_files,
                          type = "salmon",
                          tx2gene = tx2gene_df,
                          countsFromAbundance = "no",
                          ignoreTxVersion = TRUE)

# Save to RDS if requested
if (is.character(opt$rds)){
  readr::write_rds(txi, opt$rds)
}

# Make a dataframe
counts_df <- data.frame(txi$counts, stringsAsFactors = FALSE)

# Bring the sample names with
colnames(counts_df) <- sample_names

#---------------------Salmon proportion of mapped reads------------------------#
# Get the proportion of mapped reads by reading the meta files
salmon.prop.assigned <- vapply(sample_names, function(x) {
  rjson::fromJSON(file = file.path(opt$dir, x, "aux_info",
                                   "meta_info.json"))$percent_mapped/100
  }, FUN.VALUE = 1)

# Make a histogram of this information
if (is.character(opt$plot)){
  png(opt$plot)
  hist(salmon.prop.assigned, xlab = "", main = "Proportion of Mapped Reads",
      breaks = 20)
  dev.off()
}
# Filter out samples with too low of mapped reads
counts_df <- counts_df[, which(salmon.prop.assigned > opt$mapped)]

# If a sample table is provided change sample names:
if (is.character(opt$sample_file) && is.character(opt$sample_id)){
  sample_df <- readr::read_csv(opt$sample_file) |>
    dplyr::select("run_accession", opt$sample_id) |>
    dplyr::distinct()
  # rename columns in counts_df
  counts_df <- counts_df |>
    dplyr::rename_all(function(run_id){
      tibble::tibble(run_accession = run_id) |>
        dplyr::left_join(sample_df, by = "run_accession") |>
        dplyr::pull(opt$sample_id)
      })
}

# Make a gene column so read_tsv will have the info
counts_df <- counts_df |> tibble::rownames_to_column("gene")

# Save to tsv file
readr::write_tsv(counts_df, file.path(opt$outfile))
