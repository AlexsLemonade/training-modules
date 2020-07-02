# CCDL for ALSF 2018
# C. Savonen (updated by J. Shapiro)
#
# Purpose: After salmon has been run successfully on your samples, this script
# will use tximport to quantify by transcripts and put all your samples together
# in one gene matrix tsv file.

# Options:
# "-d" - Directory of where individual samples' salmon folders are located.(Optional)
# "-o" - Directory of where the output gene matrix RDS file should go.(Optional)
# "-m" - Percent mapped reads (reported as a decimal) cutoff for filtering
#        samples. Default is 0.5.
#
# Command line example:
#
# Rscript scripts/make_count_matrix.R \
# -d data/salmon_quants \
# -o data/counts.tsv \
# -m 0.5 \

#-------------------------- Get necessary packages-----------------------------#
# Attach needed libraries
library(optparse)

# Magrittr pipe
`%>%` <- dplyr::`%>%`

#--------------------------------Set up options--------------------------------#
option_list <- list(
    make_option(c("-d", "--dir"), type = "character", default = getwd(),
                help = "Directory where salmon quantification folders are located"),
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
                Any samples with less than the cutoff will be removed."),
    make_option(c("-l", "--label"), type = "character",
                default = "", 
                help = "Optional label for output files")
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Add an underscore if label is specified
if (!is.null(opt$label)){
  opt$label <- paste0(opt$label, "_")
}

#------------------------------Import Salmon reads-----------------------------#
# Get quant files
quant_files <- list.files(opt$dir, recursive = TRUE, full.names = TRUE,
                          pattern = "quant.sf")

# Get sample names
sample_names <- stringr::word(quant_files, -2, sep = "/")

# Get transcript IDs
transcripts <- read.table(quant_files[[1]], header = TRUE,
                          colClasses = c("character", rep("NULL", 4)))

# Get rid of transcript version numbers because mapIDs will not recognize them
transcripts.wout.ver <- gsub("\\.[0-9]*$", "", transcripts$Name)

# Read tx2gene file
tx2gene_df <- readr::read_tsv(opt$tx2gene)

# Do the thing
txi <- tximport::tximport(quant_files, 
                          type = "salmon", 
                          tx2gene = tx2gene_df,
                          countsFromAbundance = "no",
                          ignoreTxVersion = TRUE)

# Save to RDS file temporarily
if (is.character(opt$rds)){
  readr::write_rds(txi, opt$rds)
}

# counts_df <- readRDS("tximport_obj.RDS")

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

# Make a gene column so read_tsv will have the info
counts_df <- counts_df %>% tibble::rownames_to_column("gene")

# Save to tsv file
readr::write_tsv(counts_df, file.path(opt$outfile))
