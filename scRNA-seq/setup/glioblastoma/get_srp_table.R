#!/usr/bin/env Rscript

# Creates and prints a table of all samples for the given SRP ID and
# their RNA download links.

library(optparse)
library(magrittr)

option_list <- list( 
  make_option(c("-a", "--srp"), default = NULL,
              help = "SRP ID of the RNA-seq data you would like to download
              eg. SRP079058"),
  make_option(c("-o", "--outfile"), default = "sample_list.csv",
              help = "Location for table, which will have csv format"),
  make_option(c("-p", "--paired"), action ="store_true", 
              help = "Check that sample files are paired. 
              If there are not exactly 2 files, the sample will be dropped")
)
opts <- parse_args(OptionParser(option_list = option_list))

# gets a table of ENA ftp links from an SRA Project accession id
# returns a data frame with columns: run_accession, fastq_ftp, fastq_md5, & fastq_bytes
get_sample_df <- function(SRA_project){
  ena_url <- stringr::str_interp("http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=${SRA_project}&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes")
  ena_table <- readr::read_tsv(ena_url) %>%
    tidyr::separate_rows(fastq_ftp, fastq_md5, fastq_bytes, sep = ";")
}

sample_df <- get_sample_df(opts$srp)
if (opts$paired){
  sample_df <- sample_df %>%
    dplyr::group_by(run_accession) %>%
    dplyr::mutate(n = n()) %>%
    dplyr::filter(n == 2) %>%
    dplyr::select(-n) %>%
    dplyr::ungroup()
}

readr::write_csv(sample_df, opts$outfile)