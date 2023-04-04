#!/usr/bin/env Rscript

# Creates and prints a table of metadata from GEO for a given acession

library(optparse)

option_list <- list( 
  make_option(c("-g", "--geo"), default = NULL,
              help = "GEO ID of the RNA-seq data you would like to download
              eg. GSE84465"),
  make_option(c("-o", "--outfile"), default = "metadata.tsv",
              help = "Location for table, which will have tsv format")
)
opt <- parse_args(OptionParser(option_list = option_list))


# Get metadata from GEO
geo_meta <- GEOquery::getGEO(opt$geo)
geo_df <- data.frame(geo_meta[[1]]@phenoData@data)

readr::write_tsv(geo_df, opt$outfile)
