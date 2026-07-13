#!/usr/bin/env RScript

# Script to format annotations file for Janesick Xenium data

library(optparse)

option_list <- list(
  make_option("--input_file",
              type = "character",
              help = "Path to the 10x cell type annotations .xlsx"),
  make_option("--output_file",
              type = "character",
              help = "Path to write the formatted annotations .tsv"),
  make_option("--sheet",
              type = "integer",
              default = 4, # use annotations transferred from Flex data for rep1
              help = "Worksheet to read [default %default = rep1]")
)
opts <- parse_args(OptionParser(option_list = option_list))

# read in the appropriate sheet and rename columns for easier joining
# also we should name the cell type annotation column appropriately
cell_types_df <- openxlsx::read.xlsx(opts$input_file, sheet = opts$sheet) |>
  dplyr::rename(cell_id = Barcode,
                cell_type_annotation = Cluster)

readr::write_tsv(cell_types_df, opts$output_file)
