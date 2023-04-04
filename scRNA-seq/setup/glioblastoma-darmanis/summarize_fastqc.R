# CCDL ALSF 2018
# C. Savonen
#
# Purpose: Get a summary file of individual sequence quality reports
# from the FASTQC files.
#
# Options:
# "-d" - directory/path of where the fastqc reports have been placed.
# "-t" - path for the output summary table
# "-f" - path for the filtered output summary (only passing runs)
#
# Example bash usage:
#
# Rscript scripts/2-get_fastqc_reports.R \
# -d data/fastqc_reports \
# -t fastqc_report.csv
# -f fastqc_filtered.csv
#
library(optparse)
library(fastqcr)

# Get options using optparse
option_list <- list(
  make_option(opt_str = c("-d", "--dir"), type = "character",
              help = "Directory containing the fastqc reports"),
  make_option(opt_str = c("-t", "--table"), type = "character",
              help = "Report table path (csv)"),
  make_option(opt_str = c("-f", "--filtered"), type = "character",
              help = "Report table for passing runs only (csv)")
  )

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))


# Aggregate the reports
qc <- fastqcr::qc_aggregate(qc.dir = opt$dir)

# Write full report table to a csv file
readr::write_csv(qc_stats(qc), path = opt$table)

# Filter out samples that have failed the quality tests
qc_filtered <- data.frame(qc) |>
  dplyr::select(sample, module, status) |>
  dplyr::filter(status %in% "PASS") |>
  dplyr::arrange(sample)

# Write filtered results to a csv file
readr::write_csv(qc_filtered, path = opt$filtered)
