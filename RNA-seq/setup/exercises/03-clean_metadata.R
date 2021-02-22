# J. Taroni for ALSF CCDL 2019
# Clean metadata for SRP049821 bulk RNA-seq exercise notebook
#
# Command line usage: Rscript 03-clean_metadata.R
# Input (path hardcoded): SRA run selector table txt file for SRP049821
# Output (path hardcoded): a cleaned metadata TSV

`%>%` <- dplyr::`%>%`
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# files
sample_info_dir <- file.path(root_dir, "bulk-rnaseq", "data", "sample_metadata")
output_dir <- file.path(sample_info_dir, "leukemia")
# obtained from: https://www.ncbi.nlm.nih.gov/Traces/study/?WebEnv=NCID_1_19093743_130.14.22.33_5555_1559150230_374684744_0MetA0_S_HStore&query_key=4
input_file <- file.path(sample_info_dir, "SRP049821_SraRunTable.txt")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(output_dir, "SRP049821_metadata.tsv")

# read, clean, write
readr::read_tsv(input_file) %>%
  dplyr::select(Run, Experiment, SRA_Sample, BioSample, cell_sorting, Organism,
                cell_type, genotype_variation, strain) %>%
  dplyr::mutate(cell_sorting = gsub("-", "neg",
                                    gsub("\\+", "pos", cell_sorting))) %>%
  readr::write_tsv(output_file)
