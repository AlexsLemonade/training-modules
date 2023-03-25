# J. Taroni for ALSF CCDL 2019, modified by J Shapiro 2021
# Clean metadata for SRP049821 bulk RNA-seq exercise notebook
#
# Command line usage: Rscript 01-clean_metadata.R
# Input (path hardcoded): SRA run selector table txt file for SRP049821
# Output (path hardcoded): a cleaned metadata TSV

data_dir <- "/shared/data/training-modules/RNA-seq/data/leukemia"

# files

# obtained from: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP049821&o=acc_s%3Aa
input_file <- file.path(data_dir, "SRP049821_SraRunTable.txt")
output_file <- file.path(data_dir, "SRP049821_metadata.tsv")

# read, clean, write
metadata <- readr::read_csv(input_file)
metadata <- metadata |>
  dplyr::select(Run, Experiment, `Sample Name`, BioSample, cell_sorting, Organism,
                Cell_type, `genotype/variation`, Strain) |>
  dplyr::rename(sample_name = `Sample Name`,
                genotype_variation = `genotype/variation`) |>
  dplyr::mutate(cell_sorting = gsub("-", "neg",
                                    gsub("\\+", "pos", cell_sorting))) |>
  readr::write_tsv(output_file)
