# J. Shapiro for ALSF CCDL 2021
# Clean metadata for SRP150101 bulk RNA-seq exercise notebook
#

# Input (path hardcoded): SRA run selector table txt file for SRP150101
# Output (path hardcoded): a cleaned metadata TSV

`%>%` <- dplyr::`%>%`
data_dir <- "/shared/data/training-modules/RNA-seq/data/medulloblastoma"

# files

# obtained from: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP150101&o=acc_s%3Aa
input_file <- file.path(data_dir, "SRP150101_SraRunTable.txt")
output_file <- file.path(data_dir, "SRP150101_metadata.tsv")

# read, clean, write
metadata <- readr::read_csv(input_file)
metadata <- metadata %>%
  dplyr::select(Run, Experiment, `Sample Name`, disease_state, Organism,
                Treatment, Group) %>%
  dplyr::rename(sample_name = `Sample Name`) %>%
  dplyr::mutate(
    Treatment = dplyr::case_when(
      Treatment == "Treat with digoxin" ~ "digoxin",
      Treatment == "No treatment" ~ "none",
      TRUE ~ Treatment
      ),
    # Make group into a string so and keep it that way on reimport
    Group = stringr::str_c("Group_", Group)
    ) %>%
  readr::write_tsv(output_file)
