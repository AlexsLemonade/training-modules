# J. Shapiro for ALSF CCDL 2021
# Clean metadata for SRP150101 bulk RNA-seq exercise notebook
#

# Input (path hardcoded): SRA run selector table txt file for SRP150101
# Output (path hardcoded): a cleaned metadata TSV

`%>%` <- dplyr::`%>%`
data_dir <- "/shared/data/training-modules/RNA-seq/data/medulloblastoma"

# files

# obtained from: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP150101&o=acc_s%3Aa
sra_input_file <- file.path(data_dir, "SRP150101_SraRunTable.txt")
# obtained from refine.bio on 2021-03-06
refinebio_input_file <- file.path(data_dir, "metadata_SRP150101.tsv")
output_file <- file.path(data_dir, "SRP150101_metadata.tsv")

# two metadata data frames enter, one metadata data frame leaves
sra_metadata <- readr::read_csv(sra_input_file)
refinebio_metadata <- readr::read_tsv(refinebio_input_file)

# join the two data frames of metadata
metadata <- sra_metadata %>%
  # we'll get most of our fields from the SRA file
  dplyr::select(Run, Experiment, `Sample Name`, disease_state, Organism,
                Treatment, Group) %>%
  # the refine.bio *title* is what tells us about technical replicates
  dplyr::inner_join(dplyr::select(refinebio_metadata,
                                  refinebio_accession_code,
                                  refinebio_title),
                    by = c("Run" = "refinebio_accession_code"))

# clean!
metadata <- metadata %>%
  dplyr::rename(sample_name = `Sample Name`) %>%
  dplyr::mutate(
    # everything but the last 5 characters is the mouse ID
    mouse_id = stringr::str_sub(metadata$refinebio_title, 1, -6),
    # simplify mouse ids
    Mouse = forcats::fct_anon(as.factor(mouse_id), "M"),
    # the last four characters of the title tells us about the lane
    replicate_id = stringr::str_sub(metadata$refinebio_title, -4),
    Treatment = dplyr::case_when(
      Treatment == "Treat with digoxin" ~ "digoxin",
      Treatment == "No treatment" ~ "none",
      TRUE ~ Treatment
      ),
    # Make group into a string so and keep it that way on reimport
    Group = stringr::str_c("Group_", Group)
    ) %>%
  # we no longer need the title
  dplyr::select(-refinebio_title) %>%
  readr::write_tsv(output_file)
