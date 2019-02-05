# C. Savonen
# CCDL for ALSF 2019

# Purpose: Function to read in data from alevin output 
# Usage: to import into R environment using "source"

# Requirements: need "readr" to be installed.

# This function is based on the tutorial from the COMBINE lab
ReadAlevin <- function(base_path = NULL){
  # Function: to read in alevin output data into R
  # Args: base_path : the relative path to the directory where the alevin output
  #                   is saved.
  # Output: a dataframe of the alevin counts data where genes' data are by rows 
  #         and cells' data are by columns. 
  if (! dir.exists(base_path )){
    stop("Directory provided does not exist")
  }
  # Obtain paths to data
  barcode_loc <- file.path(base_path, "alevin", "quants_mat_rows.txt")
  gene_loc <- file.path(base_path, "alevin", "quants_mat_cols.txt")
  data_loc <- file.path(base_path, "alevin", "quants_mat.csv")
  
  if (!file.exists(barcode_loc)) {
    stop("Barcode file missing")
  }
  if (!file.exists(gene_loc)) {
    stop("Gene name file missing")
  }
  if (!file.exists(data_loc)) {
    stop("Expression matrix file missing")
  }
  # Read in the data from Alevin output
  expression_matrix <- readr::read_csv(data_loc, col_names = FALSE, 
                                       progress = FALSE)
  
  # Transpose the matrix so it is gene x cell
  expression_matrix <- t(expression_matrix[,1:ncol(expression_matrix)-1])
  
  # Apply the colnames and rownames to the dataset
  colnames(expression_matrix) <- readLines(barcode_loc)
  rownames(expression_matrix) <- readLines(gene_loc)
  
  # Make NA values into 0's
  expression_matrix[is.na(expression_matrix)] <- 0
  return(expression_matrix)
}
