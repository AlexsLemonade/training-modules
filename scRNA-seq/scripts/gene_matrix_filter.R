# C.Savonen, CCDL for ALSF
# 2018

gene_matrix_filter <- function(dataset, min_counts = 1, num_genes = 100,
                               num_samples = 100, round_opt = FALSE) {
    # This function filters a gene matrix that is a data.frame.
    # First it filters the samples/cells based on the number of genes a they
    # express (num_genes); then it filters genes based on how many samples
    # express it (num_samples)
    # If the matrix has a "genes" column that is labeled with "gene" in the column
    # name it will be kept in the output.
    # Args:
    #  dataset: a gene expression data.frame that is gene x samples, with the column
    #           that has gene information containing the words "gene"
    #  min_counts: the cutoff for minimum number of counts for a gene to be
    #              considered expressed in a particular sample
    #  num_samples: The minimum number of samples that need to express that gene
    #              to keep the gene in the set
    #  num_genes: The minimum number of genes a particular sample must express to
    #             be kept in the gene set
    #  round_opt: If you want the data to be rounded to no decimal points, then
    #             Put TRUE. Otherwise default is FALSE.
    # Returns:
    #   A filtered gene matrix as a data.frame based on the criteria given.
    #
    # Store the genes separately
    gene <- dataset %>% dplyr::select(dplyr::contains("gene"))

    # Get rid of gene column while we calculate things
    dataset <- dataset %>% dplyr::select(-dplyr::contains("gene"))

    # round the data only if the option says to
    if (round_opt == TRUE) {
      dataset <- apply(dataset, 2, round)
    }

    # Calculate how many genes a sample expresses based on min_counts
    sample.sum <- apply(dataset > min_counts, 2, sum)

    # Keep the samples we want to keep
    dataset <- dataset[, which(sample.sum > num_genes)]

    # Calculate how many samples express a particular gene
    gene.sum <- apply(dataset > min_counts, 1, sum)

    # Get indices of genes that meet criteria
    genes.keep <- which(gene.sum > num_samples)

    # Only keep those genes
    dataset <- dataset[genes.keep, ]
    gene <- gene[genes.keep, ]

    # Reattach the gene names column
    dataset <- data.frame(gene, dataset)

    return(dataset)
}
