# C.Savonen, CCDL for ALSF
# 2018

GeneMatrixFilter <- function(dataset, min_counts = 1, num_samples = 100,
                             num_genes = 100, round_opt = FALSE) {
    # This function is filters and makes dataset into ASAP format and assumes gene info
    # is the first column
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
    #   A filtered gene expression data.frame, that is able to be submitted to ASAP
    # Get rid of decimal points
    #
    # Store the genes separately
    gene <- dataset %>% dplyr::select(dplyr::contains("gene"))

    # Get rid of gene column while we play with this
    dataset <- dataset %>% dplyr::select(-dplyr::contains("gene"))

    # round the data only if the option says to
    if (round_opt == TRUE) {
      dataset <- apply(dataset, 2, round)
    }

    # Find genes that are expressed in _% of cells
    gene.sum <- apply(dataset > min_counts, 1, sum)

    # Get indices of genes that meet criteria
    genes.keep <- which(gene.sum > num_samples)

    # Only keep those genes
    dataset <- dataset[genes.keep, ]
    gene <- gene[genes.keep, ]

    # Filter samples that express at least 100 genes
    sample.sum <- apply(dataset >= min_counts, 2, sum)

    # Keep the samples we want to keep
    dataset <- dataset[, which(sample.sum > num_genes) ]

    # Need the genes to be its own column
    dataset <- data.frame(gene, dataset)

    return(dataset)
}
