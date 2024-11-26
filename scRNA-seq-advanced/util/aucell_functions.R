#' Plot an AUCell recovery code, coloring in the area under the curve (AUC)
#'
#' Adapted from https://github.com/aertslab/AUCell/blob/91753b327a39dc7a4bbed46408ec2271c485f2f0/vignettes/AUCell.Rmd#L295-L316
#'
#' @param cell_rankings Cell rankings; the output of AUCell::AUCell_buildRankings()
#' @param gene_set_collection A GSEABase GeneSetCollection object
#' @param gene_set_name The name of the gene set from the GeneSetCollection
#'   that you would like to plot a recovery curve for
#' @param barcode The cell barcode that you would like to plot a recovery
#'   curve for
#' @param auc_max_rank The maximum gene rank you would like to use for the
#'   area under the curve
#' @param seed Seed passed to set.seed(); 2024 by default
#'
#' @return Outputs a recovery curve plot with the
plot_recovery_curve <- function(cell_rankings,
                                gene_set_collection,
                                gene_set_name,
                                barcode,
                                auc_max_rank,
                                seed = 2024) {

  set.seed(seed)

  # Pull out the gene set and identify where in the cell ranks those genes
  # lie
  gene_set <- gene_set_collection[[gene_set_name]]
  gene_set_ranks <- cell_rankings[geneIds(gene_set), ]

  # Index of the cell with the barcode
  barcode_index <- which(colnames(gene_set_ranks) == barcode)

  sorted_ranks <- sort(getRanking(gene_set_ranks[, barcode_index]))
  aucCurve <- cbind(c(0, sorted_ranks, nrow(cell_rankings)), c(0:length(sorted_ranks), length(sorted_ranks)))
  plot(aucCurve,
       type="s", col="darkblue", lwd=1,
       xlab="Gene rank", ylab="# genes in the gene set",
       xlim=c(0, auc_max_rank*3), ylim=c(0, nGenes(gene_set)),
       main="Recovery curve",
       sub=paste("Cell:", colnames(gene_set_ranks)[barcode_index]))
  aucShade <- aucCurve[which(aucCurve[,1] < auc_max_rank),]
  aucShade <- rbind(aucShade, c(auc_max_rank, nrow(aucShade)))
  aucShade[,1] <-  aucShade[,1]-1
  aucShade <- rbind(aucShade, c(max(aucShade),0))
  polygon(aucShade, col="#0066aa40", border=FALSE)

  abline(v=auc_max_rank, lty=2)
  text(auc_max_rank-50, 5, "AUC")

}
