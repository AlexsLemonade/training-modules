
plot_AUC <- function(cell_rankings,
                     gene_set_collection,
                     gene_set_name,
                     barcode,
                     top_threshold = 0.05,
                     seed = 2024) {

  # Adapted from https://github.com/aertslab/AUCell/blob/91753b327a39dc7a4bbed46408ec2271c485f2f0/vignettes/AUCell.Rmd#L295-L316

  set.seed(seed)

  # Pull out the gene set and identify where in the cell ranks those genes
  # lie
  gene_set <- gene_set_collection[[gene_set_name]]
  gene_set_ranks <- cell_rankings[geneIds(gene_set), ]

  # Index of the cell with the barcode
  barcode_index <- which(colnames(gene_set_ranks) == barcode)

  # The max rank for computing the AUC will be determined by the threshold
  aucMaxRank <- nrow(cell_rankings) * top_threshold
  sorted_ranks <- sort(getRanking(gene_set_ranks[, barcode_index]))
  aucCurve <- cbind(c(0, sorted_ranks, nrow(cell_rankings)), c(0:length(sorted_ranks), length(sorted_ranks)))
  plot(aucCurve,
       type="s", col="darkblue", lwd=1,
       xlab="Gene rank", ylab="# genes in the gene set",
       xlim=c(0, aucMaxRank*2), ylim=c(0, nGenes(gene_set)),
       main="Recovery curve",
       sub=paste("Cell:", colnames(gene_set_ranks)[barcode_index]))
  aucShade <- aucCurve[which(aucCurve[,1] < aucMaxRank),]
  aucShade <- rbind(aucShade, c(aucMaxRank, nrow(aucShade)))
  aucShade[,1] <-  aucShade[,1]-1
  aucShade <- rbind(aucShade, c(max(aucShade),0))
  polygon(aucShade, col="#0066aa40", border=FALSE)

  abline(v=aucMaxRank, lty=2)
  text(aucMaxRank-50, 5, "AUC")

}
