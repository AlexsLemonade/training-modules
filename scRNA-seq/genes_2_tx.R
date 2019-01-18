# C. Savonen 
# CCDL for ALSF 
# 2019

# Purpose: Make transcript gene key list 

# Get ensembl genes
genes <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, "ENSEMBL")

# Retrieve transcripts for those
gene_2_tx <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, 
                                  keys = genes, 
                                  column = "ENSEMBLTRANS", 
                                  keytype = "ENSEMBL", 
                                  multiVals="CharacterList")

# Make into dataframe:
gene_2_tx <- reshape2::melt(gene_2_tx@listData)

# Get rid of transcripts without genes
gene_2_tx <- gene_2_tx[!is.na(gene_2_tx$value), ]

# Write to tsv file
readr::write_tsv(gene_2_tx, "genes_2_tx.tsv", col_names = FALSE)

# Read in ensembl transcriptome
genome <- readLines("data/ref_files/Homo_sapiens.GRCh38.cdna.all.fa")
genome <- genome[grep(">", genome)]

# Extract transcripts
transcripts <- stringr::word(genome, sep = " ", 1)
transcripts <- gsub(">", "", transcripts)

# Extract genes
genes <- stringr::word(genome, sep = "gene:", 2)
genes <- stringr::word(genes, sep = " ", 1)
# genes <- gsub("\\.[0-9]*$", "", genes)

readr::write_tsv(data.frame(transcripts, genes), "genes_2_tx.tsv",
                 col_names = FALSE)
