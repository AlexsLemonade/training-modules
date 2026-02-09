# Catch dependencies that renv might miss

# Don't know why markdown is missed, but it is
library(markdown)

# Recommended by ComplexHeatmap
library(magick)

# Recommended by tximport
library(eds)

# Required for DESeq2::lfcShrink()
library(apeglm)
library(ashr)

# Required for scater
library(uwot)
library(Rtsne)

# required for enrichplot::upsetplot
library(ggupset)

# Organisms we want to support but don't explicitly work with during instruction
# or exercises
library(org.Dr.eg.db)
library(org.Cf.eg.db)

# Required for vsn plots in RNA-seq exercise
library(hexbin)

# AlevinQC was not getting noticed for some reason
library(alevinQC)

# Doublet Detection for exercise notebook
library(scDblFinder)

# Loom file format functions for Single Cell data
library(LoomExperiment)

# Needed for SingleR to run de with wilcox, for cell type exercises
library(scrapper)
