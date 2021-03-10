# Catch dependencies that renv might miss

# Recommended by ComplexHeatmap
library(magick)

# Recommended by tximport
library(fishpond)

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

# SRAdb
library(SRAdb)
