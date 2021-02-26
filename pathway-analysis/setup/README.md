## Pathway analysis

This document provides background for preparing data for use in the guided notebooks of the pathway analysis module.

### Dependence on Other Modules

⚠️   **The setup of this module relies on data or results prepared as part of other modules.**

For the GSEA notebook, we use preranked values (log2 fold change) from the neuroblastoma cell line differential gene expression analysis covered in `RNA-seq/05-nb_cell_line_DESeq2`.
`05-nb_cell_line_DESeq2` relies on the tximeta output produced through the bulk RNA-seq guided exercise (`RNA-seq/04-nb_cell_line_tximeta.md`).

For the GSVA notebook, we use a subset of the OpenPBTA stranded dataset (medulloblastoma samples) that we download and process in `machine-learning/setup`. It also depends on a step, performed in `machine-learning/03-openpbta_PLIER.Rmd`, where duplicate gene symbols are collapsed to the row that has the highest mean value.

### Data Preparation

`01-prepare_NB_cell_line.Rmd` applies `DESeq2::lfcShrink()` to the neuroblastoma cell line results and saves the results table to `pathway-analysis/results/gene-metrics/nb_cell_line_mycn_amplified_v_nonamplified.tsv`.

`02-prepare_openpbta_MB_data.Rmd` takes the transformed and filtered OpenPBTA data (processed in `machine-learning/setup`) from `machine-learning/data/open-pbta/processed` and subsets the RNA-seq data and metadata files to medulloblastoma samples only.
