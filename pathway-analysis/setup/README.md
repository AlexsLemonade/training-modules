## Pathway analysis

This document provides background for preparing data for use in the guided notebooks of the pathway analysis module.

### Dependence on Other Modules

⚠️   **The setup of this module relies on data or results prepared as part of other modules.**

For the GSVA notebook, we use a subset of the OpenPBTA stranded dataset (medulloblastoma samples) that we download and process in `machine-learning/setup`. It also depends on a step, performed in `machine-learning/03-openpbta_PLIER.Rmd`, where duplicate gene symbols are collapsed to the row that has the highest mean value.

### Data Preparation

`01-leukemia_DGE.Rmd` includes two comparisons that generate differential gene expression (DGE) results used in the over-representation analysis notebook.
It uses the same dataset as the `RNA-seq/exercise_02-bulk_rnaseq.Rmd` notebook.

`02-prepare_openpbta_MB_data.Rmd` takes the transformed and filtered OpenPBTA data (processed in `machine-learning/setup`) from `machine-learning/data/open-pbta/processed` and subsets the RNA-seq data and metadata files to medulloblastoma samples only.
