# Cooking show magic

To allow participants who might get bogged down to quickly recover and continue with later notebook sections or exercises, we like to keep a set of output files on hand that they can quickly access as needed. 
This document is for tracking where these files are created and stored.

In general, the files for this module will be stored in `/shared/data/training-data/<date>/scRNA-seq-advanced/`, and should be created before each training using the latest versions of the referenced notebooks. 
The sub-path for each file should be the same as that for `training-modules/scRNA-seq-advanced/`

Some of these files will also be copied to S3 by the `sync-s3.sh` script for use by the `exercise-notebook-answers` repository for testing.

Files are listed below by the notebook that produces them:

- 01-read_filter_normalize_scRNA.Rmd
  - `data/glioblastoma/normalized/glioblastoma_normalized_sce.rds`
- 02-celltype_assignment.Rmd 
  - `analysis/PBMC-TotalSeqB/PBMC_TotalSeqB_cellinfo.tsv`
- 03-dataset_integration.Rmd 
  - `data/rms/integrated/rms_integrated_subset.rds`
- 04-differential_expression.Rmd
  - `analysis/rms/deseq/rms_myoblast_deseq_results.tsv`
  - `data/rms/integrated/rms_subset_sce.rds`
  
