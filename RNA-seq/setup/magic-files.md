# Cooking show magic

To allow participants who might get bogged down to quickly recover and continue with later notebook sections or exercises, we like to keep a set of output files on hand that they can quickly access as needed.
This document is for tracking where these files are created and stored.

In general, the files for this module will be stored in `/shared/data/training-data/<date>/RNA-seq/`, and should be created before each training using the latest versions of the referenced notebooks.
The sub-path for each file should be the same as that for `training-modules/RNA-seq`.

Files are listed below by the notebook that produces them:

- `02-gastric_cancer_tximeta.Rmd`
  - `data/gastric-cancer/txi/gastric-cancer_tximeta.rds`
- `setup/ref_notebooks/04-nb_cell_line_tximeta.Rmd` (participants create their own notebook for this step)
  - `data/NB-cell/txi/NB-cell_tximeta.rds`
- `05-nb_cell_line_DESeq2.Rmd`
  - `results/NB-cell/NB-cell_DESeq_amplified_v_nonamplified.rds`
  - `results/NB-cell/NB-cell_DESeq_amplified_v_nonamplified_results.tsv`



Note that the output files from `scripts/run_SRR585570.sh` script are already present at:

- `/shared/data/training-modules/RNA-seq/data/salmon-quant/SRR585570/`
- `/shared/data/training-modules/RNA-seq/QC/gastric-cancer/fastp/SRR585570/`
- `/shared/data/training-modules/RNA-seq/QC/gastric-cancer/fastqc/SRR585570/`

so do not need to be added to the `training-data` directory, could be needed during training.
