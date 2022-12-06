This directory holds the workflow for pre-processing integration SCE files.

For local testing, you may want to populate the `raw/` subdirectory with the following files downloaded from `ScPCA`:

- `SCPCL000478_filtered.rds`
- `SCPCL000479_filtered.rds`
- `SCPCL000480_filtered.rds`
- `SCPCL000481_filtered.rds`

In addition, the `raw/` directory should contain a file `celltypes.tsv` which was compiled based on cell types annotated in [Patel et al. (2022)](https://doi.org/10.1016/j.devcel.2022.04.003).
This file has three columns representing `ScPCA` libraries, cell barcodes, and annotated cell types.
