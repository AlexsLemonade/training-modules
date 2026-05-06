# Spatial workshop training data setup

This document describes how the training data is prepared for the spatial workshop on the RStudio Server.

## File locations

On the RStudio server the main location for the files needed for this training module will be `/shared/data/training-modules/spatial/data`.
The files are then organized by dataset.

After setup, the symlinks should be established from `/shared/data/training-modules/spatial/data` to `training-modules/spatial/data` as appropriate.
These links will usually be created using the script at `training-modules/scripts/link-data.sh`, so directories and files should be added to the relevant `link_locs` array in that script.
If no data will be written to a data directory, linking to that directory will be sufficient, but in some cases links to individual files may be required.


## Data Sets

### Mitochondrial gene lists

We use use a table of mitochondrial genes previously created in the `../../scRNA-seq-advanced/setup/mito_gene_lists.Rmd` notebook for this workshop.
This table was copied from `/shared/data/training-modules/scRNA-seq-advanced/data/reference` to `/shared/data/training-modules/spatial/data/reference`.

### Wilms tumor

We use the `SCPCL000429_spatial` dataset from [`SCPCP000006`](https://scpca.alexslemonade.org/projects/SCPCP000006).
To obtain these files, use the script `wilms-tumor/download-wilms-tumor.R`.
This script requires the [`ScPCAr` package](https://alexslemonade.github.io/ScPCAr/), which can be installed with `remotes::install_github("AlexsLemonade/ScPCAr")`.

This will create a directory `../data/wilms-tumor/SCPCS00190/` with the following files:
```console
SCPCS000190
└── spaceranger
    ├── raw_feature_bc_matrix
    │   ├── barcodes.tsv.gz
    │   ├── features.tsv.gz
    │   └── matrix.mtx.gz
    ├── SCPCL000429_spaceranger-summary.html
    └── spatial
        ├── aligned_fiducials.jpg
        ├── detected_tissue_image.jpg
        ├── scalefactors_json.json
        ├── tissue_hires_image.png
        ├── tissue_lowres_image.png
        └── tissue_positions_list.csv
```

The `SCPCS000190` directory was then copied to `/shared/data/training-modules/spatial/data/wilms-tumor/`.


### Ovarian carcinoma

This data comes from this 10x Genomics dataset: <https://www.10xgenomics.com/datasets/human-ovarian-cancer-11-mm-capture-area-ffpe-2-standard>.
This data, as well as its corresponding Visium probe set, can be downloaded with `ovarian-carcinoma/download-ovarian.R`.
This will create a directory `../data/ovarian-carcinoma/spaceranger/` with the following files:

```console
spaceranger
├── filtered_feature_bc_matrix
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── spatial
│   ├── aligned_fiducials.jpg
│   ├── aligned_tissue_image.jpg
│   ├── cytassist_image.tiff
│   ├── detected_tissue_image.jpg
│   ├── scalefactors_json.json
│   ├── spatial_enrichment.csv
│   ├── tissue_hires_image.png
│   ├── tissue_lowres_image.png
│   └── tissue_positions.csv
└── Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv
```

The `ovarian-carcinoma` directory was then copied to `/shared/data/training-modules/spatial/data/`.

### Osteosarcoma

This data comes from this `GEO` record: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8478586>, which is associated with [Reinecke et al. (2025)](https://aacrjournals.org/clincancerres/article/31/2/414/751106/Aberrant-Activation-of-Wound-Healing-Programs).
This data needs to be filtered and normalized for input to workshop notebooks.

In addition, the [OsteoCAR mouse metastasis reference](https://figshare.com/articles/dataset/OsteoCAR_A_multi-species_single-cell_atlas_of_primary_and_metastatic_osteosarcoma/31029559) also needs to be obtained and prepared for the workshop as the osteosarcoma deconvolution reference.

This code requires the `qs` (not `qs2`) package to read in the raw osteo reference.
As of January 2026, this package is not currently from CRAN (but it may be back one day), so you may need to install with `remotes::install_github("qsbase/qs")`.

To download and prepare input data and the reference for the workshop, change directories to the `setup/osteo` directory and run:

```sh
snakemake -j2
```

This will create both `../data/references/mm_mets_osteo_ref.rds` and a directory `../data/osteo/` with the following files:

```console
osteo
└── GSM8478586
    ├── normalized
    │   └── osteo_normalized_spe.rds
    └── spaceranger
        ├── filtered_feature_bc_matrix
        │   ├── barcodes.tsv.gz
        │   ├── features.tsv.gz
        │   └── matrix.mtx.gz
        └── spatial
            ├── aligned_fiducials.jpg
            ├── aligned_tissue_image.jpg
            ├── cytassist_image.tiff
            ├── detected_tissue_image.jpg
            ├── scalefactors_json.json
            ├── spatial_enrichment.csv
            ├── tissue_hires_image.png
            ├── tissue_lowres_image.png
            └── tissue_positions.csv
```

Then, `mm_mets_osteo_ref.rds` was copied to `/shared/data/training-modules/spatial/data/reference/`, and the `osteo` directory was copied to `/shared/data/training-modules/spatial/data/`. 