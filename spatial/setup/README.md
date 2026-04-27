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
Files for sample `SCPCS000190` were downloaded directly from the `ScPCA` portal, and we considered only files in `SCPCL000429_spatial/`.
The `json` and `html` files were removed, and the remaining files were organized as:
```markdown
└── SCPCS000190
    └── outs
        ├── filtered_feature_bc_matrix
        │   ├── barcodes.tsv.gz
        │   ├── features.tsv.gz
        │   └── matrix.mtx.gz
        ├── raw_feature_bc_matrix
        │   ├── barcodes.tsv.gz
        │   ├── features.tsv.gz
        │   └── matrix.mtx.gz
        └── spatial
            ├── aligned_fiducials.jpg
            ├── detected_tissue_image.jpg
            ├── scalefactors_json.json
            ├── tissue_hires_image.png
            ├── tissue_lowres_image.png
            └── tissue_positions_list.csv
```

The `SCPCS000190` directory was then placed in `/shared/data/training-modules/spatial/data/wilms-tumor/`.
No further preparation was needed.

### Ovarian carcinoma

This data comes from this 10x Genomics dataset: <https://www.10xgenomics.com/datasets/human-ovarian-cancer-11-mm-capture-area-ffpe-2-standard>.

This data, as well as its corresponding Visium probe set, were download from 10x as follows:

```sh
wget https://cf.10xgenomics.com/supp/spatial-exp/probeset/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv
wget -qO- https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma/CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma_spatial.tar.gz | tar -xz
wget -qO- https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma/CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma_filtered_feature_bc_matrix.tar.gz | tar -xz
```

The `Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv` file, the `spatial/` directory, and the `filtered_feature_bc_matrix/` directory were then placed in `/shared/data/training-modules/spatial/data/ovarian-carcinoma/`.
No further preparation was needed.


### Osteosarcoma

This data comes from this `GEO` record: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8478586>, which is associated with [Reinecke et al. (2025)](https://aacrjournals.org/clincancerres/article/31/2/414/751106/Aberrant-Activation-of-Wound-Healing-Programs).
This data needs to be filtered and normalized for input to workshop notebooks.

To download and prepare input data for the workshop, change directories to the `setup/osteo` directory and run:

```sh
snakemake -j2
```

This will place the downloaded files from GEO and the normalized SPE object into `/shared/data/training-modules/scRNA-seq-advanced/data/osteo/GSM8478586/`.



