## Machine learning training data setup

This document provides some background for setting up the data we use in the guided notebooks of the machine learning module. 
We use the [OpenPBTA](https://github.com/AlexsLemonade/OpenPBTA-analysis) stranded RNA-seq data for this module.

#### Data download

We download the `release-v16-20200320` version of the OpenPBTA data to `machine-learning/data/open-pbta/download` (relative to the root of the repository) with:

```
bash 00-data-download.sh
```

The release version can be altered with:

```
OPENPBTA_RELEASE=<different release> bash 00-data-download.sh
```

But the directory where these data are downloaded is hard coded.

#### Data cleaning

We perform variance-stabilizing transformation on the `pbta-gene-counts-rsem-expected_count.stranded.rds` data via DESeq2 in the `01-transform-rnaseq` notebook.
We also create a version of `pbta-histologies.tsv` that is filtered to only include stranded RNA-seq samples.
(Columns that are entirely missing values after filtering are removed.)

These data are saved in `machine-learning/data/open-pbta/processed`, which is the data directory used throughout the guided notebooks.
