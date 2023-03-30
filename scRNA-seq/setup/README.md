# Single Cell RNA-seq training data setup"

This document describes how the training data is prepared for the single cell RNA-seq training data on the RStudio Server.

The preprocessing for these steps is organized as snakemake workflows.
Currently, these workflows do not include conda environments or docker, as those are not fully set up on this server.

## File locations

On the RStudio server the main location for the files needed for training data is `/shared/data/training-modules/scRNA-seq`.
The files are then organized by dataset.
For users, the required files are symlinked into the appropriate locations by the `link-data.sh` script


## Smart-Seq Data

The Smart-Seq data we are using comes from the study [GSE84465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84465), which corresponds to the SRA project SRP079058. 
This is glioblastoma data that was Fluorescence-Activated Cell sorted and processed by paired-end sequencing using Smart-seq2 protocol 
[(Darmanis *et al.* 2017)](https://pubmed.ncbi.nlm.nih.gov/29091775/).

The setup workflow for this data is in the subdirectory of this module.
The files it creates will be stored in the `/shared/data/training-modules/data` directory.
To run the workflow, use the terminal to change to the `setup/glioblastoma-darmanis` directory and then run the following command:

```sh
snakemake --cores 8 --keep-going --restart-times 2
```

This is a workflow that involves downloading, QC, and quasimapping for large number (>3500) of fastq files from ENA, but each fastq file should be deleted after processing (using the snakemake `temp()` directive), so the total amount of drive space used is not too bad. 

If more cores are available, the `--cores` argument can be increased. 

`--keep-going` and `--restart-times 2` are there because I did have issues with occasional download steps failing.
This should prevent the whole workflow from requiring repetition if an individual download fails (and hopefully it will succeed on one of the retries).

## Tabula Muris data

The 10X data we use is from the Tabula Muris dataset.

The data fro this script was originally downloaded as part of this workflow, but AWS has moved the bam files files to Glacier, so it is no longer quite so easy to download the raw data. 
For now, we have all we need already downloaded, so that portion of the workflow is commented out.
If we need to revisit these downloads in the future, updates will be required.

This workflow is simpler, as it doesn't include any Salmon processing and the sample list is predefined.
Files are stored in `/shared/data/training-modules/data/tabula-muris`
Again, the procedure is the same: change directories to the `setup/tabula-muris` directory and run:

```sh
snakemake --cores 8 
```

(Keep going options here did not seem necessary, but would be fine to include!)
