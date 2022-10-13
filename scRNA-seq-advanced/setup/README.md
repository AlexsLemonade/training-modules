# Advanced Single Cell RNA-seq training data setup"

This document describes how the training data is prepared for the advanced single cell RNA-seq training data on the RStudio Server.

As new training data is added, notebooks, scripts, and/or workflows should be added to this directory to describe the steps required to recreate the required training input files. 
These files should include source locations for downloading raw data as appropriate, and all processing steps required to prepare data for use.

## File locations

On the RStudio server the main location for the files needed for this training module will be `/shared/data/training-modules/scRNA-seq-advanced/`.
The files are then organized by dataset.

After setup, the symlinks should be established from `/shared/data/training-modules/scRNA-seq-advanced/` to `training-modules/scRNA-seq-advanced/data` as appropriate.
These links will usually be created using the script at `training-modules/scripts/link-data.sh`, so directories and files should be added to the `link_locs` array in that script.
If no data will be written to a data directory, linking to that directory will be sufficient, but in some cases links to individual files may be required.


