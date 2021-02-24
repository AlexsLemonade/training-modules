# Data setup for RNA-seq module

This directory contains a snakemake workflow and two configuration files to download and preprocess the required data files for the RNA-seq module notebooks.
The config files are:
- `config-GC.yaml` for the gastric cancer data files used in the majority of the notebooks.
- `config-NB.yaml` for the neuroblastoma cell line data used in the guided exercise and beyond.
- `config-leukemia.yaml` for the mouse AML model used in the exercise notebooks
- `config-zebrafish.yaml` for the zebrafish data used by the tximeta "other species" example.

The full setup requires running the workflow twice, once with each config file, which can be done with the `setup.sh` script.

The files are saved by default to `/shared/data/training-modules/RNA-seq`, though this can be modified by adjusting the `base_dir` setting in the config files.
As the raw fastq data is not used in any notebooks, those files are set as `temp` and deleted by default after processing with `salmon` and `fastqc`.
To keep those files, `snakemake` should be run with the `--notemp` flag, as is done for the gastric cancer workflow in `setup.sh`.
