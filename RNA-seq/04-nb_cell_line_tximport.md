## Harenza, et al. neuroblastoma cell line data

In this section, we'll be working with RNA-seq data from neuroblastoma (NB) cell lines from
[Harenza, et al. _Scientific Data._ 2017.](https://doi.org/10.1038/sdata.2017.33)

The course directors have already processed the raw data using `salmon quant` and the `quant.sf` files for each sample can be found in `data/NB-cell/salmon_quant/<SAMPLE>`.

![](diagrams/rna-seq_5.png)

In the gastric cancer example, we imported Salmon-processed data with `tximeta` to then use with `DESeq2`.
We will also use `DESeq2` for these analyses, specifically for differential expression analysis.

In order to prepare the NB cell line data for differential expression analysis, we will modify the gastric cancer tximeta notebook (`02-gastric_cancer_tximeta-live.Rmd`) and save this new notebook as `nb_cell_tximeta`:

* To create a new notebook, select `File` > `New File` > `R Notebook`.
The new notebook should appear in your Source Pane in RStudio.
Save the new notebook, using Ctrl+S (Cmd+S on Mac) or `File` > `Save`, in the `training-modules/RNA-seq` directory with the name `nb_cell_line_tximeta.Rmd`.
You can add a new chunk by clicking the *Insert Chunk* button on the toolbar or by pressing *Cmd+Option+I*.

* Alter the code from `02-gastric_cancer_tximeta-live.Rmd` to use the NB cell line data.
The `quant.sf` files for each sample can be found in `data/NB-cell/salmon_quant/<SAMPLE>`.

* Save the `tximeta` output as `data/NB-cell/tximeta/NB-cell_tximeta.RDS`. Note that `data/NB-cell/tximeta/` is a new directory.
