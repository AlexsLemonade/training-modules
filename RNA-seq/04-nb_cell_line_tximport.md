## Harenza, et al. neuroblastoma cell line data

In this section, we'll be working with RNA-seq data from neuroblastoma (NB) cell lines from
[Harenza, et al. _Scientific Data._ 2017.](https://doi.org/10.1038/sdata.2017.33)

The course directors have already processed the raw data using `salmon quant` and the `quant.sf` files for each sample can be found in `data/quant/NB_cell_line/<SAMPLE>`.

In the gastric cancer example, we imported Salmon-processed data with `tximport` to then use with `DESeq2`. 
We will also use `DESeq2` for these analyses, specifically for differential expression analysis.

In order to prepare the NB cell line data for differential expression analysis, modify the gastric cancer tximport notebook (`02-gastric_cancer_tximport.Rmd`) and save this new notebook as `nb_cell_line_tximport`.

Save the `tximport` output as `data/tximport/NB_cell_line/NB_cell_line_tximport.RDS`. Note that `data/tximport/NB_cell_line` is a new directory.
