# scRNA-Seq Training Module

This CCDL-designed module covers the analysis of single-cell RNA-seq data using [Salmon's Alevin](https://salmon.readthedocs.io/en/latest/alevin.html) and [scater/scran](https://bioconductor.org/packages/devel/bioc/vignettes/scran/inst/doc/scran.html) packages.
It is designed to be taught in approximately 2.5 hours.
It depends on knowledge gained in the [Intro to R](https://github.com/AlexsLemonade/training-modules/tree/master/intro-to-R-tidyverse) module and analyses are performed within a [Docker container](https://github.com/AlexsLemonade/training-modules/tree/master/docker-install).
It covers normalization and dimension reduction methods that can be used for both tag-based and full-length single-cell data, as well as the quantification of tag-based scRNA-seq data.

The notebooks that comprise this module are:
* [Normalization using scater/scran](https://github.com/AlexsLemonade/training-modules/blob/master/scRNA-seq/01-normalizing_scRNA-seq.nb.html)
* [Quantification of tag-based data](https://alexslemonade.github.io/training-modules/scRNA-seq/02-tag-based_pre-processing_scRNA-seq.md)
* [Dimension reduction](https://alexslemonade.github.io/training-modules/scRNA-seq/03-dimension_reduction_scRNA-seq.nb.html)
* [Additional exercises](https://github.com/AlexsLemonade/training-modules/blob/master/scRNA-seq/04-scrnaseq_exercise.Rmd)
