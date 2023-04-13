# Converting to/from `SCE` objects

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [`SCE` vs `Seurat` objects](#sce-vs-seurat-objects)
- [Converting from `Seurat` to `SCE`](#converting-from-seurat-to-sce)
- [Converting from `SCE` to `Seurat`](#converting-from-sce-to-seurat)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


**This cheatsheet explains how you can convert single-cell experiment data in R between `SingleCellExperiment` (`SCE`) and `Seurat` formats.**


## `SCE` vs `Seurat` objects

When converting between `Seurat` and `SCE` objects, it's helpful to know how the different object types store and refer to similar information.
The table below shows different aspects of single-cell objects and how to access the associated data, assuming the default names for each type of single-cell object.
There are several differences between `Seurat` and `SCE` objects that are useful to be aware of when converting them.
Importantly, the term `"assay"` refers to different things in `SCE` vs. `Seurat` objects:
  - In an `SCE` object, an `assay` is a matrix of counts, with default names `"counts"` for raw counts and `"logcounts"` for normalized counts.
  - In a `Seurat` object, an `assay` instead refers to an _experiment_. The default `Seurat` assay is called `"RNA"`, and it is analogous to the "main experiment" in an `SCE` object, which is not given a particular name.
  - The `Seurat` count matrices are stored within a given assay (experiment) and have default names of `"counts"` for raw counts and `"data"` for normalized counts.

In addition, by default, `SCE` reduced dimension names are capitalized (e.g., `"PCA"`), and `Seurat` reduced dimension names are in lower case (e.g., `"pca"`).

Always bear in mind that your object(s) may be named differently from the defaults as described here!

| Data aspect | `SCE` | `Seurat` |
|------------|---------|---------|
| Raw counts matrix | `counts(sce_object)` | `seurat_obj[["RNA"]]@counts` |
| Normalized counts matrix | `logcounts(sce_object)` | `seurat_obj[["RNA"]]@data` |
| Reduced dimension: PCA matrix | `reducedDim(sce_object, "PCA)` | `seurat_obj$pca@cell.embeddings` |
| Reduced dimension: UMAP matrix | `reducedDim(sce_object, "UMAP)` | `seurat_obj$umap@cell.embeddings` |
| Cell-level metadata | `colData(sce_object)` | `seurat_obj@meta.data` |
| Feature (gene)-level metadata | `rowData(sce_object)` | `seurat_obj[["RNA"]]@meta.features`|
| Miscellaneous additional metadata | `metadata(sce_object)` | `seurat_obj@misc`|

We provide some code examples below for these conversions below.
For all code examples below, it is assumed that the `SingleCellExperiment` library has been loaded into your R environment:

```r
library(SingleCellExperiment)
```


## Converting from `Seurat` to `SCE`

The following example code assumes you have a `Seurat` object called `seurat_obj`.

```r
# Convert Seurat object to SCE object
sce_object <- Seurat::as.SingleCellExperiment(seurat_obj)
```

By default, all assays (experiments) present in the `Seurat` object will be ported into the new `SCE` object.
Recall, in `Seurat`, an assay refers to an _experiment_ which may be associated with multiple count matrices.
To only specify that certain assays are retained, you can optionally provide the argument `assay` with _`Seurat` assay names_ to retain in the `SCE` object, for example:


```r
# Convert Seurat object to SCE object, retaining only the 'RNA' experiment (assay)
sce_object <- Seurat::as.SingleCellExperiment(seurat_obj, assay = "RNA")
```

Specifying `assay` is mostly useful if there are alternative experiments, for example from CITE-Seq data, present in the `Seurat` object that you do not want to retain during `SCE` conversion.

## Converting from `SCE` to `Seurat`

The following example code assumes you are starting with an `SCE` object called `sce_object`.

The function `Seurat::as.Seurat()` can be used to convert an `SCE` object into a `Seurat` object and takes the following arguments:

- The `SCE` object to convert
- Optional named arguments with the following defaults:
  - `counts = "counts"` specifies that the `SCE` object contains a `"counts"` assay of normalized counts that should be included during conversion.
    - If there is no `"counts"` assay in the SCE object, set this argument as `counts = NULL` or rename accordingly, e.g. `counts = "whatever_assay_name_you_are_using"`.
  - `data = "logcounts"` specifies that the `SCE` object contains a `"logcounts"` assay of normalized counts that should be included during conversion.
    - If there is no `"logcounts"` assay in the SCE object, set this argument as `data = NULL` or rename accordingly, e.g. `data = "whatever_assay_name_you_are_using"`.
  - `assay = NULL` specifies that, by default, all assays (experiments) will be converted. If there are multiple assays and you wish to only convert, for example, the `"RNA"` assay, set this argument as `assay = "RNA"`.
  - `project = "SingleCellExperiment"` specifies that the `Seurat` object being created will have this associated project name. You can override this with any string of interest, e.g. `project = "sample_XYZ"`.


```r
# Convert SCE object to Seurat object, assuming both
#  `counts` and  `logcounts` assays are present
seurat_object <- Seurat::as.Seurat(sce_object)

# Convert SCE object to Seurat object, where the SCE object
#  contains a `counts` but not a `logcounts` assay
seurat_object <- Seurat::as.Seurat(sce_object, data = NULL)
```
