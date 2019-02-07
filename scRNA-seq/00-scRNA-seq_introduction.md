# Pre-processing single-cell RNA-seq data

**CCDL 2019**

#### In this section, we will be running through the basics of pre-processing single-cell RNA-seq data.

For the purposes of this tutorial, we'll summarize single-cell technologies as
being one of two groups, based on their capture methods and quantitative nature.
Depending on how the cells are sorted and what technology is used, the pre-processing steps are a bit different and the biases to look out for in post-processing also vary.

For a more extensive run through on single-cell technologies,
Kiselev et al have very [good tutorial for scRNA-seq](https://hemberg-lab.github.io/scRNA.seq.course/introduction-to-single-cell-rna-seq.html#experimental-methods).

### 1) Non-Tag-Based scRNA-seq  
*Example:* Smart-seq2 [(Picelli et al, 2014)](https://www.nature.com/articles/nprot.2014.006)   
Cells are physically separated generally into individual wells in a plate and
often also sorted by other means (eg. Fluorescence Activated Cell Sorting).
Each cell is then sequenced individual and has it's own fastq file (this will be two fastq files if this is paired-end sequencing.)
The data pre-processing steps for these types of scRNA-seq data is/can be more similar to
bulk RNA-seq methods.

#### Pros:  
- Can be paired end sequencing which has less risk for 3' bias.  
- More complete coverage of transcripts, which may be better for transcript
discovery purposes.   

#### Cons:  
- Is not very efficient (generally 96 cells per plate).  
- Takes much longer to run (Days/weeks depending on sample size).
- Is a lot more expensive.  

### 2) Tag-based scRNA-seq  
*Example:* 10X Genomics [(Zheng et al, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/28091601)  
Cells are separated by emulsion/droplets, and individual cells are given barcodes.
Everything is then sequenced.
These types of methods, because they are newer, are more likely to have
[Unique Molecular Identifiers (UMIs)](http://www.nature.com/doifinder/10.1038/nmeth.2772)
which allow you to better control for PCR amplification errors and biases.
Individual samples have two fastq files: one for the cell barcodes
and another with the individual reads.

#### Pros:  
- Can run potentially millions of cells at once.   
- A lot quicker computing wise.  
- Won't take up all your computer's storage.  
- A lot cheaper.  

#### Cons:  
- Sequencing is not bidirectional so data will likely have more intense 3' bias.  
- Coverage of these technologies generally is not as deep.  

*More sources on the comparisons and explanations of these technologies:*   
- [Zhang et al, 2018](https://doi.org/10.1016/j.molcel.2018.10.020)  
- [AlJanahi et al, 2018](https://doi.org/10.1016/j.omtm.2018.07.003)  
- [Angerer et al, 2017](http://dx.doi.org/10.1016/j.coisb.2017.07.004)  
- [Baran-Gale et al, 2018](https://doi.org/10.1093/bfgp/elx035)  
