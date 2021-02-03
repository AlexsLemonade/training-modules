# Introduction to single-cell RNA-seq

**CCDL 2020**

## Single Cell RNA-seq Technologies

For the purposes of this tutorial, we'll group single-cell technologies into two categories based on their capture methods and quantitative nature.
The pre-processing steps and the biases to look out for in post-processing vary based on technology and how the cells are sorted.

For more extensive background on single-cell experimental methods, Kiselev _et al._ have a very good [tutorial for scRNA-seq](https://scrnaseq-course.cog.sanger.ac.uk/website/index.html).

![**Overall view of full-length and tag-based workflows**](diagrams/overview_workflow.png)

### Full-length scRNA-seq  
*Example:* Smart-seq2 [(Picelli _et al._ 2014)](https://www.nature.com/articles/nprot.2014.006)    
Cells are physically separated into individual wells of a plate and are often also sorted by other means (e.g., Fluorescence Activated Cell Sorting).
Each cell is then sequenced individually and has its own fastq file.
(This will result in two fastq files if this is paired-end sequencing.)
The data preprocessing steps for these types of scRNA-seq data are similar to bulk RNA-seq methods.

#### Pros:  
- Can be paired-end sequencing which has less risk for 3' bias.  
- More complete coverage of transcripts, which may be better for transcript discovery purposes.   

#### Cons:  
- Is not very efficient (generally 96 cells per plate).  
- Takes much longer to run (days/weeks depending on sample size).
- Much more expensive.  

### Tag-based scRNA-seq  
*Example:* 10X Genomics [(Zheng _et al._ 2017)](https://www.ncbi.nlm.nih.gov/pubmed/28091601)    
Cells are separated by emulsion/droplets, and individual cells are given barcodes.
Everything is then sequenced.
These types of methods, because they are newer, are more likely to have Unique Molecular Identifiers (UMIs) [(Islam _et al._ 2014)](http://www.nature.com/doifinder/10.1038/nmeth.2772). 
These are additional barcodes added before amplification which allow the identification of PCR duplicates, allowing control for PCR amplification errors and biases.
Individual samples have two fastq files: one for the cell barcodes and another with the individual reads.

#### Pros:  
- Can run potentially millions of cells at once.   
- Much less computationally demanding.
- Won't take up all your computer's storage.  
- Much cheaper.  

#### Cons:  
- Sequencing is not bidirectional so data will likely have more intense 3' bias.  
- Coverage of these technologies is generally not as deep.  

## Resources:

- Hemberg Lab [scRNA-seq training course](https://scrnaseq-course.cog.sanger.ac.uk/website/index.html)

- [ASAP: Automated Single-cell Analysis Pipeline](https://asap.epfl.ch/) is a web server that allows you to process scRNA-seq data. ([Gardeux _et al._ 2017.](https://doi.org/10.1093/bioinformatics/btx337 ))

- Smith. [_Unique Molecular Identifiers – the problem, the solution and the proof_](https://cgatoxford.wordpress.com/2015/08/14/unique-molecular-identifiers-the-problem-the-solution-and-the-proof/) provides an excellent background on UMIs.

### Literature on the comparisons and explanations of scRNA-seq  technologies:*
- [Amezquita _et al._ (2019)](https://www.biorxiv.org/content/10.1101/590562v1)  
- [Lafzi _et al._ (2019)](https://doi.org/10.1038/s41596-018-0073-y)  
- [Luecken & Theis (2019)](https://doi.org/10.15252/msb.20188746)  
- [AlJanahi _et al._ (2018)](https://doi.org/10.1016/j.omtm.2018.07.003)  
- [Baran-Gale _et al._ (2018)](https://doi.org/10.1093/bfgp/elx035)  
- [Ziegenhain _et al._ (2018)](http://dx.doi.org/10.1016/j.molcel.2017.01.023)  
- [Zhang _et al._ (2018)](https://doi.org/10.1016/j.molcel.2018.10.020)  
- [Angerer _et al._ (2017)](http://dx.doi.org/10.1016/j.coisb.2017.07.004)  
