### Making things more reproducible: the command line

> It came with all sorts of cryptic but powerful commands, which could be invoked by typing their names, and which I learned to use only gradually. (from "In the Beginning was the Command Line" by Neal Stephenson)

Most of us spend our time using graphical user interfaces.
On our computers we point and click with mice, and on tablets and phones we use our fingers to interact with apps.
However, it is very hard to describe exactly what was done after the fact.

In contrast with GUIs, the command line interfaces that we have started to use in this set of exercises provide reproducibility.
But each step can require extensive manual entry.
As we get to know our command line, we can get reproducibility and efficiency.
We can also script commands so that we can run them again anytime.

### Download data

First, let's create a directory for our raw sequencing data (`.fastq.gz`) files.

```bash
mkdir data/fastq/gastric_cancer -p
```

We'll use the `-P` option for `wget` to download the files into our new directory.

#### SRR585570

```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR585/SRR585570/SRR585570_1.fastq.gz -P data/fastq/gastric_cancer
```
```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR585/SRR585570/SRR585570_2.fastq.gz -P data/fastq/gastric_cancer
```

#### SRR585574

```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR585/SRR585574/SRR585574_1.fastq.gz -P data/fastq/gastric_cancer
```
```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR585/SRR585574/SRR585574_2.fastq.gz -P data/fastq/gastric_cancer
```

#### Shell scripts

When you are running commands by typing them into a command line, they are being interpreted and run by a "shell" (think of this as the program that converts your typed commands into steps that the computer runs).
We've been typing things one-by-one.
We can, instead, put many of these commands together into a single file.
This is something called a "shell script."
This is the preferred way to interact with a command line, because it provides a record of the exact set of commands that were run.

Let's create our first shell script.
In RStudio go to the file menu, select New File, Text File.
Name the script `rnaseq.sh` and put it in the kitematic folder.
At the very top of the file put a [shebang](https://en.wikipedia.org/wiki/Shebang_(Unix)) like so:

```
#!/bin/bash
```

This tells the computer to execute this as a [Bash](https://en.wikipedia.org/wiki/Bash_(Unix_shell)) script.
Now, paste all of the commands that you ran to process the RNA-seq data sequentially into the `rnaseq.sh` file.

If you wanted to run the commands, you would type `bash rnaseq.sh`. 
(Ensure you are currently in the `kitematic` directory.)
This will run the script, and hopefully each step will complete successfully!
**We will not run this here in the interest of time, but you can experiment with this and the advanced topics below on your own time.**

#### Advanced topics

##### Discovering arguments

Many programs provide brief instructions - akin to a quick reference booklet.
For example, you've already run `wget` to download data.
Try running `wget --help`.
We can see that `wget` has many parameters.
Look specifically at `-i`: you can use this to download many files at once!
That could be much more convenient than typing each into the command line separately.

##### Loops

We wanted to process multiple samples at a time.
In the exercise, we typed the name of each file.
This took a lot of time, and if we had typos that really confused things.

Instead, we can write loops that work over many files.
To perform QC over all of the fastq files in a directory we can write:
```bash
mkdir <QC_DIRECTORY>

for filename in <FASTQ_DIRECTORY>/*.fastq.gz
do
  # remove leading path
  name=${filename##*/}
  # create directory for QC report
  mkdir <QC_DIRECTORY>/${name} -p
  # run fastqc
  fastqc ${filename} -o <QC_DIRECTORY>/${name}
done
```
where `QC_DIRECTORY is the directory to hold the QC reports and FASTQ_DIRECTORY is the directory that contains the `fastq` files.
The `<` and `>` indicate that this is to be replaced by you, the user, rather than run as is.

In this code, `${filename}` gets replaced with the name of each individual file.
This lets us run the QC process over each file without typing it's name, making it easier to avoid the risk of typos disrupting our analyses.

We can also create an array of sequential numbers and then run salmon over every instance of an array:
```bash
declare -a arr=("SRR585570" "SRR585574")

for samp in "${arr[@]}"
do
  echo "Processing sample $samp"
  salmon quant -i <INDEX_DIRECTORY> -l A \
        -1 <FASTQ_DIRECTORY>/${samp}_1.fastq.gz \
        -2 <FASTQ_DIRECTORY>/${samp}_2.fastq.gz \
        -o <QUANT_DIRECTORY>/${samp} \
        --gcBias --seqBias --biasSpeedSamp 5
done
```
where `INDEX_DIRECTORY` is the directory that contains the appropriate transcriptome index files, `FASTQ_DIRECTORY` is the directory that contains the `fastq` files, and `QUANT_DIRECTORY` is the directory where salmon quant output will go.

Feel free to use these techniques to improve your shell script as you write it.
