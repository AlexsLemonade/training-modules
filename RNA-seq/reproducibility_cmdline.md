### Making things more reproducible: the command line

> It came with all sorts of cryptic but powerful commands, which could be invoked by typing their names, and which I learned to use only gradually. (from "In the Beginning was the Command Line" by Neal Stephenson)

Most of us spend our time using graphical user interfaces.
On our computers we point and click with mice, and on tablets and phones we use our fingers to interact with apps.
However, it is very hard to describe exactly what was done after the fact.

In contrast with GUIs, the command line interfaces that we have started to use in this set of exercises provide reproducibility.
But each step can require extensive manual entry.
As we get to know our command line, we can get reproducibility and efficiency.
We can also script commands so that we can run them again anytime.

#### Shell scripts

When you are running commands by typing them into a command line, they are being interpreted and run by a "shell" (think of this as the program that converts your typed commands into steps that the computer runs).
We've been typing things one-by-one.
We can, instead, put many of these commands together into a single file.
This is something called a "shell script."
This is the preferred way to interact with a command line, because it provides a record of the exact set of commands that were run.

Let's create our first shell script.
In RStudio go to the file menu, select New File, Text File.
Name the script `rnaseq.sh` and put it in the kitematic folder.
Now, paste all of the commands that you ran to process the RNA-seq data sequentially into the `rnaseq.sh` file.
When you want to run the commands, change directories to the kitematic folder and type `sh rnaseq.sh`.
This will run the script, and hopefully each step will complete successfully!

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
mkdir QC

for filename in fastq/*.fastq.gz
do
  # remove leading path
  name=${filename##*/}
  # create directory for QC report
  mkdir QC/${name} -p
  # run fastqc
  fastqc ${filename} -o QC/${name}
done
```
In this code, `${filename}` gets replaced with the name of each individual file.
This lets us run the QC process over each file without typing it's name, making it easier to avoid the risk of typos disrupting our analyses.

We can also create an array of sequential numbers and then run salmon over every instance of an array:
```bash
declare -a arr=("SRR585570" "SRR585574")

for samp in "${arr[@]}"
do
  echo "Processing sample $samp"
  salmon quant -i index -l A \
        -1 fastq/${samp}_1.fastq.gz \
        -2 fastq/${samp}_2.fastq.gz \
        -p 4 -o quant/${samp} \
        --gcBias --seqBias --biasSpeedSamp 5
done
```

Feel free to use these techniques to improve your shell script as you write it.
