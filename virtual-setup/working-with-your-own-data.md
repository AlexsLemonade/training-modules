## Working with your own data

The goal of our workshop is to equip you to do initial analyses with your own data!
This guide will take you through how to get your data onto our RStudio server so you can begin analyzing your own data!

### Things to know before uploading your data:
- If you are uploading data from human samples, **please be sure that you are doing so in a manner that is consistent with participant consents and your institution’s rules**.

- Initially, we have equipped you with **50 GB of space** (if the data you would like to upload is larger than this, please consult on of the CCDL team members through Slack for assistance).

- If you don't have your own data that you are looking to analyze, but would like real transcriptomic datasets to practice with, see our recommended list here:
TODO: Add info about recommended datasets.

- As always, please Slack one of the CCDL team members if you need help (that is what we are here for!).

### Upload data that is online (from a url)

If you are retrieving your data from online, perhaps from a publicly available repository, we encourage you to use the terminal command `wget`.

**Step 1)** Go to the Terminal tab in your RStudio session.

**Step 2)** Set up your `wget` command in a script (or notebook) using this template.

The most simple `wget` command just needs the URL to pull the file from.

**Template:**
```
wget <URL>
```

**Specific example:** Here's an example of us downloading a file from ArrayExpress
```
wget https://www.ebi.ac.uk/arrayexpress/files/E-GEOD-67851/E-GEOD-67851.processed.1.zip
```

By default, the file will be saved to the current directory and the file name it had from its origin (so with the above example `E-GEOD-67851.processed.1.zip`).

Likely you will want to be more specific about where you are saving the file to and what you are calling it.
For that, we can use the `-O`, or `output` option with our `wget` command and specify a file path.  

**Template:**
```
wget -O <FILE_PATH_TO_SAVE_TO> <URL>
```

**Specific example:** Here's an example where we will download that same array express file, but instead save it to the `data` folder and call it `some_array_data.zip`.
(Best to keep the file extension consistent to avoid troubles!)

```
wget -O data/some_array_data.zip https://www.ebi.ac.uk/arrayexpress/files/E-GEOD-67851/E-GEOD-67851.processed.1.zip
```

`-O` is one of many `wget` command options.
To see the complete list of `wget` options, use the command: `wget -h` in bash.

[See more `wget` examples](https://www.tecmint.com/10-wget-command-examples-in-linux/)

### Upload data from your own computer

If the data you want to use is stored locally on your computer, here's how we recommend uploading it to the RStudio Server.

**Step 1)** We recommend you compress your data folder into a single zip file.

For most operating systems, you can right-click on your data folder, and choose `Compress` to zip up your files
- [Windows zipping](https://support.microsoft.com/en-us/help/14200/windows-compress-uncompress-zip-files)
- [Mac zipping](https://www.imore.com/how-compress-file-your-mac)

For reference, here's how you [compress files from the command line](https://coolestguidesontheplanet.com/how-to-compress-and-uncompress-files-and-folders-in-os-x-lion-10-7-using-terminal/).

**Step 2)** Once your data is compressed to a single file, [navigate to your RStudio session](./rstudio-login.md).

**Step 3)** Use the `Upload button` to choose your compressed data folder.

This button is in the lower right panel of your RStudio session:
![Upload button](screenshots/upload-button.png)

A mini screen will pop up asking you to choose the file you want to upload:

![Choose file](screenshots/upload-choose-file.png)

Choose your compressed data file, and click `OK`.
This may take some time, particularly if you have a large dataset.  

When the server is finished uploading your data, you should see your file in your `home` directory!

### Reading the data into your R environment.

This step is very dependent on the format of your data and what you are planning to do with it!
If your file is a TSV, CSV, or RDS file, follow the examples in the [`intro-to-R-tidyverse/intro-to-tidyverse.Rmd` notebook](intro-to-R-tidyverse/03-intro_to_tidyverse.Rmd).

We do NOT recommend clicking on the file in the RStudio panel to load it into your R environment, this typically won't work for anything that's not a very small file.

Plus, for reproducibility purposes, you should write the data reading step into your notebook analyses!

## Downloading files

- Where in their folders do we want our participants to be storing their data?
