# Contributing and Development Guidelines

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Instruction notebook content](#instruction-notebook-content)
  - [Learning Objectives](#learning-objectives)
  - [Style guide](#style-guide)
  - [Code Chunks](#code-chunks)
  - [References](#references)
  - [Notebook data dependencies](#notebook-data-dependencies)
- [Data file management](#data-file-management)
  - [Setup directory](#setup-directory)
  - [Server shared data folder](#server-shared-data-folder)
  - [Linking shared files](#linking-shared-files)
  - [Files stored on S3](#files-stored-on-s3)
- [Development with `renv`](#development-with-renv)
  - [Typical development workflow](#typical-development-workflow)
  - [Steps for creating renv.lock only changes pull requests](#steps-for-creating-renvlock-only-changes-pull-requests)
    - [Multiple renv.lock changes from multiple branches](#multiple-renvlock-changes-from-multiple-branches)
  - [How we use `renv` with Docker](#how-we-use-renv-with-docker)
- [Docker Image](#docker-image)
  - [Developing within the Docker container](#developing-within-the-docker-container)
  - [Testing Docker image builds via GitHub Actions](#testing-docker-image-builds-via-github-actions)
  - [Pushing to Docker Hub via GitHub Actions](#pushing-to-docker-hub-via-github-actions)
- [Automated Testing \& Rendering](#automated-testing--rendering)
  - [Spell check](#spell-check)
  - [Rendering Test](#rendering-test)
  - [Generation of live notebooks and rendering](#generation-of-live-notebooks-and-rendering)
- [Cheatsheets](#cheatsheets)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Instruction notebook content

Live instruction notebooks should contain roughly 90 minutes of live teaching material.
They should additional contain enough context to serve as stand-alone tutorials for participants to use after workshops.

### Learning Objectives

Each notebook should begin with a "Learning objectives" section.
This section contains a compact summary of the goals for the notebook, in the form of a bulleted list of objectives, preceded by the following header:

```
## Objectives

This notebook will demonstrate how to:

```
Each element of the objective list should be a specific concept and/or skill we expect learners to come away from the notebook having gained.
Phrase the objective with an "action verb" followed by a description of the skill.
The Carpentries [Curriculum Development Handbook section on learning objectives](https://cdh.carpentries.org/developing-content.html#learning-objectives) has some useful discussion about how to write effective learning objectives for courses like ours.

For most notebooks, 3-5 objectives should be sufficient.

The objectives list should be followed by a horizontal rule for visual distinction:

```
---
```

### Style guide

More to come, but for now, we should generally follow the style conventions we established in [`refinebio-examples`](https://github.com/AlexsLemonade/refinebio-examples/blob/staging/CONTRIBUTING.md#notebook-text) for text.

Some style considerations include:

* Use `::` syntax with functions (e.g. `dplyr::filter()`) even if the library is loaded, to reduce ambiguity
* When relevant, include explicit function arguments even when using the defaults (see [here for an example](https://github.com/AlexsLemonade/training-modules/blob/d2c63104c72a30273b1dd3016b03f01ffcc2e6b4/scRNA-seq/06-celltype_annotation.Rmd#L516-L517))
* Incorporate code comments that explain both what and why code is performing a certain task
* All notebooks should [conclude in a `sessionInfo()` chunk](https://github.com/AlexsLemonade/training-modules/blob/d2c63104c72a30273b1dd3016b03f01ffcc2e6b4/scRNA-seq/06-celltype_annotation.Rmd#L785-L789) to report out the compute environment

### Code Chunks

All code chunks should be named.
This is done to ease orientation and navigation during training.

Code chunks that will be blank at the start of a training session (for live coding) should be tagged with the `live = TRUE` argument.
These chunks will be stripped of code (but not full line comments) when the `-live.Rmd` version of the notebook is created by the `make-live.R` script.
This script is not usually run independently; it will usually run as part of the [Make Live Notebooks](https://github.com/AlexsLemonade/training-modules/actions/workflows/make-live.yml) GitHub Action via a[`workflow_dispatch`](https://docs.github.com/en/actions/managing-workflow-runs/manually-running-a-workflow#running-a-workflow-on-github).
That action will file a pull request with changes to any `-live.Rmd` notebooks and, optionally, rendered versions of the notebooks.


### References

We try to maintain good scholarship, citing our sources!
Often our sources are vignettes or web pages, for which we usually link directly to a web page.
When possible we should include the author of the vignette in the text to give proper credit.

In the case of journal or preprint publications, we will follow the following citation conventions:

- `([Author Date](url))` with no reference to journal except in the url (the url should be a DOI if possible).
- If there are two authors: `([Author1 and Author2 YEAR](url))`
- If 3 or more: `([Author1 _et al._ YEAR](url))`
- Citations in text: `As shown by [Author1 _et al._ (YEAR)](url)`.

### Notebook data dependencies

Each instruction notebook should run from start to finish when run after the previous notebook(s) in its module have been run.
This means that modules can have dependencies on each other, but those dependencies should be satisfied when run in order.
In general, input files that are not present in this repository will be linked as part of setting up the repository or user folder, as described in the [Data file management](#data-file-management) section of this document.

An exception is when a notebook relies on the completion of tasks run outside the notebook (salmon mapping, for example).
In that case, a Data Lab member will upload any extra required files to S3 via the [`syncup-s3.sh` script](#files-stored-on-s3) so that the notebook can run to completion during automated testing.


## Data file management

### Setup directory

Each module should contain a `setup` directory that includes instructions and scripts to download and prepare any input data files used by the notebooks in the module.
For scripts designed to be run on the RStudio server (`rstudio.ccdatalab.org`) these input files will usually be placed in the `/shared/data/training-modules` directory and then linked for use, as described below.

### Server shared data folder

When trainings are run from the RStudio server (`rstudio.ccdatalab.org`), we store large input data files in the `/shared/data/training-modules` directory.
This directory is the implicit "point of truth" for modules

The organization within this directory should mirror the arrangement of the repository, so that we can easily link files and folders from this shared directory to mirrored paths within a clone of this repository.

### Linking shared files

Linking files from the shared directory to a cloned repository is done with the `scripts/link-data.sh` bash script, so this script should be kept up to date as any needed files are added to the `/shared/data/training-modules` directory.
When possible, link to enclosing directories rather than individual files to keep links simpler and allow users to browse a realistic directory context, but see an important caveat below.

Because this script is also used to set up directories for training, the links should not include _all_ files needed for _every_ notebook:
- Files that are created during a training session should not  be included in this script.
- Directories that users will need to write to should not be links, or the user will not be able to write their own files.

### Files stored on S3

Note that only Data Lab members have access to S3.

To facilitate automated testing of training notebooks, all needed input files for training notebooks should be placed in the `ccdl-training-data` bucket on S3 and made publicly accessible.
This is facilitated by the `scripts/syncup-s3.sh` bash script, which includes the needed commands for upload/sync, and should include all directories and files needed to run the training notebooks.
You will need to set up your AWS credentials with `aws configure` before running the `syncup-s3.sh` script
As input files are added or change, those changes should be reflected in updates to `syncup-s3.sh`

## Development with `renv`

We use [`renv`](https://rstudio.github.io/renv/index.html) to manage R package dependencies for this project.
Using `renv` allows us to keep R packages in sync during _development_ in multiple scenarios – on the RStudio Server, using the project Docker image, or even locally – and generates a lockfile (`renv.lock`) that we can use to install packages on the RStudio Server for participants when it's time for a workshop.

### Typical development workflow

**For `renv` to work as intended, you will need to work within the `training-modules` project.**
Be careful not to open any module specific `.Rproj` during development, as it will disrupt the `renv` environment!

The steps for development are:

1. Open up `training-modules.Rproj`.
2. If your library isn't set up or in sync with the lockfile, you'll be prompted to run `renv::restore()`. This will happen if you first clone the project or haven't been working within the Docker container in a bit.
3. Develop as normal.
4. Run `renv::snapshot()` at the end of your session to capture the additional dependencies. *Be careful if `renv::snapshot()` suggests removing packages!* If there have been additions to `renv.lock` in another branch while you were working, you may need to run `renv::restore()` again before `renv::snapshot()`.
5. If there are dependencies you might want that are not captured automatically by `renv::snapshot()` (this may happen if a package is "recommended" by another, but not required), add them to `components/dependencies.R` with a call to `library()` and an explanatory comment. Then rerun `renv::snapshot()`
6. Commit any changes to `renv.lock` and `dependencies.R`.
7. File a pull request with _only_ the changes to the `renv.lock` and `dependencies.R`  before other changes (see next section for tips on creating this PR).
This is necessary because the automated render testing we do will fail if the Docker image has not been updated (see '[Pushing to Docker Hub via GitHub Actions](#pushing-to-docker-hub-via-github-actions)' below).

Note that when you open up the `training-modules.Rproj`, the `.Rprofile` file makes it such that the `renv` library is loaded and the repositories in the `renv.lock` file will be set with `options(repos = ...)`.

### Steps for creating renv.lock only changes pull requests

For most cases you can create your `renv.lock-only` changes PR by following these steps:

1. Create your `renv.lock-only` branch from the latest `master` branch.
2. In your `renv.lock-only` branch, checkout the renv.lock file from your development branch (where you were generally doing steps 1-6 from the previous section) using `git checkout devbranch renv.lock`.
3. Commit the renv.lock changes you just checked out.
4. Push the changes and file your renv.lock only PR.

Note that the `renv::snapshot()` command will skip adding a package to renv.lock if it isn't used in a notebook or script.

#### Multiple renv.lock changes from multiple branches

If there are changes happening on multiple branches that require renv.lock changes, you may need to follow a slightly different version of steps:

1. Create your `renv.lock-only` branch from the latest `master` branch.
2. Run `renv::restore()`.
3. Install the packages needed on both branches (`install.packages()` or etc).
4. Add those packages to `components/dependencies.R`.
5. Run `renv::snapshot()`.
6. Only commit the `renv.lock` changes to your branch.
7. Push the changes and file your renv.lock only PR.

### How we use `renv` with Docker

We use the `renv.lock` file in this repository to install R dependencies on the image per the [`renv` documentation for creating Docker images](https://rstudio.github.io/renv/articles/docker.html#creating-docker-images-with-renv-1).
Specifically, we copy the `renv.lock` file and run `renv::restore()`, which installs the packages in the lockfile into the Docker image's system library.

In practice, this means that you will not need to add individual R packages to the `Dockerfile`, but you may have to add _system dependencies_ (e.g., via `apt-get`) required for the installation of those R packages.

## Docker Image

### Developing within the Docker container

To use the Docker image for development, pull from Docker Hub with:

```
docker pull ccdl/training_rstudio:edge
```

To run the container and mount a local volume, use the following from the root of this repository:

```
docker run \
  --mount type=bind,target=/home/rstudio/training-modules,source=$PWD \
  -e PASSWORD=<PASSWORD> \
  -p 8787:8787 \
  ccdl/training_rstudio:edge
```

Replacing `<PASSWORD>` with the password of your choice.
You can then navigate to `localhost:8787` in your browser and login with username `rstudio` and the password you just set via `docker run`.

To work on the project, you should then open `training-modules/training-modules.Rproj` on the RStudio server, then proceed as in the [typical development workflow](#typical-development-workflow).

### Testing Docker image builds via GitHub Actions

When a pull request changes either the `Dockerfile` or `renv.lock`, a GitHub Actions workflow (`build-docker.yml`) will be run to test that the image will successfully build.

### Pushing to Docker Hub via GitHub Actions

When a pull request is merged into `master`, the `build-docker.yml` GitHub Actions workflow will be triggered.
The project Docker image will be rebuilt and pushed as `ccdl/training_rstudio:edge`.

When a new Git tag is created, the `build-docker.yml` GitHub Actions workflow will also be triggered, and will push a version of the image to Docker Hub as `ccdl/training_rstudio:<tag>`.

The most recent tagged version of the image will also be tagged as as `ccdl/training_rstudio:latest`.

## Automated Testing & Rendering

### Spell check

We perform spell checking for every pull request to `master` as part of a GitHub Actions workflow (`spell-check.yml`); it is designed to fail if any spelling errors are detected.
This workflow uses the [`AlexsLemonade/spellcheck` GitHub Action](https://github.com/alexslemonade/spellcheck).
You can see what errors are detected by downloading the `spell_check_errors` artifact from the completed workflow.

We maintain a custom at `components/dictionary.txt` where words that should not be flagged as typos can be added.

### Rendering Test

Every pull request to `master` that changes `.Rmd` files (or one of the rendering scripts) will be tested via a GitHub Actions workflow (`render-rmds.yml`) to ensure that the `.Rmd` files can be run successfully.
This action first downloads input files for the notebooks from S3, so if there are changes to the input files, these should be made first, with associated changes as needed to `syncup-s3.sh` (see [Files stored on S3](#files-stored-on-s3)).


### Generation of live notebooks and rendering

After a pull request with changes to notebook files has been merged to master, we use the `make-live.yml` workflow to render current versions of the notebooks to html and to make the `-live.Rmd` versions of the files for training sessions.
This workflow then files a PR to `master` with the rendered and live files.
`make-live.yml` is currently manually triggered, but will likely change to running automatically on each PR with changes to notebook files in the near future.


## Cheatsheets

Training modules have corresponding cheatsheets in `module-cheatsheets`.
When choosing documentation links to incorporate in cheatsheets, we prefer to use [`https://rdrr.io/`](https://rdrr.io/) when possible for Base R and Bioconductor, and we prefer to use [`https://www.tidyverse.org/`](https://www.tidyverse.org/) for `tidyverse` functions.

Cheatsheets are written in plain markdown and are converted to a shareable PDF format using the `Node.js` package [`mdpdf`](https://github.com/BlueHatbRit/mdpdf), with the default PDF style.

To render these packages, you will therefore have to first install `npm`, the `Node.js` package manager.
You can install `Node.js` and `npm` using Homebrew with `brew install node`, or you can install into a Conda environment with `conda install nodejs`
With `Node.js` installed, you can install the `mdpdf` package into the default `Node.js` library with:

```
npm install -g mdpdf
```

In addition, cheatsheet table of contents are created with the `Node.js` library [`doctoc`](https://github.com/thlorenz/doctoc), which can similarly be globally installed with:

```
npm install -g doctoc
```

To re-render a cheatsheet to PDF after making desired changes in its markdown, take the following steps:

* Navigate in terminal to the `module-cheatsheets` directory
* Run `doctoc` on the file to update its table of contents: `doctoc cheatsheet-file.md`
  * Alternatively, you can use a tool such as the VS Code Extension [Markdown All in One](https://marketplace.visualstudio.com/items?itemName=yzhang.markdown-all-in-one).
  Either way, please ensure that the table of contents is not duplicated.
* Convert the markdown file to an updated PDF version: `mdpdf cheatsheet-file.md`
