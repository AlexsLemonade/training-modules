# Contributing and Development Guidelines

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Development with `renv`](#development-with-renv)
  - [Typical development workflow](#typical-development-workflow)
  - [How we use `renv` with Docker](#how-we-use-renv-with-docker)
- [Docker Image](#docker-image)
  - [Developing within the Docker container](#developing-within-the-docker-container)
  - [Testing Docker image builds via GitHub Actions](#testing-docker-image-builds-via-github-actions)
  - [Pushing to Dockerhub via GitHub Actions](#pushing-to-dockerhub-via-github-actions)
- [Spell check](#spell-check)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


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
4. Run `renv::snapshot()` at the end of your session to capture the additional dependencies.
5. Commit any changes to `renv.lock`.

Note that when you open up the `training-modules.Rproj`, the `.Rprofile` file makes it such that the `renv` library is loaded and the repositories in the `renv.lock` file will be set with `options(repos = ...)`.

### How we use `renv` with Docker

We use the `renv.lock` file in this repository to install R dependencies on the image per the [`renv` documentation for creating Docker images](https://rstudio.github.io/renv/articles/docker.html#creating-docker-images-with-renv-1).
Specifically, we copy the `renv.lock` file and run `renv::restore()`, which installs the packages in the lockfile into the Docker image's system library.

In practice, this means that you will not need to add individual R packages to the `Dockerfile`, but you may have to add _system dependencies_ (e.g., via `apt-get`) required for the installation of those R packages.

## Docker Image

### Developing within the Docker container

To use the Docker image for development, pull from Dockerhub with:

```
docker pull ccdl/training_dev:latest
```

To run the container and mount a local volume, use the following from the root of this repository:

```
docker run \
  --mount type=bind,target=/home/rstudio/training-modules,source=$PWD \
  -e PASSWORD=<PASSWORD> \
  -p 8787:8787 \
  ccdl/training_dev:latest
```

Replacing `<PASSWORD>` with the password of your choice.
You can then navigate to `localhost:8787` in your browser and login with username `rstudio` and the password you just set via `docker run`.

To work on the project, you should then open `training-modules/training-modules.Rproj` on the RStudio server, then proceed as in the [typical development workflow](#typical-development-workflow).

### Testing Docker image builds via GitHub Actions

When a pull request changes either the `Dockerfile` or `renv.lock`, a GitHub Actions workflow (`build-docker.yml`) will be run to test that the image will successfully build.

### Pushing to Dockerhub via GitHub Actions

When a pull request is merged into `master`, the `build-and-push-docker.yml` GitHub Actions workflow will be triggered. 
The project Docker image will be rebuilt and pushed as `ccdl/training_dev:latest`.

## Spell check

We perform spell checking for every pull request to `master` as part of a GitHub Actions workflow (`spell-check.yml`); it is designed to fail if more than 3 spelling errors are detected.
You can see what errors are detected in `stdout` for the `Run spell check` step of the workflow.
This workflow uses a script, `scripts/spell-check.R`, to spell check `.md` and completed `.Rmd` files.
The custom dictionary we use for the project can be found at `components/dictionary.txt`. 

To run spell check locally, you can run the following from the root of the repository:

```
Rscript --vanilla scripts/spell-check.R
```

The spelling errors will be listed in `spell_check_errors.tsv` in the root of the repo; this file is ignored by git.
