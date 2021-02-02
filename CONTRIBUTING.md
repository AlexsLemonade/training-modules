# Contributing and Development Guidelines

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

```
TODO: Pulling the image once #303 goes in
```

### Testing Docker image builds via GitHub Actions

When a pull request changes either the `Dockerfile` or `renv.lock`, a GitHub Action workflow (`build-docker.yml`) will be run to test that the image will successfully build.
