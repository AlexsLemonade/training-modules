# Preparing for a Data Lab workshop release

These instructions cover the steps needed to prepare for a Data Lab workshop, including preparing the user `skel` directory, building the Docker image, testing, and preparing the release tag.

We tag a release of the training modules repository for each training workshop that the Data Lab conducts.
These tagged versions are used both as a reference and to create and push a Docker image specific to that training session.

Before tagging a release, please ensure that the following steps have been completed.

## Prepare teaching and exercise notebooks

- Run the [Make Live Notebooks](https://github.com/AlexsLemonade/training-modules/actions/workflows/make-live.yml) GitHub Action to create up-to-date versions of the  `-live.Rmd` notebooks and the rendered HTML files.
- In the [`exercise-notebook-answers`](https://github.com/AlexsLemonade/exercise-notebook-answers) repository, run the [Unsolve exercise notebooks](https://github.com/AlexsLemonade/exercise-notebook-answers/actions/workflows/unsolve.yml) and [Copy exercises to training](https://github.com/AlexsLemonade/exercise-notebook-answers/actions/workflows/exercises-to-training.yml) GitHub Actions workflows to ensure that the exercise notebooks are up-to-date and that the answers are copied to the `training-modules` repository.

## Select training modules in `current-modules.json`

First, update the `current-modules.json` file in the `training-modules` repository to reflect the modules that will be used in the training workshop:
  - Set the `release-tag` field to the tag you plan to create for this training workshop.
    We have usually used the format `2025-june`.
  - Set the `modules` to a list of the modules that will be used in the training workshop.
    These modules will be copied to the Docker image `/etc/skel` directory with only the `-live.Rmd` notebooks, ready for participants to use.
  - Set the `reference-modules` to a list of the modules that will be used as reference material for the training workshop.
    For example, when we teach the `scRNA-seq-advanced` module, we included the `scRNA-seq` module as a reference of the material we are building from.
    These modules will be copied to the Docker image `/etc/skel` directory with completed notebooks.

An example `current-modules.json` file (set for an advanced single cell workshop) might look like this:

```json
{
  "release-tag": "2025-dev",
  "modules": ["scRNA-seq-advanced"],
  "reference-modules": ["scRNA-seq"]
}
```

Once the `current-modules.json` file is set, file a pull request with these changes to `master`.
  - Merging the pull request will trigger the [Build Docker Image](https://github.com/AlexsLemonade/training-modules/actions/workflows/build-docker.yml) Action, which will build and push a Docker image to `ccdl/training_rstudio:edge` with the specified modules in the `/etc/skel` directory (created by the `setup-skel.py` script).

## Test the Docker image

Once the [Build Docker Image](https://github.com/AlexsLemonade/training-modules/actions/workflows/build-docker.yml) Action is complete, you can test that the Docker image has been built correctly and that the modules are available as expected in the `/etc/skel` directory that will be used as a template when creating user accounts.

One important thing to note is that there are a number of symbolic links in the `/etc/skel` directory that point to shared data files in the `/shared` directory on the Data Lab RStudio server.
If running the Docker image locally, these links will be broken, but they should work correctly on the RStudio server if the `/shared` directory is mounted.

### Testing on the RStudio Server

To test the Docker image on the Data Lab RStudio server, follow the instructions in the `rstudio-server` README: <https://github.com/AlexsLemonade/rstudio-server/tree/dev?tab=readme-ov-file#testing-edge>

Briefly, you will need to create a new user of the `edge` session type, which uses the `ccdl/training_rstudio:edge` Docker image and deploy that change to production.
Then, log in as that user via RStudio server and check that the `home` directory contains the expected contents from the `/etc/skel` directory, which will look something like (assuming the user is `dlab`):

```bash
$ ls -lah ~
drwxr-x--- 6 dlab dlab 4.0K Nov 19 20:16 .
drwxr-xr-x 1 root root 4.0K Nov 19 20:13 ..
-rw-r--r-- 1 dlab dlab  220 Jan  6  2022 .bash_logout
-rw-r--r-- 1 dlab dlab 3.7K Jan  6  2022 .bashrc
drwxr-xr-x 3 dlab dlab 4.0K Nov 19 20:16 .config
drwxr-xr-x 3 dlab dlab 4.0K Nov 19 20:16 .local
-rw-r--r-- 1 dlab dlab  807 Jan  6  2022 .profile
drwxr-xr-x 3 dlab dlab 4.0K Nov 19 20:16 R
lrwxrwxrwx 1 dlab dlab   12 Nov 18 22:28 shared-data -> /shared/data
drwxr-xr-x 4 dlab dlab 4.0K Nov 18 22:28 training-modules
```

Within the `training-modules` directory, you should see the modules that you specified in the `current-modules.json` file, and you can check that the files are present as expected.
You should also be able to check that the symbolic links to the `/shared/data/` directory are working correctly.
For example, links within modules' `data/` directories should point to the correct files (usually in `/shared/data/training-modules/`).

After verifying that the Docker image is working as expected, you can cleanup the test user as described in the `rstudio-server` README.

### Testing locally

You can also test the Docker image locally, which may be more convenient for testing via the RStudio GUI, but note that the symbolic links to the `/shared/data` directory will not work correctly (unless you have a local copy of the data you can mount to the same path).

To test the Docker image locally, you can use the following commands to pull the Docker image and run launch the container with an RStudio server session:

```bash
docker pull --platform linux/amd64 ccdl/training_rstudio:edge
docker run \
  -e PASSWORD={PASSWORD} \
  -p 8787:8787 \
  ccdl/training_rstudio:edge
```

Replace `{PASSWORD}` with a password of your choice.
Then you can open a web browser and go to `http://localhost:8787` to access the RStudio server, logging in with the username `rstudio` and the password you specified.

As noted, the symbolic links to the `/shared/data` directory will not work correctly in this context.
For example, if you run `ls /etc/skel` in the terminal, you will see a visual indicator (e.g., red text on a black background) that the `shared-data` link is broken.
Non-link files within the `/etc/skel/training-modules/` directory, however, should be present as expected.

## Tagging the release

Once you have tested the Docker image and are satisfied that it contains all of the expected modules and files for the coming workshop, [create a release in the `training-modules` repository](https://github.com/AlexsLemonade/training-modules/releases) using the same tag name that you specified in the `current-modules.json` file.
This will again trigger the [Build Docker Image](https://github.com/AlexsLemonade/training-modules/actions/workflows/build-docker.yml) action, tagging a Docker release to be associated with the training workshop.
