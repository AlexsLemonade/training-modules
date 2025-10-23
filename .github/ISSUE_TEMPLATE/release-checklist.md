---
name: Release checklist
about: Prepare for a new workshop release of training-modules
title: "Prepare for workshop release: YYYY-month"
labels: release
assignees: ''

---

## Steps for a workshop release of `training-modules`
Use this checklist to prepare for a new release of `training-modules` for a training event.
Update the title of this issue to reflect the version you are releasing (e.g. "Prepare for workshop release: 2025-august").

Detailed instructions can be found in the [workshop releases documentation](https://github.com/AlexsLemonade/training-modules/blob/master/workshop-releases.md).

### Preparing for the release and testing

- [ ] Are all of the updates planned for this training resolved? If there are any issues that are unresolved, mark this issue as blocked by those.
- [ ] Run the [Make Live Notebooks](https://github.com/AlexsLemonade/training-modules/actions/workflows/make-live.yml) GitHub Action to create up-to-date versions of the  `-live.Rmd` notebooks and the rendered HTML files.
- [ ] In the [`exercise-notebook-answers`](https://github.com/AlexsLemonade/exercise-notebook-answers) repository, run the [Unsolve exercise notebooks](https://github.com/AlexsLemonade/exercise-notebook-answers/actions/workflows/unsolve.yml) GitHub Action to create up-to-date versions of the `-exercise.Rmd` notebooks without answers. This should automatically trigger the [Copy exercises to training](https://github.com/AlexsLemonade/exercise-notebook-answers/actions/workflows/exercises-to-training.yml) GitHub Action workflows to ensure that the exercise notebooks are up-to-date and that the answers are copied to the `training-modules` repository
- [ ] Update [current-modules.json](https://github.com/AlexsLemonade/training-modules/blob/master/current-modules.json) with the release tag and the modules that will be included in this training event.
- [ ] Test that the Docker image has been built correctly and that the modules are available as expected in the `/etc/skel` directory.

### Creating a release

- [ ] On the [releases page](https://github.com/AlexsLemonade/scpcatools/releases), choose `Draft a new release`.
- [ ] Create a release using the same tag name that you specified in the `current-modules.json` file.
This will again trigger the [Build Docker Image](https://github.com/AlexsLemonade/training-modules/actions/workflows/build-docker.yml) action, tagging a Docker release to be associated with the training workshop.
- [ ] Publish the release!
