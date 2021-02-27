#! /bin/bash
set -euo pipefail

# This script is used to fetch exercise notebooks from the exercise-notebook-answers repo
# Because the exercise-notebook-answers repo is private, the script requires that a 
# GH_TOKEN environment variable be set with token for an account with access to that repo.

# the exercise-notebook-answers url where files will be downloaded from
base_url=https://${GH_TOKEN}@raw.githubusercontent.com/AlexsLemonade/exercise-notebook-answers/master

# url for exercise-list.txt
list_url=${base_url}/components/exercise-list.txt

exercise_paths=$(curl -s ${list_url})

for path in $exercise_paths
do
  # get exercise files and put them in the right place
  curl -s ${base_url}/${path} > ${path}
done
