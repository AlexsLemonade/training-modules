library(renv)
source("renv/activate.R")

# Set the repos using the renv.lock file
renv_json <- jsonlite::read_json("renv.lock")
renv_r_repos <- renv_json$R$Repositories

# Extract the names
repo_names <- purrr::flatten_chr(
  purrr::map(renv_r_repos,
             ~ .x$Name)
)

# Extract the URLs
repo_urls <- purrr::flatten_chr(
  purrr::map(renv_r_repos,
             ~ .x$URL)
  )

# Set the repo names
names(repo_urls) <- repo_names

# Set the options
options(repos = repo_urls)
  
# Remove all these objects
rm(renv_json, renv_r_repos, repo_names, repo_urls)
