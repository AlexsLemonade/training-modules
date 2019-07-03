# The ready script

# C. Savonen
# CCDL for ALSF 2019

# This script is meant to be run in a training workshop to test if everything is
# set up.

# Usage: source("kitematic/ready_script.R")

# Read in the image
image <- png::readPNG("kitematic/ready.png")

# Raster the image
grid::grid.raster(image)

# Remove image object from environment
rm(image)
