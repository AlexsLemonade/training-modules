# The ready script

# C. Savonen
# CCDL for ALSF 2019

# This script is meant to be run in a training workshop to test if everything is
# set up.

# Usage: source(file.path("kitematic", "ready_script.R"))

# Read in the image
image <- png::readPNG(file.path("kitematic", "ready.png"))

# Raster the image
grid::grid.raster(image)

# Remove image object from environment
rm(image)
