
# Read & download data, load & install packages.
source("R/load.R")

# Make figures as PDF (see output/)
# Conversion to png relies on an installation of ImageMagick (and GhostScript); may fail on your system.
source("R/make_figures.R")

# Convert manuscript (output is manuscript.docx).
rmarkdown::render("manuscript.Rmd", encoding = "utf8")
