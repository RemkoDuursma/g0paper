
# Read & download data, load & install packages.
message("Loading packages, reading data...")
source("R/load.R")

# Make figures as PDF (see output/)
message("Making figures...")
source("R/make_figures.R")

# Convert manuscript (output is manuscript.docx & manuscript_suppinfo.docx).
message("Rendering documents...")
rmarkdown::render("manuscript.Rmd", encoding = "utf8", quiet=TRUE)
rmarkdown::render("manuscript_suppinfo.Rmd", encoding = "utf8", quiet=TRUE)
