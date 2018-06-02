[![DOI](https://zenodo.org/badge/84373768.svg)](https://zenodo.org/badge/latestdoi/84373768)

# Reproducible workflow for Tansley review on the minimum conductance

This repository contains code needed to reproduce the article:

**On the minimum leaf conductance: its role in models of plant water use, and ecological and environmental controls**

by 
Remko A. Duursma, Chris J. Blackman, Rosana Lopéz, Nicolas K. Martin-StPaul2, Hervé Cochard, Belinda E. Medlyn. In review at *New Phytologist* as a 'Tansley review'.


## Instructions

From [R](https://www.r-project.org/), you can run all analyses, make figures, and compile the manuscript into docx, with this command:

```
source("run.R")
```

The workflow uses `rmarkdown` for the manuscript, and produces `manuscript.docx` (and a supporting information document), and figures as PDF in `output/`. 

One dependency (the `speciesmap` package) is not on CRAN, and has to be installed with: 

```
devtools::install_github("remkoduursma/speciesmap")
```

On Windows, you need to install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) first.

The workflow also depends on an installation of [pandoc](https://pandoc.org/), but if you use Rstudio, it comes with an installation of pandoc already.

The first time you run the complete workflow, it will be quite slow because species observation data and climate data will be downloaded (you do need an internet connection the first run) and processed.


## Zenodo

This repository is archived at Zenodo, with the doi: [10.5281/zenodo.1258003](doi.org/10.5281/zenodo.1258003)
