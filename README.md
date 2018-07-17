[![DOI](https://zenodo.org/badge/84373768.svg)](https://zenodo.org/badge/latestdoi/84373768)

# Reproducible workflow for Tansley review on the minimum conductance

This repository contains code needed to reproduce the article:

**On the minimum leaf conductance: its role in models of plant water use, and ecological and environmental controls**

by 
Remko A. Duursma, Chris J. Blackman, Rosana Lopéz, Nicolas K. Martin-StPaul2, Hervé Cochard, Belinda E. Medlyn. In review at *New Phytologist* as a 'Tansley review'.


## Instructions

This workflow depends on an installation of [R](https://www.r-project.org/). Missing packages will be installed on the fly. The workflow also depends on an installation of [pandoc](https://pandoc.org/), but if you use Rstudio, you do not need to install it again.

After cloning this repository, you can run all analyses, make figures (in `output/`), and compile the manuscript into docx, with this command:

```
source("run.R")
```

(or `Rscript run.R` from the command line).


If you want to re-generate the species climate envelopes (see `R/load.R` for details), you will need the `speciesmap` package, which can be installed with: 

```
devtools::install_github("remkoduursma/speciesmap")
```

On Windows, you need to install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) first.

## Data

All raw data files are included in this repository, except those that can be downloaded from a permanent URL (at figshare). Please read the source file `R/load.R` carefully to find the raw data files, how they are read, processed, etc.

The literature compilation of minimum conductance data is [available from this link](https://www.github.com/remkoduursma/gmindatabase).

## Zenodo

This repository is [archived at Zenodo](https://zenodo.org/account/settings/github/repository/RemkoDuursma/g0paper)
