# Reproducible workflow for Tansley review on the minimum conductance

This repository contains code needed to reproduce the article:

**Remko A. Duursma1, Chris J. Blackman1, Rosana Lopéz1, Nicolas K. Martin-StPaul2, Hervé Cochard3, Belinda E. Medlyn1** On the minimum leaf conductance: its role in models of plant water use, and ecological and environmental controls. New Phytologist, in review.


## Instructions

Simply source `run.R` to install dependencies, read data, compile manuscript and figures. The workflow uses `rmarkdown` for the manuscript, and produces `manuscript.docx` (and a supporting information document), and figures in `output/`. 

One dependency (the `speciesmap` package) is not on CRAN, and has to be installed with: 

```
devtools::install_github("remkoduursma/speciesmap")
```

The workflow also depends on an installation of pandoc (or use Rstudio if you don't).

The first time you run this, it will be quite slow because species observation data and climate data will be downloaded and processed.
