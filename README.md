[![R build status](https://github.com/hubentu/MACSr/workflows/R-CMD-check/badge.svg)](https://github.com/hubentu/MACSr/actions)

# MACSr
MACS3 R/BioC wrapper

## Installation
``` r
devtools::install_github("macs3-project/MACSr")
```

From Bioconductor:
``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("MACSr")
```

The package is built on
[basilisk](https://bioconductor.org/packages/release/bioc/html/basilisk.html). The
dependent python library
[macs3](https://github.com/macs3-project/MACS) will be installed
automatically inside its conda environment.

## User Guide
``` r
vignette(package = "MACSr")
```
