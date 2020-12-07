[![R build status](https://github.com/hubentu/MACSr/workflows/R-CMD-check/badge.svg)](https://github.com/hubentu/MACSr/actions)

# MACSr
MACS3 R/BioC wrapper

## Installation
The package is dependent on the python library `macs3`. Please install it first.
```sh
pip3 install macs3
```

Then it can be installed.
```{r}
devtools::install_github("macs3-project/MACSr")
```

## User Guide
``` r
vignette(package = "MACSr")
```
