
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rq

The goal of rq is to calculate the Relative Quantity of qPCR sets.

## Installation

You can install the latest version of rq by downloading the latest
release from this repository and executing:

``` r
# install.packages('rq_<version>.tar.gz', repos = NULL, type = 'source', dependencies = TRUE)
```

or directly from R with an internet connection using `devtools`:

``` r
# devtools::install_github('godrigos/rq', dependencies = TRUE)
```

## Example

This is an example on how to use rq:

``` r
library(rq)
RQ <- rq(qpcr, "pUC18", "WT", "ITS2")
RQ
```

For more details consult the package help files inside R.
