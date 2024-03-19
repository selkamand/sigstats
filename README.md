
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sigstats

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/sigstats)](https://CRAN.R-project.org/package=sigstats)
<!-- badges: end -->

**sigstats** enables common mathematical operations / transformations to
be applied to **sigverse** style signatures / catalogues

## Installation

You can install the development version of sigstats like so:

``` r
if (!require("pak", quietly = TRUE))
    install.packages("pak")

pak::pak("selkamand/sigstats")
```

## Quick Start

``` r
library(sigstats)
library(sigstash)
library(sigvis)

# Load a signature collection
signatures <- sig_load("COSMIC_v3.3.1_SBS_GRCh38")

# Create a model that represents a mix of SBS2 (40%) and SBS13 (60%)
model <- c(SBS2 = 0.4, SBS13 = 0.6)

# Add selected signatures to the combined model
combined_signatures <- sig_combine(signatures, model)

# Visualise result
sig_visualise(
  combined_signatures, 
  class = "model",
  title = "Model",
  subtitle = "Created by combining SBS2 (40%) and SBS13 (60%)"
  )
#> ✔ All channels matched perfectly to set [sbs_96]. Using this set for sort order
#> ✔ All types matched perfectly to set [sbs_type]. Using this set for sort order
#> ✔ Types matched perfectly to palette [snv_type]
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />
