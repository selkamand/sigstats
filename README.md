
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sigstats

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/sigstats)](https://CRAN.R-project.org/package=sigstats)
<!-- badges: end -->

**sigstats** enables common mathematical operations / transformations to
be applied to **sigverse** style signatures / decompositions

## Installation

You can install the development version of sigstats like so:

``` r
# install.packages('remotes')
remotes::install_github('selkamand/sigstats')
```

## Quick Start

``` r
library(sigstats)
library(sigstash)
library(sigvis)

# Load a signature collection
signatures <- sig_load("COSMIC_v3.3.1_SBS_GRCh38")

# Create a model that represents a mix of SBS1 (40%) and SBS2 (60%)
model <- c(SBS1 = 0.4, SBS2 = 0.6)

# Add selected signatures to the combined model
combined_signatures <- sig_combine(signatures, model)

# Visualise result
sig_visualise(combined_signatures, type = "model")
```
