
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MalReBay

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/SwissTPH/MalReBay/graph/badge.svg)](https://app.codecov.io/gh/SwissTPH/MalReBay)
[![R-CMD-check](https://github.com/SwissTPH/MalReBay/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/SwissTPH/MalReBay/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**MalReBay** is an R package for Bayesian molecular correction of
malaria therapeutic efficacy studies (TES). Given paired genotyping data
from a patient’s Day 0 infection and a later recurrence, it estimates
the posterior probability of **recrudescence** (treatment failure)
versus **reinfection** (new parasite). It supports length-polymorphic
markers (microsatellites, MSP, GLURP) and amplicon sequencing data, and
uses Stan HMC-NUTS via `cmdstanr` for inference.

For a full walkthrough see the [package
vignette](https://swisstph.github.io/MalReBay/).

## Installation

MalReBay is installed in two steps. The second step is a one-time setup
and only needs to be done once, even across future R sessions.

**Install MalReBay** from GitHub:

``` r
# install.packages("remotes")  # run this line first if you don't have remotes
remotes::install_github("SwissTPH/MalReBay")
```

**Install CmdStan** (the statistical engine MalReBay uses for Bayesian
inference):

``` r
cmdstanr::install_cmdstan()
```

> CmdStan is a standalone statistical computing toolkit that MalReBay
> uses behind the scenes to run its Bayesian model. Think of it like a
> calculator engine MalReBay provides the interface, CmdStan does the
> heavy computation. This download takes around 5 minutes and only needs
> to be done once on your computer.

## Quick Start

``` r
library(MalReBay)

results <- MalReBay(
  filepath        = "path/to/genotype_data.xlsx",
  marker_filepath = "path/to/marker_info.xlsx",
  output_folder   = "my_results"
)

head(results$posterior_probabilities)
```
