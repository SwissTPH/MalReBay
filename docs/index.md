# MalReBay

## Overview

**MalReBay** is an R package for Bayesian molecular correction of
malaria therapeutic efficacy studies (TES). Given paired genotyping data
from a patient’s initial infection (Day 0) and a later recurrence, it
estimates the posterior probability that the recurrence is a
**recrudescence** (treatment failure, same parasite persisting) versus a
**reinfection** (new parasite acquired after treatment).

Unlike rule-based allele-matching methods, MalReBay:

- models sequencing error, allelic dropout, and allele loss explicitly
- handles **polyclonal infections** (MOI \> 1) at every locus
- incorporates **local allele frequency** information to weigh evidence
- uses a full **Bayesian MCMC** engine (Stan HMC-NUTS) for principled
  uncertainty quantification
- supports both **length-polymorphic markers** (microsatellites, MSP,
  GLURP) and **amplicon sequencing** (haplotype) data

------------------------------------------------------------------------

## Key Features

| Feature | Detail |
|----|----|
| **Data types** | Length-polymorphic markers and amplicon sequencing |
| **MCMC engine** | Stan HMC-NUTS via `rstan` — fast, well-mixing, gradient-based |
| **Polyclonal infections** | Analytically marginalises over within-host clone multiplicity |
| **Error model** | Estimates `q_mismatch`, `q_loss`, `q_dropout` from the data |
| **Prior** | Recrudescence probability fixed at 0.5 (no directional bias) |
| **Convergence** | Rank-normalised R̂, bulk/tail ESS, Geweke, Gelman-Rubin |
| **Output** | Per-patient posterior probability + match-counting comparison table |

------------------------------------------------------------------------

## Installation

``` r
# Install the development version from GitHub
remotes::install_github("SwissTPH/MalReBay")
```

> **Note:** MalReBay links against `rstan` and `StanHeaders`. On first
> use the Stan model is compiled once and cached. This may take a few
> minutes but only happens once per R installation.

------------------------------------------------------------------------

## Quick Start

``` r
library(MalReBay)

# Run the full pipeline on the bundled Angola 2021 TES example data
results <- MalReBay(output_folder = "my_results")

# View posterior probabilities
head(results$posterior_probabilities)
```

------------------------------------------------------------------------

## Workflow

MalReBay wraps a four-step pipeline that can also be called step by
step:

``` R
import_data()
    ↓
classify_infections()
    ↓
summarise_results()
    ↓
save_results()
```
