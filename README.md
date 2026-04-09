
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MalReBay

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/SwissTPH/MalReBay/graph/badge.svg)](https://app.codecov.io/gh/SwissTPH/MalReBay)
[![R-CMD-check](https://github.com/SwissTPH/MalReBay/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/SwissTPH/MalReBay/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

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

    import_data()
        ↓
    classify_infections()
        ↓
    summarise_results()
        ↓
    save_results()

### Step 1 — Import data

``` r
imported_data <- import_data(
  filepath        = "path/to/genotype_data.xlsx",
  marker_filepath = "path/to/marker_details.xlsx"
)
```

`import_data()` reads and validates your Excel file, detects the data
type (`length_polymorphic` or `ampseq`) automatically from the marker
metadata, and returns a structured list used by all downstream
functions.

### Step 2 — Run MCMC classification

``` r
mcmc_config <- list(
  n_chains     = 4,       # independent chains (run in parallel)
  iter         = 2000,    # total iterations per chain (including warmup)
  burn_in_frac = 0.5,     # fraction used as warmup
  adapt_delta  = 0.90,    # Stan step-size adaptation target (0.8–0.99)
  random_seed  = 42
)

mcmc_results <- classify_infections(
  imported_data = imported_data,
  mcmc_config   = mcmc_config
)
```

### Step 3 — Summarise and save results

``` r
summary_results <- summarise_results(
  mcmc_results  = mcmc_results,
  imported_data = imported_data,
  output_folder = "my_results"
)
```

------------------------------------------------------------------------

## Input Data Format

### Length-polymorphic markers

The Excel file must contain:

| Column | Description |
|----|----|
| `Sample.ID` | Patient identifier with suffix `" Day 0"` or `" recurrence"` |
| `Site` | Study site name |
| `MarkerName_1`, `_2`, … | Fragment lengths in base pairs (one column per allele slot) |

Example marker names: `313_1`, `313_2`, `TA1_1`, `TA1_2`, `MSP1_1` …

### Amplicon sequencing data

Same structure but marker columns contain **haplotype strings** instead
of fragment lengths:

| Column | Description |
|----|----|
| `Sample.ID` | As above |
| `Site` | As above |
| `MarkerName_allele_1`, `_allele_2`, … | Haplotype labels (e.g., `SERA2-1`, `SERA2-2`) |

### Marker metadata file

A separate Excel file (`makers_details.xlsx`) specifies one row per
marker:

| Column           | Description                                            |
|------------------|--------------------------------------------------------|
| `marker_id`      | Matches the prefix of marker columns in the data file  |
| `binning_method` | `microsatellite`, `msp_glurp`, or `exact` (for ampseq) |
| `repeatlength`   | Repeat unit length (for microsatellites)               |

------------------------------------------------------------------------

## MCMC Configuration Parameters

| Parameter | Default | Description |
|----|----|----|
| `n_chains` | 4 | Number of independent chains |
| `iter` | 2000 | Total iterations per chain (warmup + sampling) |
| `burn_in_frac` | 0.5 | Fraction of iterations discarded as warmup |
| `adapt_delta` | 0.85 | Stan HMC step-size target; increase to 0.95 if divergent transitions occur |
| `random_seed` | 42 | Random seed for reproducibility |

A default configuration file is bundled with the package and used
automatically when no `mcmc_config` argument is supplied.

------------------------------------------------------------------------

## Outputs

`MalReBay()` returns a named list and writes the following files to
`output_folder`:

| File | Content |
|----|----|
| `posterior_probabilities.csv` | Per-patient posterior probability of recrudescence, number of comparable loci |
| `bayesian_match_counting_comparison.csv` | Side-by-side comparison of Bayesian probabilities and simple allele-match counts |
| `mcmc_convergence_summary.csv` | Per-site R̂, bulk ESS, tail ESS diagnostics |
| `recrudescence_probability_histogram.png` | Distribution of posterior probabilities across all patients |
| `convergence_diagnosis/` | Gelman-Rubin, trace, and ACF plots per site |

### Interpreting posterior probabilities

| Posterior probability | Interpretation                     |
|-----------------------|------------------------------------|
| \> 0.9                | Strong evidence of recrudescence   |
| 0.5 – 0.9             | Moderate evidence of recrudescence |
| 0.1 – 0.5             | Moderate evidence of reinfection   |
| \< 0.1                | Strong evidence of reinfection     |

A threshold of 0.5 is the standard decision boundary, consistent with
the equal (0.5 / 0.5) prior used by the model.

------------------------------------------------------------------------

## Example Datasets

Two example datasets are bundled with the package:

| Dataset | File | Type | Patients | Sites | Markers |
|----|----|----|----|----|----|
| Angola 2021 TES | `Angola_2021_TES_7NMS.xlsx` | Length-polymorphic | 70 | 3 | 7 microsatellites |
| Amplicon sequencing | `Amplicon_Sequencing.xlsx` | Ampseq | — | — | SERA2, TRAP3, cpp2, csp2 |

``` r
# Length-polymorphic example
lp_file     <- system.file("extdata", "Angola_2021_TES_7NMS.xlsx", package = "MalReBay")
marker_file <- system.file("extdata", "makers_details.xlsx",       package = "MalReBay")
results_lp  <- MalReBay(filepath = lp_file, marker_filepath = marker_file,
                         output_folder = "results_lp")

# Amplicon sequencing example
ampseq_file  <- system.file("extdata", "Amplicon_Sequencing.xlsx", package = "MalReBay")
results_amps <- MalReBay(filepath = ampseq_file, marker_filepath = marker_file,
                          output_folder = "results_ampseq")
```

------------------------------------------------------------------------

## Citation

If you use MalReBay in your research, please cite:

> Ochieng V.A., Plucinski M., Golumbeanu M. *MalReBay: Bayesian
> classification for malaria recurrences.* Swiss Tropical and Public
> Health Institute. <https://github.com/SwissTPH/MalReBay>

The statistical methodology builds on:

> Plucinski M.M. et al. (2015). Robust algorithm for systematic
> classification of malaria late treatment failures as recrudescence or
> reinfection using microsatellite genotyping. *Antimicrobial Agents and
> Chemotherapy*, 59(10).

------------------------------------------------------------------------

## Authors

- **Veronica Adhiambo Ochieng** — AIMS Rwanda / Swiss TPH
  (`veronica.adhiambo@aims.ac.rw`)
- **Mateusz Plucinski** — CDC Malaria Branch
- **Monica Golumbeanu** — Swiss Tropical and Public Health Institute
  (`monica.golumbeanu@swisstph.ch`)

------------------------------------------------------------------------

## Links

- Package website: <https://swisstph.github.io/MalReBay>
- Source code: <https://github.com/SwissTPH/MalReBay>
- Bug reports: <https://github.com/SwissTPH/MalReBay/issues>
