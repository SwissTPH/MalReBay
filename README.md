
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

## Installation

Installing MalReBay takes three short steps and about 10 minutes in
total, almost all of it hands-off waiting. You only need to do this
**once** per computer — after that, `library(MalReBay)` is all you
need, in every future R session.

> **Why so many steps?** MalReBay’s statistical engine is written in
> [Stan](https://mc-stan.org/), a language for Bayesian modelling that
> needs to be compiled into fast machine code. Step 1 installs that
> compiler toolkit (CmdStan) and step 2 checks it’s ready to go. Step 3
> installs MalReBay itself, which compiles its model against the
> CmdStan from step 1 — that’s why **the order below matters**: CmdStan
> has to be in place *before* you install MalReBay.

### Prerequisites

- **R version 4.1.0 or later.** Check yours by running
  `R.version.string` in the R console.
- **A C++ compiler**, needed to build CmdStan and MalReBay’s model.
  This is usually a one-time, one-command install:
  - **Windows:** install
    [Rtools](https://cran.r-project.org/bin/windows/Rtools/) (pick the
    version matching your R version).
  - **macOS:** open Terminal and run `xcode-select --install`.
  - **Linux (Debian/Ubuntu):** open a terminal and run
    `sudo apt-get install build-essential`.

### Step 1 — Install CmdStan

``` r
# install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()
```

This downloads and compiles CmdStan, MalReBay’s statistical engine — it
takes around 5 minutes and shows its own progress as it goes.

### Step 2 — Check your setup (optional but recommended)

``` r
cmdstanr::check_cmdstan_toolchain()
```

This confirms your C++ compiler is correctly set up *before* you try
to install MalReBay, so any problem is caught here with a clear
message rather than a confusing error later. If it reports a problem,
follow its instructions (usually pointing back to the Prerequisites
step above), then run it again.

### Step 3 — Install MalReBay

``` r
# install.packages("remotes")  # run this line first if you don't have remotes
remotes::install_github("SwissTPH/MalReBay")
```

This step compiles MalReBay’s Stan model against the CmdStan from step
1, so it takes a little longer than a typical package install (roughly
a minute) — that one-time cost is also why the order matters: if
CmdStan isn’t installed yet, this step will fail.

**That’s it!** Open a fresh R session and confirm everything worked:

``` r
library(MalReBay)
```

If this loads without a `NOTE: CmdStan is not installed...` message,
you’re ready to go — see [Quick Start](#quick-start) below.

### Troubleshooting

- **“unused argument” or a compiler error during step 3:** re-run
  `cmdstanr::check_cmdstan_toolchain()` from step 2 — it usually
  points directly at the missing piece (most often a missing C++
  compiler).
- **Installed once already, now reinstalling after making changes to
  the package:** every reinstall recompiles the Stan model from
  scratch, so step 3 will always take about a minute, even for small
  changes elsewhere in the package.
- **Still stuck?** Open an issue at
  <https://github.com/SwissTPH/MalReBay/issues> with the exact error
  message — it’s the fastest way for us to help.

## Quick Start

Try MalReBay right away on its bundled example dataset — no files of
your own needed yet:

``` r
library(MalReBay)

results <- MalReBay()

head(results$posterior_probabilities)
```

Once you’re ready to use your own data, point `MalReBay()` at your
files:

``` r
results <- MalReBay(
  filepath        = "path/to/genotype_data.xlsx",
  marker_filepath = "path/to/marker_info.xlsx",
  output_folder   = "my_results"
)

head(results$posterior_probabilities)
```

For a full walkthrough of the input file formats and a worked example,
see the [package vignette](https://swisstph.github.io/MalReBay/).
