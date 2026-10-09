# Contributing to MalReBay

Thanks for considering a contribution. This document covers how the package
is organised internally and the mechanics of making a change — for the
statistical model itself and how to *use* the package, see the
[tutorial vignette](vignettes/MalReBay.qmd) instead.

## Getting set up

Follow the [README's installation instructions](README.md#installation)
first (CmdStan, then MalReBay) — you need a working CmdStan install to run
almost anything in this package, including the test suite.

Then clone the repo and install the package from source:

```r
# install.packages("devtools")
devtools::install_deps(dependencies = TRUE)
devtools::load_all()
```

## Package structure

MalReBay classifies each recurrence in four stages — `import_data()` →
`classify_infections()` → `summarise_results()` → `save_results()`, wrapped
together by `MalReBay()` (see the vignette's workflow diagram). Internally,
the code is organised by what part of that pipeline it belongs to:

### Data preparation (`R/import_data.R`, `R/allele_utils.R`, `R/quality_control.R`)

- `import_data()` reads the genotyping Excel/CSV and marker metadata,
  auto-detects length-polymorphic vs. AmpSeq data from the first allele
  value, and validates the input structure.
- `define_alleles()` / `recodeallele()` group raw length-polymorphic
  fragment sizes into allele bins, using each marker's `repeatlength` and
  `binning_method`.
- `detect_msp_variants()` identifies which loci belong to the MSP1
  (K1/MAD20/RO33) and MSP2 (3D7/FC27/IC) families, needed for the classic
  WHO 2/3-3/3 trio rule.
- `compute_locus_comparability()` and `data_quality_check()` compute
  per-patient locus availability and run basic sanity checks on the input.

### Stan data prep and Bayesian inference (`R/stan_data_prep.R`, `R/stan_interface.R`, `src/stan/`)

- `prepare_stan_data()` (+ `validate_stan_data()`, `stan_data_only()`) turns
  the imported data frame into the integer-coded arrays the Stan model
  expects: recoded alleles, hidden-allele flags for polyclonal infections,
  the locus-comparability matrix, and per-locus allele-distance matrices.
- `run_stan_sites()` / `extract_stan_results()` load the precompiled model
  via `instantiate::stan_package_model()` and run `$sample()` once per TES
  site.
- `src/stan/malrebay_model.stan` is the actual model — a single Stan file
  shared by all three marker types; per-locus dispatch on `binning_method`
  (microsatellite distance decay / MSP-GLURP family clustering / AmpSeq
  exact match) happens inside the model via the `method_int` array, not by
  branching in R. It's precompiled at package install time (see
  `src/install.libs.R`), which ships a portable `stanc` binary and a stub
  CmdStan directory so installation doesn't require a full CmdStan
  toolchain — only *using* the package does (see README).

### Traditional match counting (`R/match_counting.R`, WHO helpers in `R/plot_utils.R`)

- `perform_match_counting()` (+ `assign_clusters()`) implements the
  deterministic allele match-counting algorithm — per-locus `R`/`NI`/`IND`/`ERR`
  calls — used as the classical comparison against MalReBay's Bayesian
  probability.
- `build_who_table()`, `apply_who_rule()`, and `combine_msp_variants()`
  turn those per-locus calls into the WHO loose/strict pass-fail
  classification: the 2/3-3/3 rule for MSP1/MSP2 trio panels, or a
  proportional 70%/100% rule otherwise.

### Results and output (`R/MalReBay.R`, `R/plot_utils.R`, `R/quality_control.R`)

- `classify_infections()`, `summarise_results()`, `save_results()`, and
  `MalReBay()` are the four pipeline stages described above.
- Plotting functions in `R/plot_utils.R` — `plot_moi()`,
  `plot_markers_diversity()`, `plot_allele_distribution()`,
  `plot_probability_histogram()`, `plot_comparison_heatmap()`, and
  `plot_likelihood_diagnostics()` (MCMC convergence diagnostics) — all
  return ggplot2 objects rather than drawing directly, so callers can
  inspect, recombine, or save them.
- `check_mcmc_diagnostics()` (`R/quality_control.R`) prints plain-language
  guidance when divergences, treedepth hits, or E-BFMI look problematic
  after a Stan run.

### Package-level bits

- `R/MalReBay-package.R` — package-level documentation and
  `utils::globalVariables()` declarations for the non-standard-evaluation
  column names used throughout the dplyr pipelines.
- `R/zzz.R` — `.onAttach()` warns once at load time if CmdStan isn't found.
- `R/utils-pipe.R` — re-exports the `%>%` pipe.

## Tests

Tests live in `tests/testthat/` and run with:

```r
devtools::test()
# or a single file:
testthat::test_file("tests/testthat/test_MalReBay.R")
```

Most tests use in-memory mock data; a few use a real Angola TES dataset
pre-imported and saved as `inst/extdata/imported_data.rds`, to avoid
`import_data()` file-path issues before the package is installed. MCMC-based
tests are slow (they run real Stan sampling) — expect the full suite to take
a few minutes.

## Documentation

Roxygen comments in `R/*.R` are the source of truth for `man/*.Rd` and
`NAMESPACE`. After changing any `#'` roxygen block, regenerate both with:

```r
devtools::document()
# or: roxygen2::roxygenise()
```

Do not hand-edit files under `man/` or `NAMESPACE` — they're marked
"Generated by roxygen2: do not edit by hand" for a reason; a hand edit will
just get silently overwritten (or drift out of sync) the next time someone
runs `devtools::document()`.

## The tutorial vignette

`vignettes/MalReBay.qmd` is the single source for the package tutorial —
it's a runnable Quarto notebook (open it directly in RStudio/Positron to
execute it) *and* the vignette that gets rendered onto the
[pkgdown site](https://swisstph.github.io/MalReBay/) as an Article. There is
no separate notebook copy anywhere else in the repo — edit this file
directly, and nowhere else, when updating the tutorial.

Building it requires the [Quarto CLI](https://quarto.org/docs/get-started/)
plus the `quarto` R package (`install.packages("quarto")`). To check it
still renders after an edit:

```r
quarto::quarto_render("vignettes/MalReBay.qmd")
# or, to also rebuild the full site:
pkgdown::build_site()
```

## Submitting changes

- Open an issue first for anything non-trivial, so the approach can be
  discussed before you invest time in it.
- Keep pull requests focused — one logical change per PR is easier to
  review than several bundled together.
- Make sure `devtools::check()` and `devtools::test()` pass locally before
  opening a PR; CI (`R-CMD-check.yaml`) runs the same checks on GitHub.
- If your change touches exported functions, update the relevant roxygen
  docs (see above) and, if it changes user-facing behaviour, the tutorial
  vignette.
