# Classify malaria infections using Bayesian MCMC

Runs the Bayesian MCMC engine across all sites to classify patient
samples as recrudescence or reinfection. Called automatically by
[`MalReBay`](https://swisstph.github.io/MalReBay/reference/MalReBay.md),
but can also be used directly for a step-by-step workflow.

## Usage

``` r
classify_infections(
  imported_data,
  mcmc_config = system.file("extdata", "default_mcmc_config.rds", package = "MalReBay"),
  n_workers = 1,
  verbose = TRUE,
  suppress_warnings = TRUE
)
```

## Arguments

- imported_data:

  A list returned by
  [`import_data`](https://swisstph.github.io/MalReBay/reference/import_data.md).

- mcmc_config:

  Path to an MCMC configuration Excel file, or a named list of
  parameters. Defaults to the bundled configuration.

- n_workers:

  Number of parallel workers. Defaults to `1`.

- verbose:

  Logical. Print progress messages. Defaults to `TRUE`.

- suppress_warnings:

  Logical. If `TRUE` (default), silences the MCMC sampler's own console
  warnings about divergent transitions, treedepth, and E-BFMI (the "N of
  M transitions ended with a divergence" messages cmdstanr prints as
  soon as sampling finishes, independently of `verbose`). Set to `FALSE`
  to see them – e.g. while tuning `adapt_delta` in `mcmc_config` for a
  problematic site. This only affects what gets printed; the same
  diagnostics are always available afterwards via
  [`summarise_results`](https://swisstph.github.io/MalReBay/reference/summarise_results.md)'s
  `convergence` table.

## Value

A named list of raw MCMC results per site, or `NULL` if no valid results
are produced. Pass this to
[`summarise_results`](https://swisstph.github.io/MalReBay/reference/summarise_results.md).

## See also

[`import_data`](https://swisstph.github.io/MalReBay/reference/import_data.md),
[`summarise_results`](https://swisstph.github.io/MalReBay/reference/summarise_results.md),
[`MalReBay`](https://swisstph.github.io/MalReBay/reference/MalReBay.md)

## Examples

``` r
if (FALSE) { # \dontrun{
imported <- import_data()
results  <- classify_infections(imported)
} # }
```
