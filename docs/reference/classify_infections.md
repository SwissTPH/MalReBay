# Classify malaria infections using Bayesian MCMC

Classify malaria infections using Bayesian MCMC

## Usage

``` r
classify_infections(
  imported_data,
  mcmc_config = system.file("extdata", "default_mcmc_config.xlsx", package = "MalReBay"),
  n_workers = 1,
  verbose = TRUE
)
```

## Arguments

- imported_data:

  List returned by
  [`import_data()`](https://swisstph.github.io/MalReBay/reference/import_data.md).

- mcmc_config:

  Path to MCMC configuration Excel file.

- n_workers:

  Parallel workers. Ignored for length-polymorphic data (Stan handles
  parallelism internally). Used for ampseq.

- verbose:

  Print progress messages.

## Value

List of MCMC results passed to
[`summarise_results()`](https://swisstph.github.io/MalReBay/reference/summarise_results.md).
