# Run the MalReBay malaria recrudescence classification pipeline

Imports genotype data, runs Bayesian MCMC classification to distinguish
recrudescent from reinfection malaria treatment failures, summarises
results, and saves output plots and tables.

## Usage

``` r
MalReBay(
  filepath = system.file("extdata", "Dataset_microsatellite_panel.xlsx", package = "MalReBay"),
  marker_filepath = system.file("extdata", "makers_details.xlsx", package = "MalReBay"),
  additional_filepath = NULL,
  mcmc_config = system.file("extdata", "default_mcmc_config.rds", package = "MalReBay"),
  output_folder = NULL,
  n_workers = 1,
  verbose = TRUE,
  suppress_warnings = TRUE,
  prob_threshold = 0.5,
  plots = TRUE
)
```

## Arguments

- filepath:

  Path to the genotype data Excel file. Defaults to the package example
  dataset.

- marker_filepath:

  Path to the marker metadata Excel file. Defaults to the package
  example marker file.

- mcmc_config:

  Path to the MCMC configuration Excel file. Defaults to the package
  default configuration.

- output_folder:

  Path to a folder for saving results and plots. Set to `NULL` to skip
  saving (results are still returned invisibly).

- n_workers:

  Number of parallel workers. For length-polymorphic data this is
  handled automatically by the sampler; the argument is retained for
  compatibility with amplicon-sequencing data.

- verbose:

  If `TRUE`, print progress messages to the console.

- suppress_warnings:

  Logical. If `TRUE` (default), silences the MCMC sampler's own console
  warnings about divergent transitions, treedepth, and E-BFMI. See
  [`classify_infections`](https://swisstph.github.io/MalReBay/reference/classify_infections.md)
  for details – these diagnostics remain available afterwards via
  `summary_results$convergence` regardless of this setting.

- prob_threshold:

  Numeric in `[0, 1]`. The posterior probability at or above which a
  recurrence is classified as `"Recrudescence"` (below it,
  `"Reinfection"`). Defaults to `0.5`, the natural cutoff under this
  model's equal 50/50 prior. See
  [`summarise_results`](https://swisstph.github.io/MalReBay/reference/summarise_results.md)
  for where this is applied.

- plots:

  Logical. If `FALSE`, skips all plots (descriptive, result and
  convergence diagnostic plots); convergence diagnostics are still
  computed and printed. Defaults to `TRUE`.

## Value

A list of per-site summary results (invisibly). See
[`summarise_results()`](https://swisstph.github.io/MalReBay/reference/summarise_results.md)
for details of the list structure.

## Examples

``` r
if (FALSE) { # \dontrun{
# Run on the bundled example data with default settings
results <- MalReBay()

# Run on your own data and save output to a folder
results <- MalReBay(
  filepath      = "path/to/your_data.xlsx",
  output_folder = "path/to/output"
)
} # }
```
