# Summarise MCMC Classification Results

Post-processes raw MCMC output from
[`classify_infections`](https://swisstph.github.io/MalReBay/reference/classify_infections.md)
into interpretable summaries including posterior probabilities,
convergence diagnostics, and a match counting comparison table. Called
automatically by
[`MalReBay`](https://swisstph.github.io/MalReBay/reference/MalReBay.md).

## Usage

``` r
summarise_results(
  mcmc_results,
  imported_data,
  output_folder = NULL,
  verbose = TRUE,
  prob_threshold = 0.5,
  plots = TRUE
)
```

## Arguments

- mcmc_results:

  A list returned by
  [`classify_infections`](https://swisstph.github.io/MalReBay/reference/classify_infections.md).

- imported_data:

  A list returned by
  [`import_data`](https://swisstph.github.io/MalReBay/reference/import_data.md).

- output_folder:

  Path for saving convergence diagnostic plots. `NULL` skips saving.

- verbose:

  Logical. Print progress messages.

- prob_threshold:

  Numeric in `[0, 1]`. The posterior probability at or above which a
  recurrence is classified as `"Recrudescence"` (below it,
  `"Reinfection"`); used to add the `Classification` column to
  `posterior_probabilities` (and `comparison`). Defaults to `0.5`, the
  natural cutoff under this model's equal 50/50 prior (see
  @sec-interpretation-probability in the analysis notebook).

- plots:

  Logical. If `FALSE`, convergence diagnostic plots are not saved to
  `output_folder`; the diagnostics themselves are still computed and
  returned in `convergence`. Defaults to `TRUE`.

## Value

A named list with `posterior_probabilities`, `comparison`,
`convergence`, and `mcmc_loglikelihoods`. Both `posterior_probabilities`
and `comparison` include a `Classification` column derived from
`prob_threshold`.

## See also

[`classify_infections`](https://swisstph.github.io/MalReBay/reference/classify_infections.md),
[`save_results`](https://swisstph.github.io/MalReBay/reference/save_results.md),
[`MalReBay`](https://swisstph.github.io/MalReBay/reference/MalReBay.md)

## Examples

``` r
if (FALSE) { # \dontrun{
imported <- import_data()
results  <- classify_infections(imported)
summary  <- summarise_results(results, imported)
} # }
```
