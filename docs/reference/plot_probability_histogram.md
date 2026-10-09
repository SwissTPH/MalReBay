# Plot Posterior Probability Histogram

Creates and saves a histogram of posterior probabilities of
recrudescence for all patients. This is an internal function called
automatically by
[`MalReBay`](https://swisstph.github.io/MalReBay/reference/MalReBay.md).

## Usage

``` r
plot_probability_histogram(
  summary_results,
  output_folder = NULL,
  verbose = TRUE
)
```

## Arguments

- summary_results:

  A list returned by
  [`summarise_results`](https://swisstph.github.io/MalReBay/reference/summarise_results.md),
  containing a `posterior_probabilities` data frame with a `Probability`
  column.

- output_folder:

  A string specifying the directory where the histogram PNG will be
  saved. Defaults to `"results"`.

- verbose:

  Logical. If `TRUE`, prints a message when the file is saved. Defaults
  to `TRUE`.

## Value

When `output_folder` is `NULL`, displays the histogram on the current
graphics device and returns `invisible(NULL)`. When `output_folder` is
provided, saves a PNG and returns the file path invisibly. Returns
`invisible(NULL)` if no probabilities are available.

## Examples

``` r
summary_results <- list(
  posterior_probabilities = data.frame(
    Patient_ID  = c("P1", "P2", "P3"),
    Probability = c(0.12, 0.87, 0.45)
  )
)
plot_probability_histogram(summary_results)
#> INFO: No 'Site' column in posterior_probabilities; skipping per-site histogram.

```
