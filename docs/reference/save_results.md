# Save Classification Results and Generate Plots

Writes MCMC results to CSV files and generates descriptive and
result-based plots (Diversity, MOI, and Posterior Histograms).

## Usage

``` r
save_results(
  summary_results,
  imported_data = NULL,
  output_folder = NULL,
  verbose = TRUE,
  plots = TRUE
)
```

## Arguments

- summary_results:

  List from `summarise_results`.

- imported_data:

  List from `import_data`. When `NULL`, data-dependent plots (diversity,
  MOI) are skipped.

- output_folder:

  Path to save files. `NULL` prints plots to the R plot window and skips
  CSV saving.

- verbose:

  Logical.

- plots:

  Logical. If `FALSE`, no plots are generated, shown or saved; only the
  CSV files are written. Defaults to `TRUE`.

## Value

A named character vector of paths to the files that were written, or
`invisible(NULL)` when `output_folder` is `NULL`.

## See also

[`summarise_results`](https://swisstph.github.io/MalReBay/reference/summarise_results.md),
[`MalReBay`](https://swisstph.github.io/MalReBay/reference/MalReBay.md)

## Examples

``` r
if (FALSE) { # \dontrun{
imported <- import_data()
results  <- classify_infections(imported)
summary  <- summarise_results(results, imported)
save_results(summary, imported, output_folder = "my_results")
} # }
```
