# Generate and Plot MCMC Convergence Diagnostics

Calculates standard MCMC convergence diagnostics and generates
diagnostic plots in a multi-stage process: 1) Gelman-Rubin, 2) Traceplot
& Histogram, 3) All ACF plots.

## Usage

``` r
generate_likelihood_diagnostics(
  all_chains_loglikelihood,
  site_name,
  save_plot = FALSE,
  output_folder = NULL,
  verbose = TRUE
)
```

## Arguments

- all_chains_loglikelihood:

  A list where each element is a numeric vector representing the
  log-likelihood history of one MCMC chain.

- site_name:

  A character string for labeling plots.

- save_plot:

  A logical. If `TRUE`, plots are saved as PNG files.

- output_folder:

  A character string specifying the path to save plots.

- verbose:

  A logical. If `TRUE`, prints diagnostic summaries.

## Value

An invisible list containing the calculated diagnostic results.
