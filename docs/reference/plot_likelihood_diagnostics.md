# Plot MCMC Likelihood Diagnostics

Calculates standard MCMC convergence diagnostics and builds four ggplot2
diagnostic panels: Gelman-Rubin shrink factor, Traceplot, Log-Posterior
Histogram, and Autocorrelation. The panels are returned as ggplot
objects (`$traceplot`, `$gelman_rubin`, `$log_posterior`,
`$autocorrelation`) alongside the numeric diagnostics, so they can be
printed, saved, or recombined like any other plot in this package –
there's no need to save to PNG and reload just to view them (e.g. in a
notebook).

## Usage

``` r
plot_likelihood_diagnostics(
  all_chains_loglikelihood = NULL,
  site_name,
  stan_fit = NULL,
  save_plot = TRUE,
  output_folder = NULL,
  verbose = TRUE,
  combine_plots = FALSE
)
```

## Arguments

- all_chains_loglikelihood:

  A list where each element is a numeric vector representing the
  log-likelihood history of one MCMC chain.

- site_name:

  A character string for labeling plots.

- stan_fit:

  An optional `CmdStanMCMC` object (from cmdstanr) for additional
  diagnostics.

- save_plot:

  A logical. If `TRUE`, panels are also saved as PNG files.

- output_folder:

  A character string specifying the path to save plots.

- verbose:

  A logical. If `TRUE`, prints diagnostic summaries.

- combine_plots:

  A logical. If `TRUE`, the returned list also includes `$combined`: all
  four panels arranged in a single 2x2
  [`patchwork::wrap_plots()`](https://patchwork.data-imaginist.com/reference/wrap_plots.html)
  grid.

## Value

An invisible list with the numeric diagnostics (`gelman`, `ess`,
`rhat_rank`, `ess_bulk`, `ess_tail`, `geweke`) and the ggplot2 panels
(`traceplot`, `gelman_rubin` – `NULL` with a single chain,
`log_posterior`, `autocorrelation`, and `combined` when
`combine_plots = TRUE`).

## Examples

``` r
if (FALSE) { # \dontrun{
  chains <- list(
    c(-10.2, -9.8, -10.5, -9.9, -10.1),
    c(-9.9,  -10.3, -10.0, -9.7, -10.4)
  )
  diag <- plot_likelihood_diagnostics(
    all_chains_loglikelihood = chains,
    site_name  = "TestSite",
    save_plot  = FALSE,
    verbose    = FALSE
  )
  diag$traceplot
} # }
```
