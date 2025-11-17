# Generate and Plot MCMC Convergence Diagnostics

Calculates standard MCMC convergence diagnostics (Gelman-Rubin R-hat,
Effective Sample Size) and generates diagnostic plots (traceplot, Gelman
plot, histogram, ACF) for the log-likelihood.

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

  A character string for labeling plots and creating the output
  subdirectory when saving.

- save_plot:

  A logical. If `TRUE`, plots are saved as PNG files. If `FALSE` (the
  default), plots are displayed in the active R graphics device.

- output_folder:

  A character string specifying the path to a directory where plots
  should be saved. Only used if `save_plot = TRUE`.

- verbose:

  A logical. If `TRUE` (the default), prints the Gelman-Rubin and ESS
  summaries to the console.

## Value

An invisible list containing the calculated diagnostic results:

- `gelman`: The Gelman-Rubin diagnostic result.

- `ess`: The Effective Sample Size result.

- `n_valid`: The number of valid (finite) samples per chain.

## Details

This function takes the log-likelihood history from multiple MCMC
chains, cleans the data, and uses the `coda` package to compute
diagnostics. It can either display plots in the active R graphics device
or save them as PNG files to a specified directory.

## Examples

``` r
# --- Example 1: Interactive Plotting ---
# To view plots in the RStudio Plots pane (will ask before showing each plot)
if (FALSE) { # \dontrun{
  generate_likelihood_diagnostics(all_chains, site_name = "Test Site")
} # }

# --- Example 2: Saving Plots to a File ---
# Create a temporary directory for the output
temp_dir <- tempdir()
# Simulate some dummy MCMC chains for demonstration
all_chains <- list(
  chain1 = rnorm(1000, mean = -150, sd = 5),
  chain2 = rnorm(1000, mean = -150, sd = 5)
)
# Run the function to save plots
generate_likelihood_diagnostics(all_chains,
                                site_name = "Test Site",
                                save_plot = TRUE,
                                output_folder = temp_dir)
#> --- Convergence Diagnostics for: Test Site ---
#> Gelman-Rubin R-hat:
#> Potential scale reduction factors:
#> 
#>               Point est. Upper C.I.
#> loglikelihood          1          1
#> 
#> 
#> Effective Sample Size (ESS):
#> loglikelihood 
#>          2000 

# Check that the files were created
list.files(file.path(temp_dir, "convergence_diagnosis", "Test Site"))
#> [1] "gelman_loglikelihood.png"    "loglikelihood_acf.png"      
#> [3] "loglikelihood_histogram.png" "loglikelihood_traceplot.png"
```
