# Classify Infections Using a Bayesian MCMC Framework

This is the main function of the `MalReBay` package. It runs the
complete Bayesian analysis workflow to classify parasite infections as
either recrudescence or reinfection based on genotyping data. The
function takes a processed data object from
[`import_data()`](https://swisstph.github.io/MalReBay/reference/import_data.md),
runs MCMC simulation with automatic convergence checking, and summarizes
the final results.

## Usage

``` r
classify_infections(
  imported_data,
  mcmc_config,
  output_folder,
  n_workers = 1,
  verbose = TRUE
)
```

## Arguments

- imported_data:

  A list object returned by the
  [`import_data()`](https://swisstph.github.io/MalReBay/reference/import_data.md)
  function. This list must contain `$late_failures` (data.frame),
  `$additional` (data.frame), `$marker_info` (data.frame), and
  `$data_type` (character string).

- mcmc_config:

  A list of MCMC configuration parameters. See Details for a full list
  of options and their defaults.

- output_folder:

  A character string specifying the path to the folder where all results
  (summary CSVs, diagnostic plots) will be saved. The folder will be
  created if it does not exist. **Note:** If a `convergence_diagnosis`
  sub-folder exists from a previous run, it will be automatically
  removed to ensure that all diagnostic plots are from the current
  analysis.

- n_workers:

  An integer specifying the number of parallel cores to use. Defaults to
  1 (sequential). If `n_workers > 1`, the analysis will be parallelized
  by site using the 'future' package.

- verbose:

  A logical. If `TRUE` (the default), the function will print progress
  messages, warnings about missing markers, and other information to the
  console. If `FALSE`, it will run silently, only showing critical
  errors.

## Value

A list containing three elements:

- summary:

  A data frame summarizing the main results for each patient...

- marker_details:

  A detailed data frame providing the mean likelihood ratio...

- mcmc_loglikelihoods:

  A list containing the raw MCMC log-likelihood chains for each site,
  which can be passed to
  [`generate_likelihood_diagnostics()`](https://swisstph.github.io/MalReBay/reference/generate_likelihood_diagnostics.md).

CSV files and diagnostic plots are also saved to the `output_folder`.

## Details

The function operates by first importing and validating the data from
the provided Excel file. It automatically detects the data type
(length-polymorphic or amplicon sequencing) and selects the appropriate
MCMC engine. The MCMC simulation is run in parallel across multiple
chains and proceeds in chunks, stopping automatically when convergence
criteria are met or the maximum number of iterations is reached.

The `mcmc_config` list provides fine-grained control over the
simulation. Key parameters include:

- **`n_chains`**: Number of parallel chains to run. (Default: 4).

- **`chunk_size`**: Number of iterations per chunk. After each chunk,
  convergence is assessed. (Default: 10000).

- **`max_iterations`**: The maximum total number of iterations before
  the simulation is forcibly stopped, even if not converged. (Default:
  100000).

- **`burn_in_frac`**: The fraction of initial samples to discard from
  each chain before summarizing results. (Default: 0.25).

- **`rhat_threshold`**: The Gelman-Rubin diagnostic threshold for
  convergence. (Default: 1.01).

- **`ess_threshold`**: The Effective Sample Size threshold for
  convergence. (Default: 400).

- **`record_hidden_alleles`**: A logical flag. If `TRUE`, the full state
  of imputed hidden alleles is saved. This can generate very large
  output files and is mainly for debugging. (Default: FALSE).

Any parameters not specified in the list will use the default values.

## Examples

``` r
# 1. Get the path to the example data file
example_file <- system.file("extdata", "Angola_2021_TES_7NMS.xlsx",
                            package = "MalReBay")

# 2. Create a temporary directory for the results
temp_output_dir <- file.path(tempdir(), "MalReBay_example")
dir.create(temp_output_dir, showWarnings = FALSE)

# 3. Define a minimal MCMC configuration for a quick run
# NOTE: For a real analysis, use more iterations.
quick_mcmc_config <- list(
  n_chains = 2,
  max_iterations = 500,
  chunk_size = 500
)

if (FALSE) { # \dontrun{
# 4. Import and process the data first
processed_data <- import_data(filepath = example_file)

# 5. Run the analysis using the processed data object and 2 cores
results_list <- classify_infections(
  imported_data = processed_data,
  mcmc_config = quick_mcmc_config,
  output_folder = temp_output_dir,
  n_workers = 2,
  verbose = TRUE
)

# 6. View the top rows of the main summary table
print(head(results_list$summary))
} # }
```
