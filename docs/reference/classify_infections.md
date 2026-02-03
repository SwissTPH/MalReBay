# Classify Infections Using a Bayesian MCMC Framework

This is the main function of the `MalReBay` package. It runs the
complete Bayesian analysis workflow to classify parasite infections as
either recrudescence or reinfection based on genotyping data. The
function takes a processed data object, runs MCMC simulation with
automatic convergence checking, and summarizes the final results.

## Usage

``` r
classify_infections(
  imported_data,
  mcmc_config = list(),
  output_folder = "results",
  n_workers = 1,
  verbose = TRUE
)
```

## Arguments

- imported_data:

  A list object (typically from an import helper) containing:

  - `late_failures`: Data frame of paired Day 0/Day Failure samples.

  - `additional`: Data frame of background/population samples.

  - `marker_info`: Data frame with marker names and repeat lengths.

  - `data_type`: Character, either "length_polymorphic" or "ampseq".

- mcmc_config:

  A list of MCMC configuration parameters. See Details.

- output_folder:

  A character string specifying the path to the folder where all results
  (summary CSVs, diagnostic plots) will be saved. The folder will be
  created if it does not exist. If a `convergence_diagnosis` sub-folder
  exists from a previous run, it will be automatically updated.

- n_workers:

  An integer specifying the number of parallel cores to use. Defaults
  to 1. If `n_workers > 1`, the analysis is parallelized by site.

- verbose:

  A logical. If `TRUE`, the function prints progress messages.

## Value

A list containing three elements:

- summary:

  A data frame summarizing the classification probability for each
  patient.

- marker_details:

  A detailed data frame providing the mean likelihood ratio and distance
  per marker.

- comparison:

  A comparison table matching Bayesian results with basic match-counting
  logic.

## Details

The function operates by validating the data from the provided input
object. It automatically detects the data type (length-polymorphic or
amplicon sequencing) and selects the appropriate MCMC engine. The MCMC
simulation is run in parallel across multiple chains and proceeds in
chunks, stopping automatically when convergence criteria are met or the
maximum number of iterations is reached.

The `mcmc_config` list provides fine-grained control over the
simulation. Key parameters include:

- **`n_chains`**: Number of parallel chains to run. (Default: 4).

- **`chunk_size`**: Number of iterations per chunk. After each chunk,
  convergence is assessed. (Default: 1000).

- **`max_iterations`**: The maximum total number of iterations before
  the simulation is forcibly stopped. (Default: 10000).

- **`burn_in_frac`**: The fraction of initial samples to discard from
  each chain before summarizing results. (Default: 0.25).

- **`record_hidden_alleles`**: A logical flag. If `TRUE`, the full state
  of imputed hidden alleles is saved. (Default: FALSE).

## Examples

``` r
if (FALSE) { # \dontrun{
# 1. Import and process the data
processed_data <- import_data(filepath = "my_data.xlsx")

# 2. Run the analysis
results <- classify_infections(
  imported_data = processed_data,
  n_workers = 2
)
} # }
```
