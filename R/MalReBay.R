#' Run the Full MalReBay Pipeline
#'
#' A wrapper function that runs the complete Bayesian malaria recrudescence
#' classification pipeline. It handles data import, MCMC processing, 
#' results summarization, and output generation (plots and CSVs).
#'
#' The pipeline executes the following steps in order:
#' \enumerate{
#'   \item Import and validate genotype data via \code{\link{import_data}}
#'   \item Run the Bayesian MCMC classification via \code{\link{classify_infections}}
#'   \item Summarise MCMC results via \code{\link{summarise_results}}
#'   \item Save outputs and generate plots via \code{\link{save_results}}
#' }
#'
#' @param filepath A string path to the input Excel file containing genotype
#'   data. Defaults to the built-in example dataset.
#' @param marker_filepath A string path to the marker information Excel file.
#' @param mcmc_config A path to a configuration Excel file or a named list of parameters.
#' @param output_folder A string path to the directory where all output files
#'   (CSVs and Plots) will be saved. If \code{NULL}, results are not saved.
#' @param n_workers An integer specifying the number of parallel workers.
#' @param verbose Logical. If \code{TRUE}, prints progress messages.
#'
#' @return Invisibly returns a named list containing:
#'   \describe{
#'     \item{\code{posterior_probabilities}}{Data frame with results and locus summaries.}
#'     \item{\code{comparison}}{Data frame combining match counting and Bayesian results.}
#'     \item{\code{convergence}}{Diagnostics (R-hat and ESS).}
#'     \item{\code{mcmc_loglikelihoods}}{Raw MCMC chains.}
#'     \item{\code{chain_results}}{Full posterior samples.}
#'   }
#'
#' @export
MalReBay <- function(
    filepath        = system.file("extdata", "Angola_2021_TES_7NMS.xlsx",    package = "MalReBay"),
    marker_filepath = system.file("extdata", "makers_details.xlsx", package = "MalReBay"),
    mcmc_config     = system.file("extdata", "default_mcmc_config.xlsx", package = "MalReBay"),
    output_folder   = NULL,
    n_workers       = 1,
    verbose         = TRUE
) {
  
  if (verbose) message("Starting MalReBay Pipeline...")  
  
  # Import data
  imported_data <- import_data(
    filepath        = filepath,
    marker_filepath = marker_filepath,
    verbose         = verbose
  )
  
  # MCMC Classification
  mcmc_results <- classify_infections(
    imported_data = imported_data,
    mcmc_config   = mcmc_config,
    n_workers     = n_workers,
    verbose       = verbose
  )
  
  # Summarise results
  summary_results <- summarise_results(
    mcmc_results  = mcmc_results,
    imported_data = imported_data,
    output_folder = output_folder,
    verbose       = verbose
  )
  
  # Save results and plots
  save_results(
    summary_results = summary_results, 
    imported_data   = imported_data, 
    output_folder   = output_folder, 
    verbose         = verbose
  ) 
  
  invisible(summary_results)
}