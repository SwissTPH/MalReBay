#' Run the Full MalReBay Pipeline
#'
#' A wrapper function that runs the complete Bayesian malaria recrudescence
#' classification pipeline from data import through to saved outputs. When
#' called without arguments, it runs on the built-in example dataset, making
#' it useful for verifying that the package is installed and working correctly.
#'
#' The pipeline executes the following steps in order:
#' \enumerate{
#'   \item Import and validate genotype data via \code{\link{import_data}}
#'   \item Generate descriptive plots via \code{\link{plot_markers_diversity}}
#'     and \code{\link{plot_moi}}
#'   \item Run the Bayesian MCMC classification via
#'     \code{\link{classify_infections}}
#'   \item Summarise MCMC results via \code{\link{summarise_results}}
#'   \item Plot posterior probability histogram
#'   \item Save all outputs to CSV via \code{\link{save_results}}
#' }
#'
#' @param filepath A string path to the input Excel file containing genotype
#'   data. Defaults to the built-in example dataset included with the package.
#' @param marker_filepath A string path to the marker information Excel file.
#'   If \code{NULL}, the function looks for a \code{marker_info} sheet within
#'   \code{filepath}. Defaults to the built-in example marker file.
#' @param mcmc_config A named list of MCMC configuration parameters. Valid
#'   entries are:
#'   \describe{
#'     \item{\code{n_chains}}{Number of MCMC chains. Default \code{2}.}
#'     \item{\code{chunk_size}}{Iterations per convergence check. Default
#'       \code{500}.}
#'     \item{\code{max_iterations}}{Maximum total iterations. Default
#'       \code{1000}.}
#'     \item{\code{burn_in_frac}}{Fraction of samples discarded as burn-in.
#'       Default \code{0.25}.}
#'     \item{\code{rhat_threshold}}{Gelman-Rubin convergence threshold.
#'       Default \code{1.1}.}
#'     \item{\code{ess_threshold}}{Minimum effective sample size for
#'       convergence. Default \code{200}.}
#'     \item{\code{record_hidden_alleles}}{Logical. Whether to record imputed
#'       hidden alleles. Default \code{FALSE}.}
#'   }
#' @param output_folder A string path to the directory where all output files
#'   will be saved. Defaults to a temporary directory via \code{tempdir()}.
#' @param n_workers An integer specifying the number of parallel workers for
#'   MCMC. Defaults to \code{1} (sequential). Set higher to use multiple
#'   cores.
#' @param verbose Logical. If \code{TRUE}, prints progress messages at each
#'   pipeline step. Defaults to \code{TRUE}.
#'
#' @return Invisibly returns a named list containing:
#'   \describe{
#'     \item{\code{posterior_probabilities}}{A data frame with per-patient
#'       posterior probabilities of recrudescence and locus-level summaries.}
#'     \item{\code{comparison}}{A data frame combining match counting results
#'       with posterior probabilities for each patient.}
#'     \item{\code{convergence}}{A data frame with per-site Gelman-Rubin
#'       R-hat and effective sample size diagnostics. \code{NULL} if
#'       diagnostics are unavailable.}
#'     \item{\code{mcmc_loglikelihoods}}{A named list of per-site MCMC
#'       log-likelihood chains for further diagnostics if needed.}
#'   }
#'
#' @examples
#' \donttest{
#'   # Run on built-in example data
#'   results <- MalReBay()
#'
#'   # Run on your own data
#'   results <- MalReBay(
#'     filepath        = "path/to/your/data.xlsx",
#'     marker_filepath = "path/to/your/markers.xlsx",
#'     mcmc_config     = list(n_chains = 4, max_iterations = 5000),
#'     output_folder   = "my_results",
#'     n_workers       = 4
#'   )
#'
#'   # Access results
#'   head(results$posterior_probabilities)
#' }
#'
#' @seealso
#'   \code{\link{import_data}} for data import and validation,
#'   \code{\link{classify_infections}} for the MCMC engine,
#'   \code{\link{summarise_results}} for post-processing,
#'   \code{\link{save_results}} for saving outputs,
#'   \code{\link{plot_markers_diversity}} and \code{\link{plot_moi}} for
#'   descriptive plots.
#'
#' @export
MalReBay <- function(
    filepath        = system.file("extdata", "Angola_2021_TES_7NMS.xlsx",    package = "MalReBay"),
    marker_filepath = system.file("extdata", "makers_details.xlsx", package = "MalReBay"),
    mcmc_config     = list(n_chains = 2, chunk_size = 500, max_iterations = 1000),
    output_folder   = tempdir(),
    n_workers       = 1,
    verbose         = TRUE
) {

  if (verbose) message("Running analysis")  
  # Import data
  imported_data <- import_data(
    filepath        = filepath,
    marker_filepath = marker_filepath,
    verbose         = verbose
  )
  
  # Run descriptive plots
  all_data <- dplyr::bind_rows(imported_data$late_failures, imported_data$additional)
  if (nrow(all_data) > 0) {
    plot_markers_diversity(all_data, imported_data$data_type, imported_data$marker_info, output_folder = output_folder)
    plot_moi(all_data, output_folder = output_folder)
  }
  
  # MCMC
  mcmc_results <- classify_infections(
    imported_data = imported_data,
    mcmc_config   = mcmc_config,
    n_workers     = n_workers,
    verbose       = verbose
  )
  
  # Summarise MCMC results
  summary_results <- summarise_results(
    mcmc_results  = mcmc_results,
    imported_data = imported_data,
    output_folder = output_folder,
    verbose       = verbose
  )
  
  # PLot the posterior probabilities
  plot_probability_histogram(summary_results, output_folder = output_folder, verbose = verbose)
  
  # Save the results
  save_results(summary_results, output_folder = output_folder, verbose = verbose) 
  if (verbose) message("Analysis is completed and results are save in", output_folder)

  invisible(summary_results)
}