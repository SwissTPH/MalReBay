#' Run the MalReBay malaria recrudescence classification pipeline
#'
#' Imports genotype data, runs Bayesian MCMC classification to distinguish
#' recrudescent from reinfection malaria treatment failures, summarises
#' results, and saves output plots and tables.
#'
#' @param filepath        Path to the genotype data Excel file.
#'                        Defaults to the package example dataset.
#' @param marker_filepath Path to the marker metadata Excel file.
#'                        Defaults to the package example marker file.
#' @param mcmc_config     Path to the MCMC configuration Excel file.
#'                        Defaults to the package default configuration.
#' @param output_folder   Path to a folder for saving results and plots.
#'                        Set to \code{NULL} to skip saving (results are
#'                        still returned invisibly).
#' @param n_workers       Number of parallel workers. For
#'                        length-polymorphic data this is handled
#'                        automatically by the sampler; the argument is
#'                        retained for compatibility with amplicon-sequencing
#'                        data.
#' @param verbose         If \code{TRUE}, print progress messages to the
#'                        console.
#'
#' @return A list of per-site summary results (invisibly). See
#'   \code{summarise_results()} for details of the list structure.
#'
#' @examples
#' \dontrun{
#' # Run on the bundled example data with default settings
#' results <- MalReBay()
#'
#' # Run on your own data and save output to a folder
#' results <- MalReBay(
#'   filepath      = "path/to/your_data.xlsx",
#'   output_folder = "path/to/output"
#' )
#' }
#'
#' @export
MalReBay <- function(
    filepath        = system.file("extdata", "Angola_2021_TES_7NMS.xlsx",
                                  package = "MalReBay"),
    marker_filepath = system.file("extdata", "makers_details.xlsx",
                                  package = "MalReBay"),
    mcmc_config     = system.file("extdata", "default_mcmc_config.xlsx",
                                  package = "MalReBay"),
    output_folder   = NULL,
    n_workers       = 1,
    verbose         = TRUE
) {

  if (verbose) message("Starting MalReBay pipeline...")

  imported_data <- import_data(
    filepath        = filepath,
    marker_filepath = marker_filepath,
    verbose         = verbose
  )

  mcmc_results <- classify_infections(
    imported_data = imported_data,
    mcmc_config   = mcmc_config,
    n_workers     = n_workers,
    verbose       = verbose
  )

  summary_results <- summarise_results(
    mcmc_results  = mcmc_results,
    imported_data = imported_data,
    output_folder = output_folder,
    verbose       = verbose
  )

  save_results(
    summary_results = summary_results,
    imported_data   = imported_data,
    output_folder   = output_folder,
    verbose         = verbose
  )

  invisible(summary_results)
}
