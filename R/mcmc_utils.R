#' Run Bayesian MCMC Classification
#'
#' Runs the Bayesian MCMC engine across all sites in the imported data to
#' classify patient samples as recrudescence or reinfection. This is an
#' internal function called automatically by \code{\link{MalReBay}}.
#'
#' @param imported_data A list returned by \code{\link{import_data}}.
#' @param mcmc_config A named list of MCMC parameters. See
#'   \code{\link{MalReBay}} for valid entries and defaults.
#' @param n_workers An integer specifying the number of parallel workers.
#'   Defaults to \code{1} (sequential).
#' @param verbose Logical. If \code{TRUE}, prints progress messages.
#'   Defaults to \code{TRUE}.
#'
#' @return A named list of raw MCMC results per site, containing
#'   classifications, log-likelihoods, locus likelihood ratios, locus
#'   distances, locus summaries, patient IDs, and loci names. Returns
#'   \code{NULL} if no valid results are produced.
#'
#' @seealso \code{\link{summarise_results}} for post-processing the output.
#'
#' @keywords internal
classify_infections <- function(imported_data, mcmc_config = list(), n_workers = 1, verbose = TRUE) {
  
  # Validate imported_data structure
  required_elements <- c("late_failures", "additional", "marker_info", "data_type")
  if (!is.list(imported_data) || !all(required_elements %in% names(imported_data))) {
    stop(
      "'imported_data' must be a valid list object returned by the import_data() function.",
      call. = FALSE
    )
  }
  
  # Build config with defaults
  config <- utils::modifyList(list(
    n_chains              = 4,
    chunk_size            = 500,
    max_iterations        = 1000,
    burn_in_frac          = 0.25,
    rhat_threshold        = 1.1,
    ess_threshold         = 200,
    record_hidden_alleles = FALSE
  ), mcmc_config)
  
  # Extract data
  late_failures <- imported_data$late_failures
  additional    <- imported_data$additional
  marker_info   <- imported_data$marker_info
  data_type     <- imported_data$data_type
  
  # Parallelization setup
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  
  avail_cores    <- future::availableCores()
  workers_to_use <- min(n_workers, avail_cores)
  
  if (workers_to_use > 1) {
    if (verbose) message("INFO: Running in parallel on ", workers_to_use, " cores.")
    future::plan(future::multisession, workers = workers_to_use)
  } else {
    future::plan(future::sequential)
  }
  
  # Run MCMC
  results <- run_all_sites(
    late_failures = late_failures,
    additional    = additional,
    marker_info   = marker_info,
    mcmc_config   = config,
    data_type     = data_type,
    verbose       = verbose
  )
  
  if (is.null(results) || length(results$ids) == 0) {
    warning("MCMC results are empty.")
    return(NULL)
  }
  
  return(results)
}


#' Summarise MCMC Classification Results
#'
#' Post-processes raw MCMC output from \code{\link{classify_infections}} into
#' interpretable summaries including posterior probabilities, convergence
#' diagnostics, and a match counting comparison table. This is an internal
#' function called automatically by \code{\link{MalReBay}}.
#'
#' @param mcmc_results A list returned by \code{\link{classify_infections}}.
#' @param imported_data A list returned by \code{\link{import_data}}.
#' @param output_folder A string path to the directory where convergence
#'   diagnostic plots will be saved. If \code{NULL}, plots are not saved.
#'   Defaults to \code{NULL}.
#' @param verbose Logical. If \code{TRUE}, prints progress messages.
#'   Defaults to \code{TRUE}.
#'
#' @return A named list containing:
#'   \describe{
#'     \item{\code{posterior_probabilities}}{A data frame with per-patient
#'       posterior probabilities of recrudescence, and per-site locus
#'       availability summaries.}
#'     \item{\code{comparison}}{A data frame combining allele match counting
#'       results with posterior probabilities for each sample.}
#'     \item{\code{convergence}}{A data frame with per-site Gelman-Rubin
#'       R-hat and effective sample size diagnostics. \code{NULL} if
#'       diagnostics are unavailable.}
#'     \item{\code{mcmc_loglikelihoods}}{A named list of per-site MCMC
#'       log-likelihood chains.}
#'   }
#'
#' @seealso \code{\link{classify_infections}} for the MCMC engine,
#'   \code{\link{save_results}} for saving the output to CSV.
#'
#' @keywords internal
summarise_results <- function(mcmc_results, imported_data, output_folder = NULL, verbose = TRUE) {
  
  # Validate inputs
  if (is.null(mcmc_results) || length(mcmc_results$ids) == 0) {
    stop("'mcmc_results' is empty. Check classify_infections() ran successfully.", call. = FALSE)
  }
  
  required_elements <- c("late_failures", "additional", "marker_info", "data_type")
  if (!is.list(imported_data) || !all(required_elements %in% names(imported_data))) {
    stop("'imported_data' must be a valid list object returned by import_data().", call. = FALSE)
  }
  
  late_failures <- imported_data$late_failures
  marker_info   <- imported_data$marker_info
  
  # Create posterior probability summary and convergence diagnostics
  summary_list     <- list()
  convergence_list <- list()
  
  for (site in names(mcmc_results$ids)) {
    
    summary_list[[site]] <- data.frame(
      Site        = site,
      Sample.ID   = mcmc_results$ids[[site]],
      Probability = rowMeans(mcmc_results$classifications[[site]], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    
    # Record diagnostics of the loglikelihoods
    diag_vals <- NULL
    if (!is.null(mcmc_results$all_chains_loglikelihood[[site]])) {
      diag_vals <- plot_likelihood_diagnostics(
        mcmc_results$all_chains_loglikelihood[[site]],
        site,
        output_folder = output_folder,
        verbose       = verbose
      )
    }
    
    # Store R-hat and ESS into the convergence_list
    if (!is.null(diag_vals)) {
      convergence_list[[site]] <- data.frame(
        Site                  = site,
        Gelman_Rubin_Rhat     = if (!is.null(diag_vals$gelman)) round(diag_vals$gelman$psrf[1, 1], 4) else NA,
        Effective_Sample_Size = if (!is.null(diag_vals$ess))    round(as.numeric(diag_vals$ess), 2)    else NA,
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Build final summary with locus info
  posterior_probabilities <- dplyr::bind_rows(summary_list) %>%
    dplyr::left_join(
      dplyr::bind_rows(mcmc_results$locus_summary, .id = "Site") %>%
        dplyr::rename(
          Sample.ID         = patient_id,
          N_Available_D0    = n_available_d0,
          N_Available_DF    = n_available_df,
          N_Comparable_Loci = n_comparable_loci
        ),
      by = c("Sample.ID", "Site")
    )
  
  # Build comparison table
  match_results    <- perform_match_counting(late_failures, marker_info)
  comparison_table <- late_failures %>%
    dplyr::mutate(Base.ID = trimws(gsub(" Day 0| recurrence", "", Sample.ID))) %>%
    dplyr::left_join(match_results, by = c("Base.ID" = "Sample.ID")) %>%
    dplyr::left_join(
      dplyr::select(posterior_probabilities, Sample.ID, Probability, N_Comparable_Loci),
      by = c("Base.ID" = "Sample.ID")
    ) %>%
    dplyr::select(-Base.ID)
  
  # Build convergence summary
  convergence_summary <- if (length(convergence_list) > 0) {
    dplyr::bind_rows(convergence_list)
  } else {
    NULL
  }
  
  return(list(
    posterior_probabilities = posterior_probabilities,
    comparison              = comparison_table,
    convergence             = convergence_summary,
    mcmc_loglikelihoods     = mcmc_results$all_chains_loglikelihood
  ))
}

#' Save Classification Results to CSV
#'
#' Writes the summarised MCMC classification results to CSV files in the
#' specified output folder. This is an internal function called automatically
#' by \code{\link{MalReBay}}.
#'
#' The following files are written:
#' \describe{
#'   \item{\code{posterior_probabilities.csv}}{Per-patient posterior
#'     probabilities of recrudescence with locus availability summaries.}
#'   \item{\code{bayesian_match_counting_comparison.csv}}{Combined allele
#'     match counting and posterior probability results per sample.}
#'   \item{\code{mcmc_convergence_summary.csv}}{Per-site Gelman-Rubin R-hat
#'     and effective sample size diagnostics. Only written if convergence
#'     diagnostics are available.}
#' }
#'
#' @param summary_results A list returned by \code{\link{summarise_results}}.
#' @param output_folder A string path to the directory where CSV files will
#'   be saved. Created recursively if it does not exist. Defaults to
#'   \code{"results"}.
#' @param verbose Logical. If \code{TRUE}, prints a message confirming where
#'   files were saved. Defaults to \code{TRUE}.
#'
#' @return Invisibly returns a named character vector of file paths written:
#'   \describe{
#'     \item{\code{posterior_probabilities}}{Path to the posterior
#'       probabilities CSV.}
#'     \item{\code{comparison}}{Path to the match counting comparison CSV.}
#'     \item{\code{convergence}}{Path to the convergence summary CSV.
#'       Only present if convergence data is available.}
#'   }
#'
#' @seealso \code{\link{summarise_results}} for generating the input to this
#'   function.
#'
#' @keywords internal
save_results <- function(summary_results, output_folder = "results", verbose = TRUE) {
  
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  
  if (!is.list(summary_results) ||
      !all(c("posterior_probabilities", "comparison", "convergence") %in% names(summary_results))) {
    stop("'summary_results' must be a valid list returned by summarise_results().", call. = FALSE)
  }
  
  paths <- c()
  
  # Posterior probabilities
  path_pp <- file.path(output_folder, "posterior_probabilities.csv")
  utils::write.csv(summary_results$posterior_probabilities, path_pp, row.names = FALSE)
  paths["posterior_probabilities"] <- path_pp
  
  # Comparison table
  path_ct <- file.path(output_folder, "bayesian_match_counting_comparison.csv")
  utils::write.csv(summary_results$comparison, path_ct, row.names = FALSE)
  paths["comparison"] <- path_ct
  
  # Convergence summary
  if (!is.null(summary_results$convergence)) {
    path_cv <- file.path(output_folder, "mcmc_convergence_summary.csv")
    utils::write.csv(summary_results$convergence, path_cv, row.names = FALSE)
    paths["convergence"] <- path_cv
  }
  
  if (verbose) message("INFO: Results saved to: ", output_folder)
  
  invisible(paths)
}