#' Classify malaria infections using Bayesian MCMC
#'
#' Runs the Bayesian MCMC engine across all sites to classify patient samples
#' as recrudescence or reinfection. Called automatically by
#' \code{\link{MalReBay}}, but can also be used directly for a step-by-step
#' workflow.
#'
#' @param imported_data A list returned by \code{\link{import_data}}.
#' @param mcmc_config   Path to an MCMC configuration Excel file, or a named
#'   list of parameters. Defaults to the bundled configuration.
#' @param n_workers     Number of parallel workers. Defaults to \code{1}.
#' @param verbose       Logical. Print progress messages. Defaults to \code{TRUE}.
#'
#' @return A named list of raw MCMC results per site, or \code{NULL} if no
#'   valid results are produced. Pass this to \code{\link{summarise_results}}.
#'
#' @seealso \code{\link{import_data}}, \code{\link{summarise_results}},
#'   \code{\link{MalReBay}}
#'
#' @examples
#' \dontrun{
#' imported <- import_data()
#' results  <- classify_infections(imported)
#' }
#'
#' @export
classify_infections <- function(imported_data,
                                mcmc_config = system.file(
                                                  "extdata", 
                                                  "default_mcmc_config.xlsx",
                                                  package = "MalReBay"),
                                n_workers = 1,
                                verbose   = TRUE) {
  
  required_elements <- c("late_failures", "additional", "marker_info", "data_type")
  if (!is.list(imported_data) ||
      !all(required_elements %in% names(imported_data))) {
    stop("'imported_data' must be a valid list returned by import_data().",
         call. = FALSE)
  }
  
  if (is.list(mcmc_config)) {
    config <- mcmc_config
  } else {
    cfg_df           <- as.data.frame(readxl::read_excel(mcmc_config))
    cfg_df$parameter <- trimws(cfg_df$parameter)
    config           <- stats::setNames(as.list(cfg_df$value), cfg_df$parameter)
  }
  
  late_failures <- imported_data$late_failures
  additional    <- imported_data$additional
  marker_info   <- imported_data$marker_info
  data_type     <- imported_data$data_type
  
  if (data_type == "length_polymorphic") {
    results <- run_stan_sites(
      late_failures = late_failures,
      additional    = additional,
      marker_info   = marker_info,
      mcmc_config   = config,
      data_type     = data_type,
      verbose       = verbose
    )
  }  else if (data_type == "ampseq") {
    results <- run_stan_sites_ampseq(
      late_failures = imported_data$late_failures,
      additional    = imported_data$additional,
      marker_info   = imported_data$marker_info,
      mcmc_config   = config,
      verbose       = verbose
    )
  }
  
  if (!exists("results")) {
    stop("Unsupported data_type: '", data_type,
         "'. Must be 'length_polymorphic' or 'ampseq'.", call. = FALSE)
  }
  
  if (is.null(results) || length(results$ids) == 0) {
    warning("MCMC results are empty.", call. = FALSE)
    return(NULL)
  }
  
  return(results)
}


#' Summarise MCMC Classification Results
#'
#' Post-processes raw MCMC output from \code{\link{classify_infections}} into
#' interpretable summaries including posterior probabilities, convergence
#' diagnostics, and a match counting comparison table. Called automatically
#' by \code{\link{MalReBay}}.
#'
#' @param mcmc_results A list returned by \code{\link{classify_infections}}.
#' @param imported_data A list returned by \code{\link{import_data}}.
#' @param output_folder Path for saving convergence diagnostic plots.
#'   \code{NULL} skips saving.
#' @param verbose Logical. Print progress messages.
#'
#' @return A named list with \code{posterior_probabilities}, \code{comparison},
#'   \code{convergence}, and \code{mcmc_loglikelihoods}.
#'
#' @seealso \code{\link{classify_infections}}, \code{\link{save_results}},
#'   \code{\link{MalReBay}}
#'
#' @examples
#' \dontrun{
#' imported <- import_data()
#' results  <- classify_infections(imported)
#' summary  <- summarise_results(results, imported)
#' }
#'
#' @export
summarise_results <- function(mcmc_results,
                              imported_data,
                              output_folder = NULL,
                              verbose       = TRUE) {
  
  if (is.null(mcmc_results) || length(mcmc_results$ids) == 0)
    stop("'mcmc_results' is empty. Check classify_infections() ran successfully.",
         call. = FALSE)
  
  required_elements <- c("late_failures", "additional", "marker_info", "data_type")
  if (!is.list(imported_data) ||
      !all(required_elements %in% names(imported_data)))
    stop("'imported_data' must be a valid list returned by import_data().",
         call. = FALSE)
  
  late_failures <- imported_data$late_failures
  marker_info   <- imported_data$marker_info
  
  summary_list     <- list()
  convergence_list <- list()
  
  for (site in names(mcmc_results$ids)) {
    
    cls  <- mcmc_results$classifications[[site]]
    nids <- length(mcmc_results$ids[[site]])
    
    if (is.null(cls) || length(cls) == 0) {
      probs <- rep(NA_real_, nids)
    } else if (is.null(dim(cls))) {
      probs <- mean(cls, na.rm = TRUE)
    } else if (nrow(cls) == nids) {
      probs <- rowMeans(cls, na.rm = TRUE) 
    } else {
      probs <- colMeans(cls, na.rm = TRUE) 
    }
    
    summary_list[[site]] <- data.frame(
      Site        = site,
      Sample.ID   = mcmc_results$ids[[site]],
      Probability = probs,
      stringsAsFactors = FALSE
    )
    
    # Convergence checks
    stan_fit  <- mcmc_results$stan_fits[[site]]
    loglik_ch <- mcmc_results$all_chains_loglikelihood[[site]]
    save_plot <- !is.null(output_folder)
    
    diag_vals <- NULL
    if (!is.null(stan_fit) || !is.null(loglik_ch)) {
      diag_vals <- tryCatch(
        plot_likelihood_diagnostics(
          all_chains_loglikelihood = loglik_ch,
          site_name                = site,
          stan_fit                 = stan_fit,
          save_plot                = save_plot,
          output_folder            = output_folder,
          verbose                  = verbose
        ),
        error = function(e) {
          if (verbose)
            message("WARNING: Diagnostics failed for site '", site,
                    "': ", e$message)
          NULL
        }
      )
    }
    
    if (!is.null(diag_vals)) {
      convergence_list[[site]] <- data.frame(
        Site = site,
        Gelman_Rubin_Rhat     = if (!is.null(diag_vals$gelman))
          round(diag_vals$gelman$psrf[1, 1], 4) else NA_real_,
        Rank_Rhat             = if (!is.na(diag_vals$rhat_rank))
          round(diag_vals$rhat_rank, 4) else NA_real_,
        ESS_Bulk              = if (!is.na(diag_vals$ess_bulk))
          round(diag_vals$ess_bulk, 1) else NA_real_,
        ESS_Tail              = if (!is.na(diag_vals$ess_tail))
          round(diag_vals$ess_tail, 1) else NA_real_,
        Effective_Sample_Size = if (!is.null(diag_vals$ess))
          round(as.numeric(diag_vals$ess), 2) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  
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
  
  match_results    <- perform_match_counting(late_failures, marker_info)
  comparison_table <- late_failures %>%
    dplyr::mutate(
      Base.ID = trimws(gsub(" Day 0| recurrence", "", Sample.ID))
    ) %>%
    dplyr::left_join(match_results,
                     by = c("Base.ID" = "Sample.ID")) %>%
    dplyr::left_join(
      dplyr::select(posterior_probabilities,
                    Sample.ID, Probability, N_Comparable_Loci),
      by = c("Base.ID" = "Sample.ID")
    ) %>%
    dplyr::select(-Base.ID)
  
  convergence_summary <- if (length(convergence_list) > 0)
    dplyr::bind_rows(convergence_list) else NULL
  
  list(
    posterior_probabilities = posterior_probabilities,
    comparison              = comparison_table,
    convergence             = convergence_summary,
    mcmc_loglikelihoods     = mcmc_results$all_chains_loglikelihood
  )
}


#' Save Classification Results and Generate Plots
#'
#' Writes MCMC results to CSV files and generates descriptive and
#' result-based plots (Diversity, MOI, and Posterior Histograms).
#'
#' @param summary_results List from \code{summarise_results}.
#' @param imported_data   List from \code{import_data}. When \code{NULL},
#'   data-dependent plots (diversity, MOI) are skipped.
#' @param output_folder   Path to save files. \code{NULL} prints plots to
#'   the R plot window and skips CSV saving.
#' @param verbose         Logical.
#' @return A named character vector of paths to the files that were written,
#'   or \code{invisible(NULL)} when \code{output_folder} is \code{NULL}.
#'
#' @seealso \code{\link{summarise_results}}, \code{\link{MalReBay}}
#'
#' @examples
#' \dontrun{
#' imported <- import_data()
#' results  <- classify_infections(imported)
#' summary  <- summarise_results(results, imported)
#' save_results(summary, imported, output_folder = "my_results")
#' }
#'
#' @export
save_results <- function(summary_results,
                         imported_data = NULL,
                         output_folder = NULL,
                         verbose       = TRUE) {
  
  required_keys <- c("posterior_probabilities", "comparison")
  if (!is.list(summary_results) ||
      !all(required_keys %in% names(summary_results))) {
    stop("'summary_results' must be a valid list returned by summarise_results().",
         call. = FALSE)
  }
  
  if (!is.null(output_folder) && !dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
    if (verbose) message("INFO: Created output folder: ", output_folder)
  }
  
  # Generate descriptive plots — display on screen when no output_folder,
  # save to disk when output_folder is provided
  if (!is.null(imported_data)) {
    all_data <- dplyr::bind_rows(imported_data$late_failures,
                                 imported_data$additional)
    
    if (nrow(all_data) > 0) {
      p_div <- plot_markers_diversity(
        all_data,
        imported_data$data_type,
        imported_data$marker_info,
        output_folder = output_folder
      )
      if (is.null(output_folder) && !is.null(p_div)) print(p_div)
      
      p_moi <- plot_moi(all_data, output_folder = output_folder)
      if (is.null(output_folder) && !is.null(p_moi)) {
        for (p in p_moi) print(p)
      }
    }
  }
  
  plot_probability_histogram(
    summary_results,
    output_folder = output_folder,
    verbose       = verbose
  )
  
  # Stop here if no output folder — plots already shown above
  if (is.null(output_folder)) {
    if (verbose) message("INFO: No output_folder provided. Skipping file saving.")
    return(invisible(NULL))
  }
  
  # Write CSVs
  saved_paths <- character(0)
  
  pp_path <- file.path(output_folder, "posterior_probabilities.csv")
  utils::write.csv(summary_results$posterior_probabilities,
                   pp_path, row.names = FALSE)
  saved_paths["posterior_probabilities"] <- pp_path
  
  comp_path <- file.path(output_folder, "bayesian_match_counting_comparison.csv")
  utils::write.csv(summary_results$comparison,
                   comp_path, row.names = FALSE)
  saved_paths["comparison"] <- comp_path
  
  if (!is.null(summary_results$convergence)) {
    cv_path <- file.path(output_folder, "mcmc_convergence_summary.csv")
    utils::write.csv(summary_results$convergence,
                     cv_path, row.names = FALSE)
    saved_paths["convergence"] <- cv_path
  }
  
  invisible(saved_paths)
}

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
    verbose       = verbose
  )

  if (is.null(mcmc_results)) {
    if (verbose) message("classify_infections() returned no results. Returning NULL.")
    return(NULL)
  }

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
