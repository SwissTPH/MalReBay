#' Classify malaria infections using Bayesian MCMC
#'
#' @param imported_data List returned by \code{import_data()}.
#' @param mcmc_config   Path to MCMC configuration Excel file.
#' @param verbose       Print progress messages.
#' @return List of MCMC results passed to \code{summarise_results()}.
#' @keywords internal
classify_infections <- function(imported_data,
                                mcmc_config = system.file(
                                  "extdata", "default_mcmc_config.xlsx",
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
#' @keywords internal
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
#' @keywords internal
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

  # If no folder is provided, we skip all saving logic
  if (is.null(output_folder)) {
    if (verbose) message("INFO: No output_folder provided. Skipping file saving.")
    return(invisible(NULL))
  }

  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
    if (verbose) message("INFO: Created output folder: ", output_folder)
  }

  if (!is.null(imported_data)) {
    all_data <- dplyr::bind_rows(imported_data$late_failures,
                                 imported_data$additional)

    if (nrow(all_data) > 0) {
      plot_markers_diversity(
        all_data,
        imported_data$data_type,
        imported_data$marker_info,
        output_folder = output_folder
      )

      plot_moi(all_data, output_folder = output_folder)
    }
  }

  plot_probability_histogram(
    summary_results,
    output_folder = output_folder,
    verbose       = verbose
  )

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