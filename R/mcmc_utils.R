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
#' @importFrom foreach foreach %dopar% registerDoSEQ
#' @importFrom doParallel registerDoParallel
#'
#' @keywords internal
classify_infections <- function(imported_data, 
                                mcmc_config = system.file("extdata", "default_mcmc_config.xlsx", package = "MalReBay"),
                                n_workers = 1, 
                                verbose = TRUE
                                ) {
  
  # Validate imported_data structure
  required_elements <- c("late_failures", "additional", "marker_info", "data_type")
  if (!is.list(imported_data) || !all(required_elements %in% names(imported_data))) {
    stop(
      "'imported_data' must be a valid list object returned by the import_data() function.",
      call. = FALSE
    )
  }
  
  # mcmc configuration
  if (is.list(mcmc_config)) {
    config <- mcmc_config
  } else {
    cfg_df           <- as.data.frame(readxl::read_excel(mcmc_config))
    cfg_df$parameter <- trimws(cfg_df$parameter)
    config           <- setNames(as.list(cfg_df$value), cfg_df$parameter)
  }
  
  # Extract data
  late_failures <- imported_data$late_failures
  additional    <- imported_data$additional
  marker_info   <- imported_data$marker_info
  data_type     <- imported_data$data_type
  
  # Parallelization setup
  avail_cores    <- parallel::detectCores(logical = FALSE)
  workers_to_use <- min(n_workers, avail_cores)

  if (workers_to_use > 1) {
    cl <- parallel::makeCluster(workers_to_use)
    doParallel::registerDoParallel(cl)
    on.exit({ parallel::stopCluster(cl); foreach::registerDoSEQ() }, add = TRUE)
  } else {
    foreach::registerDoSEQ()
  }
  # Run MCMC
  results <- run_all_sites(
    late_failures = late_failures,
    additional    = additional,
    marker_info   = marker_info,
    mcmc_config   = config,
    output_folder = NULL,  
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
      Probability = {
        if (is.null(dim(mcmc_results$classifications[[site]]))) cls
        else rowMeans(mcmc_results$classifications[[site]], na.rm = TRUE)
      },
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

#' Save Classification Results and Generate Plots
#'
#' Writes MCMC results to CSV files and generates descriptive and 
#' result-based plots (Diversity, MOI, and Posterior Histograms).
#'
#' @param summary_results List from \code{summarise_results}.
#' @param imported_data List from \code{import_data}.
#' @param output_folder Path to save files. If \code{NULL}, plots are printed
#'   to the R plot window and CSVs are not saved.
#' @param verbose Logical.
#'
#' @keywords internal
save_results <- function(summary_results, imported_data = NULL, output_folder = NULL, verbose = TRUE) {

  # Validate summary_results
  required_elements <- c("posterior_probabilities", "comparison")
  if (!is.list(summary_results) || !all(required_elements %in% names(summary_results))) {
    stop("summary_results must be a valid list with posterior_probabilities and comparison elements.",
         call. = FALSE)
  }

  # Generate descriptive plots when imported_data is provided
  if (!is.null(imported_data)) {
    all_data <- dplyr::bind_rows(imported_data$late_failures, imported_data$additional)

    if (nrow(all_data) > 0) {
      p_div <- plot_markers_diversity(
        all_data,
        imported_data$data_type,
        imported_data$marker_info,
        output_folder = output_folder
      )
      if (!is.null(p_div)) print(p_div)

      p_moi <- plot_moi(all_data, output_folder = NULL)
      if (!is.null(p_moi) && length(p_moi) > 0) {
        print(p_moi)
        if (!is.null(output_folder)) {
          if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
          combined_moi <- ggpubr::ggarrange(plotlist = p_moi, ncol = 1)
          ggplot2::ggsave(
            file.path(output_folder, "moi_per_marker_by_site.png"),
            combined_moi,
            width  = 12,
            height = 6 * length(p_moi)
          )
        }
      }
    }
  }

  # Generate posterior probabilities histogram
  p_hist <- plot_probability_histogram(summary_results, output_folder = output_folder, verbose = verbose)
  if (!is.null(p_hist)) print(p_hist)

  # Save CSV Data
  if (is.null(output_folder)) {
    if (verbose) message("INFO: No output folder provided. CSV files not saved.")
    return(invisible(NULL))
  }

  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

  # Save CSVs
  pp_path <- file.path(output_folder, "posterior_probabilities.csv")
  ct_path <- file.path(output_folder, "bayesian_match_counting_comparison.csv")

  utils::write.csv(summary_results$posterior_probabilities, pp_path, row.names = FALSE)
  utils::write.csv(summary_results$comparison,              ct_path, row.names = FALSE)

  paths <- c(posterior_probabilities = pp_path, comparison = ct_path)

  if (!is.null(summary_results$convergence)) {
    cv_path <- file.path(output_folder, "mcmc_convergence_summary.csv")
    utils::write.csv(summary_results$convergence, cv_path, row.names = FALSE)
    paths <- c(paths, convergence = cv_path)
  }

  invisible(paths)
}