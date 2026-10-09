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
#' @param suppress_warnings Logical. If \code{TRUE} (default), silences the
#'   MCMC sampler's own console warnings about divergent transitions,
#'   treedepth, and E-BFMI (the "N of M transitions ended with a divergence"
#'   messages cmdstanr prints as soon as sampling finishes, independently of
#'   \code{verbose}). Set to \code{FALSE} to see them -- e.g. while tuning
#'   \code{adapt_delta} in \code{mcmc_config} for a problematic site. This
#'   only affects what gets printed; the same diagnostics are always
#'   available afterwards via \code{\link{summarise_results}}'s
#'   \code{convergence} table.
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
                                                  "default_mcmc_config.rds",
                                                  package = "MalReBay"),
                                n_workers         = 1,
                                verbose           = TRUE,
                                suppress_warnings = TRUE) {

  required_elements <- c("late_failures", "additional", "marker_info", "data_type")
  if (!is.list(imported_data) ||
      !all(required_elements %in% names(imported_data))) {
    stop("'imported_data' must be a valid list returned by import_data().",
         call. = FALSE)
  }

  if (is.list(mcmc_config)) {
    config <- mcmc_config
  } else {
    cfg_df           <- as.data.frame(readRDS(mcmc_config))
    cfg_df$parameter <- trimws(cfg_df$parameter)
    config           <- stats::setNames(as.list(cfg_df$value), cfg_df$parameter)
  }

  late_failures <- imported_data$late_failures
  additional    <- imported_data$additional
  marker_info   <- imported_data$marker_info

  results <- run_stan_sites(
    late_failures     = late_failures,
    additional        = additional,
    marker_info       = marker_info,
    mcmc_config       = config,
    verbose           = verbose,
    suppress_warnings = suppress_warnings
  )

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
#' @param prob_threshold Numeric in \verb{[0, 1]}. The posterior probability
#'   at or above which a recurrence is classified as \code{"Recrudescence"}
#'   (below it, \code{"New infection"}); used to add the \code{Classification}
#'   column to \code{posterior_probabilities} (and \code{comparison}).
#'   Defaults to \code{0.5}, the natural cutoff under this model's equal
#'   50/50 prior (see @sec-interpretation-probability in the analysis
#'   notebook).
#' @param plots Logical. If \code{FALSE}, convergence diagnostic plots are
#'   not saved to \code{output_folder}; the diagnostics themselves are still
#'   computed and returned in \code{convergence}. Defaults to \code{TRUE}.
#'
#' @return A named list with \code{posterior_probabilities}, \code{comparison},
#'   \code{convergence}, \code{mcmc_loglikelihoods}, and the
#'   \code{prob_threshold} used. \code{posterior_probabilities} includes a
#'   \code{Classification} column derived from \code{prob_threshold}.
#'   \code{comparison} has one row per patient: the per-marker match-counting
#'   calls, the two overall match-counting calls (e.g.
#'   \code{Match_counting_2of3}/\code{Match_counting_3of3} for MSP panels,
#'   \code{Match_counting_5of7}/\code{Match_counting_7of7} for a 7-marker
#'   panel), \code{MalReBay_probability} and \code{MalReBay_classification}.
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
                              output_folder  = NULL,
                              verbose        = TRUE,
                              prob_threshold = 0.5,
                              plots          = TRUE) {

  if (!is.numeric(prob_threshold) || length(prob_threshold) != 1 ||
      is.na(prob_threshold) || prob_threshold < 0 || prob_threshold > 1) {
    stop("'prob_threshold' must be a single number between 0 and 1.", call. = FALSE)
  }
  
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
    save_plot <- plots && !is.null(output_folder)
    
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
          Sample.ID            = patient_id,
          N_Markers_Day0       = n_available_d0,
          N_Markers_Recurrence = n_available_df,
          N_Markers_Compared   = n_comparable_loci
        ),
      by = c("Sample.ID", "Site")
    ) %>%
    dplyr::mutate(
      Classification = ifelse(Probability >= prob_threshold, "Recrudescence", "New infection")
    )

  # One row per patient: the per-marker match-counting calls (R/NI/IND/ERR),
  # the overall match-counting calls (2/3 & 3/3 for MSP panels, else ~70% &
  # 100% of markers, e.g. 5/7 & 7/7), then MalReBay's probability and
  # classification -- no raw alleles
  match_results    <- perform_match_counting(late_failures, marker_info)
  marker_result_cols <- setdiff(colnames(match_results),
                                c("Sample.ID", "Number_Matches", "Number_Loci_Compared"))
  patient_sites    <- late_failures %>%
    dplyr::transmute(Sample.ID = trimws(gsub(" Day 0| recurrence", "", Sample.ID)), Site) %>%
    dplyr::distinct()
  comparison_table <- patient_sites %>%
    dplyr::inner_join(match_results, by = "Sample.ID") %>%
    dplyr::left_join(
      dplyr::select(posterior_probabilities, Site, Sample.ID, Probability, Classification),
      by = c("Site", "Sample.ID")
    ) %>%
    dplyr::select(Site, Sample.ID, dplyr::all_of(marker_result_cols), Probability, Classification)

  # WHO_loose/WHO_strict are 1/0/NA (NA = too few markers scored to apply the rule)
  who_result   <- build_who_table(comparison_table, marker_info)
  call_label   <- function(x) dplyr::case_when(x == 1 ~ "Recrudescence", x == 0 ~ "New infection",
                                                TRUE ~ NA_character_)
  # "Match counting 5/7" -> "Match_counting_5of7"
  col_name     <- function(label) gsub("/", "of", gsub(" ", "_", label))
  loose_col    <- col_name(who_result$who_loose_label)
  strict_col   <- col_name(who_result$who_strict_label)
  # MSP panels also get the combined msp1/msp2 family calls the 2/3 rule is applied to
  family_cols  <- intersect(c("msp1", "msp2"), setdiff(colnames(who_result$table), marker_result_cols))
  comparison_table <- who_result$table
  comparison_table[[loose_col]]  <- call_label(comparison_table$WHO_loose)
  comparison_table[[strict_col]] <- call_label(comparison_table$WHO_strict)
  comparison_table <- comparison_table %>%
    dplyr::select(Site, Sample.ID, dplyr::all_of(c(marker_result_cols, family_cols, loose_col, strict_col)),
                  MalReBay_probability = Probability, MalReBay_classification = Classification)
  
  convergence_summary <- if (length(convergence_list) > 0)
    dplyr::bind_rows(convergence_list) else NULL
  
  list(
    posterior_probabilities = posterior_probabilities,
    comparison              = comparison_table,
    convergence             = convergence_summary,
    mcmc_loglikelihoods     = mcmc_results$all_chains_loglikelihood,
    prob_threshold          = prob_threshold
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
#' @param plots           Logical. If \code{FALSE}, no plots are generated,
#'   shown or saved; only the CSV files are written. Defaults to \code{TRUE}.
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
                         verbose       = TRUE,
                         plots         = TRUE) {
  
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
  if (plots && !is.null(imported_data)) {
    all_data <- dplyr::bind_rows(imported_data$late_failures,
                                 imported_data$additional)
    
    if (nrow(all_data) > 0) {
      p_div <- plot_markers_diversity(
        all_data,
        imported_data$data_type,
        imported_data$marker_info,
        output_folder = output_folder
      )
      if (is.null(output_folder) && !is.null(p_div)) {
        print(p_div$all_sites)
        for (p in p_div$by_site) print(p)
      }
      
      p_moi <- plot_moi(all_data, output_folder = output_folder)
      if (is.null(output_folder) && !is.null(p_moi)) {
        for (p in p_moi) print(p)
      }
    }
    
    plot_comparison_heatmap(
      summary_results = summary_results,
      marker_info     = imported_data$marker_info,
      output_folder   = output_folder,
      verbose         = verbose
    )
    
  }
  
  if (plots) {
    plot_probability_histogram(
      summary_results,
      output_folder = output_folder,
      verbose       = verbose
    )
  }
  
  
  # Stop here if no output folder — plots (if any) already shown above
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
#' @param suppress_warnings Logical. If \code{TRUE} (default), silences the
#'   MCMC sampler's own console warnings about divergent transitions,
#'   treedepth, and E-BFMI. See \code{\link{classify_infections}} for
#'   details -- these diagnostics remain available afterwards via
#'   \code{summary_results$convergence} regardless of this setting.
#' @param prob_threshold  Numeric in \verb{[0, 1]}. The posterior probability
#'   at or above which a recurrence is classified as \code{"Recrudescence"}
#'   (below it, \code{"New infection"}). Defaults to \code{0.5}, the natural
#'   cutoff under this model's equal 50/50 prior. See
#'   \code{\link{summarise_results}} for where this is applied.
#' @param plots           Logical. If \code{FALSE}, skips all plots
#'   (descriptive, result and convergence diagnostic plots); convergence
#'   diagnostics are still computed and printed. Defaults to \code{TRUE}.
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
    filepath        = system.file("extdata", "Dataset_microsatellite_panel.xlsx",
                                  package = "MalReBay"),
    marker_filepath = system.file("extdata", "makers_details.xlsx",
                                  package = "MalReBay"),
    additional_filepath = NULL,
    mcmc_config     = system.file("extdata", "default_mcmc_config.rds",
                                  package = "MalReBay"),
    output_folder   = NULL,
    n_workers       = 1,
    verbose         = TRUE,
    suppress_warnings = TRUE,
    prob_threshold  = 0.5,
    plots           = TRUE
) {

  if (verbose) message("Starting MalReBay pipeline...")

  imported_data <- import_data(
    filepath             = filepath,
    additional_filepath  = additional_filepath,
    marker_filepath      = marker_filepath,
    verbose              = verbose
  )

  mcmc_results <- classify_infections(
    imported_data     = imported_data,
    mcmc_config       = mcmc_config,
    verbose           = verbose,
    suppress_warnings = suppress_warnings
  )

  if (is.null(mcmc_results)) {
    if (verbose) message("classify_infections() returned no results. Returning NULL.")
    return(NULL)
  }

  summary_results <- summarise_results(
    mcmc_results   = mcmc_results,
    imported_data  = imported_data,
    output_folder  = output_folder,
    verbose        = verbose,
    prob_threshold = prob_threshold,
    plots          = plots
  )

  save_results(
    summary_results = summary_results,
    imported_data   = imported_data,
    output_folder   = output_folder,
    verbose         = verbose,
    plots           = plots
  )

  invisible(summary_results)
}
