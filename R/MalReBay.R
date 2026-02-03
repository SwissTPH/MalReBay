#' Classify Infections Using a Bayesian MCMC Framework
#'
#' @description
#' This is the main function of the `MalReBay` package. It runs the
#' complete Bayesian analysis workflow to classify parasite infections as either
#' recrudescence or reinfection based on genotyping data. The function takes a processed data 
#' object, runs MCMC simulation with automatic convergence 
#' checking, and summarizes the final results.
#'
#' @details
#' The function operates by validating the data from the
#' provided input object. It automatically detects the data type (length-polymorphic
#' or amplicon sequencing) and selects the appropriate MCMC engine. The MCMC simulation
#' is run in parallel across multiple chains and proceeds in chunks, stopping
#' automatically when convergence criteria are met or the maximum number of
#' iterations is reached.
#'
#' The `mcmc_config` list provides fine-grained control over the simulation.
#' Key parameters include:
#' \itemize{
#'   \item \strong{`n_chains`}: Number of parallel chains to run. (Default: 4).
#'   \item \strong{`chunk_size`}: Number of iterations per chunk. After each chunk,
#'     convergence is assessed. (Default: 1000).
#'   \item \strong{`max_iterations`}: The maximum total number of iterations before
#'     the simulation is forcibly stopped. (Default: 10000).
#'   \item \strong{`burn_in_frac`}: The fraction of initial samples to discard
#'     from each chain before summarizing results. (Default: 0.25).
#'   \item \strong{`record_hidden_alleles`}: A logical flag. If `TRUE`, the full
#'     state of imputed hidden alleles is saved. (Default: FALSE).
#' }
#'
#' @param imported_data A list object (typically from an import helper) containing:
#'   \itemize{
#'     \item \code{late_failures}: Data frame of paired Day 0/Day Failure samples.
#'     \item \code{additional}: Data frame of background/population samples.
#'     \item \code{marker_info}: Data frame with marker names and repeat lengths.
#'     \item \code{data_type}: Character, either "length_polymorphic" or "ampseq".
#'   }
#' @param mcmc_config A list of MCMC configuration parameters. See Details.
#' @param output_folder A character string specifying the path to the folder where
#'   all results (summary CSVs, diagnostic plots) will be saved. The folder
#'   will be created if it does not exist. If a \code{convergence_diagnosis}
#'   sub-folder exists from a previous run, it will be automatically updated.
#' @param n_workers An integer specifying the number of parallel cores to use.
#'   Defaults to 1. If \code{n_workers > 1}, the analysis is parallelized by site.
#' @param verbose A logical. If \code{TRUE}, the function prints progress messages.
#'
#' @return A list containing three elements:
#' \item{summary}{A data frame summarizing the classification probability for each patient.}
#' \item{marker_details}{A detailed data frame providing the mean likelihood ratio and distance per marker.}
#' \item{comparison}{A comparison table matching Bayesian results with basic match-counting logic.}
#' 
#' @importFrom utils modifyList write.csv
#' @importFrom dplyr bind_rows rename left_join select mutate
#' @importFrom future plan multisession sequential availableCores
#' @importFrom ggplot2 ggplot aes geom_jitter geom_violin geom_text scale_fill_brewer facet_grid labs theme_classic ggsave
#' @importFrom tidyr pivot_longer crossing replace_na
#' @importFrom grDevices png dev.off
#' @importFrom graphics hist
#'
#' @export
#' 
#' @examples
#' \dontrun{
#' # 1. Import and process the data
#' processed_data <- import_data(filepath = "my_data.xlsx")
#'
#' # 2. Run the analysis
#' results <- classify_infections(
#'   imported_data = processed_data,
#'   n_workers = 2
#' )
#' }
classify_infections <- function(imported_data, mcmc_config = list(), output_folder = "results", n_workers = 1, verbose = TRUE){
  
  # Validation and Directory Setup
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  
  required_elements <- c("late_failures", "additional", "marker_info", "data_type")
  if (!is.list(imported_data) || !all(required_elements %in% names(imported_data))) {
    stop(
      "'imported_data' must be a valid list object returned by the import_data() function.",
      call. = FALSE
    )
  }
  
  convergence_dir <- file.path(output_folder, "convergence_diagnosis")
  if (dir.exists(convergence_dir)) {
    if (verbose) message("INFO: Removing existing convergence diagnosis folder to ensure fresh results.")
    unlink(convergence_dir, recursive = TRUE)
  }
  
  config <- utils::modifyList(list(
    n_chains = 4, chunk_size = 1000, max_iterations = 10000, 
    burn_in_frac = 0.25, record_hidden_alleles = FALSE
  ), mcmc_config)
  
  # Extract original data
  genotypedata_latefailures <- imported_data$late_failures
  additional_genotypedata   <- imported_data$additional
  marker_information        <- imported_data$marker_info
  data_type                 <- imported_data$data_type
  
  # Marker identification
  marker_suffix_regex <- "(_allele_|_)\\d+$"
  allele_cols <- colnames(genotypedata_latefailures)[3:ncol(genotypedata_latefailures)]
  
  # Define the variables exactly as you have them
  base_names_in_data <- gsub(marker_suffix_regex, "", allele_cols)
  markers_to_use     <- intersect(marker_information$marker_id, unique(base_names_in_data))
  
  if (length(markers_to_use) == 0) stop("No matching markers found between metadata and data.")
  
  # Select the valid allele columns
  id_cols <- allele_cols[base_names_in_data %in% markers_to_use]
  
  # Create subsets
  latefailures_subset <- genotypedata_latefailures[, c(colnames(genotypedata_latefailures)[1:2], id_cols)]
  additional_subset  <- if (nrow(additional_genotypedata) > 0) {
    additional_genotypedata[, c(colnames(additional_genotypedata)[1:2], id_cols)]
  } else {
    additional_genotypedata
  }
  
  final_marker_info <- marker_information[marker_information$marker_id %in% markers_to_use, ]
  
  if (verbose) message("INFO: Using markers: ", paste(markers_to_use, collapse = ", "))
  
  # MCMC parallelization setup
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  
  # Determine worker count
  avail_cores <- future::availableCores()
  workers_to_use <- min(n_workers, avail_cores)
  
  if (workers_to_use > 1) {
    if (verbose) message("INFO: Running in parallel on ", workers_to_use, " cores.")
    future::plan(future::multisession, workers = workers_to_use)
  } else {
    future::plan(future::sequential)
  }
  
  # MCMC engine
  results <- run_all_sites(
    genotypedata_latefailures = latefailures_subset,
    additional_genotypedata   = additional_subset,
    marker_info_subset        = final_marker_info,
    mcmc_config               = config,
    data_type                 = data_type,
    output_folder             = output_folder,
    verbose                   = verbose
  )
  
  if (is.null(results) || length(results$ids) == 0) {
    warning("MCMC results are empty.")
    return(NULL)
  }
  
  # Processing output files
  summary_list <- list()
  details_list <- list()
  convergence_list <- list() # This collects the R-hat/ESS for the final table
  
  for (site in names(results$ids)) {
    # Calculate Probabilities
    summary_list[[site]] <- data.frame(
      Site = site, Sample.ID = results$ids[[site]],
      Probability = rowMeans(results$classifications[[site]], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    
    # Capture Diagnostics
    diag_vals <- NULL 
    if (!is.null(results$all_chains_loglikelihood[[site]])) {
      diag_vals <- plot_likelihood_diagnostics(
        results$all_chains_loglikelihood[[site]], 
        site, 
        output_folder = output_folder,
        verbose = verbose
      )
    }
    
    # Store R-hat and ESS into the convergence_list
    if (!is.null(diag_vals)) {
      convergence_list[[site]] <- data.frame(
        Site = site,
        Gelman_Rubin_Rhat = if(!is.null(diag_vals$gelman)) round(diag_vals$gelman$psrf[1,1], 4) else NA,
        Effective_Sample_Size = if(!is.null(diag_vals$ess)) round(as.numeric(diag_vals$ess), 2) else NA,
        stringsAsFactors = FALSE
      )
    }
    
    # Marker Details
    lrs  <- apply(results$locus_lrs[[site]], c(1, 2), mean, na.rm = TRUE)
    dsts <- apply(results$locus_dists[[site]], c(1, 2), mean, na.rm = TRUE)
    site_details <- as.data.frame(as.table(lrs))
    colnames(site_details) <- c("Sample.ID", "Marker", "Mean_LR")
    site_details$Mean_Distance <- as.vector(dsts)
    site_details$Site <- site
    details_list[[site]] <- site_details
  }
  
  final_summary <- dplyr::bind_rows(summary_list) %>%
    dplyr::left_join(dplyr::bind_rows(results$locus_summary, .id = "Site") %>%
                       dplyr::rename(Sample.ID=patient_id, N_Available_D0=n_available_d0, 
                                     N_Available_DF=n_available_df, N_Comparable_Loci=n_comparable_loci),
                     by = c("Sample.ID", "Site"))
  
  marker_details_df <- dplyr::bind_rows(details_list) %>%
    dplyr::mutate(Interpretation = ifelse(Mean_LR > 1, "Recrudescence", "Reinfection/Error"))
  
  # Plots & Diversity
  all_genotypedata_for_plots <- dplyr::bind_rows(latefailures_subset, additional_subset)
  if (nrow(all_genotypedata_for_plots) > 0) {
    plot_markers_diversity(all_genotypedata_for_plots, data_type, final_marker_info, output_folder = output_folder)
    plot_moi(all_genotypedata_for_plots, output_folder = output_folder)
  }
  
  # Posteriror Probability Histogram
  if (nrow(final_summary) > 0) {    
    grDevices::png(file.path(output_folder, "recrudescence_probability_histogram.png"), 
                   width = 8, height = 6, units = "in", res = 300)
    graphics::hist(as.numeric(as.character(final_summary$Probability)), 
                   breaks = seq(0, 1, by = 0.05), # Force bins from 0 to 1
                   col = "skyblue",
                   main = "Posterior Probability Distribution",
                   xlab = "Probability of Recrudescence", ylab = "Number of Patients")
    grDevices::dev.off()
    if(verbose) message("INFO: Probability histogram saved to ", output_folder)
  }
  
  
  # Match counting comparison
  match_results <- perform_match_counting(latefailures_subset, final_marker_info)
  final_comparison_table <- genotypedata_latefailures %>%
    dplyr::mutate(Patient.ID = gsub(" (Day 0|Day Failure)$", "", Sample.ID)) %>%
    dplyr::left_join(match_results, by = c("Patient.ID" = "Sample.ID")) %>%
    dplyr::left_join(dplyr::select(final_summary, Sample.ID, Probability, N_Comparable_Loci), 
                     by = c("Patient.ID" = "Sample.ID"))
  
  # Export the Global MCMC Convergence Summary
  if (length(convergence_list) > 0) {
    mcmc_summary_df <- dplyr::bind_rows(convergence_list)
    utils::write.csv(mcmc_summary_df, file.path(output_folder, "mcmc_convergence_summary.csv"), row.names = FALSE)
    if (verbose) message("INFO: Global MCMC convergence summary saved to: ", output_folder)
  }
  
  # Write Files
  utils::write.csv(final_summary, file.path(output_folder, "posterior_probabilities.csv"), row.names = FALSE)
  utils::write.csv(marker_details_df, file.path(output_folder, "markers_classification_analysis.csv"), row.names = FALSE)
  utils::write.csv(final_comparison_table, file.path(output_folder, "bayesian_match_counting_comparison.csv"), row.names = FALSE)
  
  return(list(
    summary = final_summary, 
    marker_details = marker_details_df, 
    comparison = final_comparison_table,
    mcmc_loglikelihoods = results$all_chains_loglikelihood 
  ))
}