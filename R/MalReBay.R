#' Classify Infections Using a Bayesian MCMC Framework
#'
#' @description
#' This is the main user-facing function of the `MalReBay` package. It runs the
#' complete Bayesian analysis workflow to classify parasite infections as either
#' recrudescence or reinfection based on genotyping data. The function handles
#' data import, MCMC simulation with automatic convergence checking, and the
#' summarization and saving of final results.
#'
#' @details
#' The function operates by first importing and validating the data from the
#' provided Excel file. It automatically detects the data type (length-polymorphic
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
#'     convergence is assessed. (Default: 10000).
#'   \item \strong{`max_iterations`}: The maximum total number of iterations before
#'     the simulation is forcibly stopped, even if not converged. (Default: 100000).
#'   \item \strong{`burn_in_frac`}: The fraction of initial samples to discard
#'     from each chain before summarizing results. (Default: 0.25).
#'   \item \strong{`rhat_threshold`}: The Gelman-Rubin diagnostic threshold for
#'     convergence. (Default: 1.01).
#'   \item \strong{`ess_threshold`}: The Effective Sample Size threshold for
#'     convergence. (Default: 400).
#'   \item \strong{`record_hidden_alleles`}: A logical flag. If `TRUE`, the full
#'     state of imputed hidden alleles is saved. This can generate very large
#'     output files and is mainly for debugging. (Default: FALSE).
#' }
#' Any parameters not specified in the list will use the default values.
#'
#' @param input_filepath A character string. The full path to the input Excel file
#'   containing the cleaned genotyping data.
#' @param mcmc_config A list of MCMC configuration parameters. See Details for
#'   a full list of options and their defaults.
#' @param output_folder A character string specifying the path to the folder where
#'   all results (summary CSVs, diagnostic plots) will be saved. The folder
#'   will be created if it does not exist.
#' @param verbose A logical. If `TRUE` (the default), the function will print
#'   progress messages, warnings about missing markers, and other information
#'   to the console. If `FALSE`, it will run silently, only showing critical
#'   errors.
#'
#' @return
#' A list containing two data frames:
#' \item{summary}{A data frame summarizing the main results for each patient,
#'   including the posterior probability of recrudescence and a summary of the
#'   genetic data available for comparison.}
#' \item{marker_details}{A detailed data frame providing the mean likelihood ratio
#'   and mean allele distance for each genetic marker for each patient, offering
#'   locus-specific insights into the classification.}
#' Two CSV files corresponding to these data frames are also saved to the
#' `output_folder`.
#'
#' @importFrom utils modifyList write.csv
#' @importFrom dplyr bind_rows rename left_join full_join
#'
#' @export
#'
#' @examples
#' # Get the path to the example data file included with the package
#' example_file <- system.file("extdata", "Angola_2021_TES_7NMS.xlsx",
#'                             package = "MalReBay")
#'
#' # Create a temporary directory to save the results for this example
#' temp_output_dir <- file.path(tempdir(), "MalReBay_example_silent")
#' dir.create(temp_output_dir)
#'
#' # Define a minimal MCMC configuration for a quick and illustrative run
#' # NOTE: These settings are for demonstration only and are not sufficient
#' # for a real analysis.
#' quick_mcmc_config <- list(
#'   n_chains = 2,
#'   chunk_size = 500,
#'   max_iterations = 500, # Set max_iterations equal to chunk_size for one run
#'   rhat_threshold = 1.1, # Relaxed for the example
#'   ess_threshold = 50   # Relaxed for the example
#' )
#'
#' \dontrun{
#' # Run the analysis silently
#' results_list <- classify_infections(
#'   input_filepath = example_file,
#'   mcmc_config = quick_mcmc_config,
#'   output_folder = temp_output_dir,
#'   verbose = FALSE # Set to FALSE to suppress messages
#' )
#'
#' # Run the analysis with full messages (the default behavior)
#' # results_list_verbose <- classify_infections(
#' #   input_filepath = example_file,
#' #   mcmc_config = quick_mcmc_config,
#' #   output_folder = temp_output_dir,
#' #   verbose = TRUE
#' # )
#'
#' # View the top rows of the main summary table
#' print(head(results_list$summary))
#' }
#'
classify_infections <- function(input_filepath,
                                mcmc_config,
                                output_folder,
                                verbose = TRUE) {
  
  # Data Import and Validation ---
  imported_data <- import_data(filepath = input_filepath)
  genotypedata_latefailures <- imported_data$late_failures
  additional_genotypedata <- imported_data$additional
  marker_information <- imported_data$marker_info
  data_type <- imported_data$data_type
  
  # Clean column names
  colnames(genotypedata_latefailures) <- trimws(colnames(genotypedata_latefailures))
  if (nrow(additional_genotypedata) > 0) {
    colnames(additional_genotypedata) <- trimws(colnames(additional_genotypedata))
  }
  colnames(genotypedata_latefailures) <- gsub("R033_", "RO33_", colnames(genotypedata_latefailures))
  colnames(additional_genotypedata) <- gsub("R033_", "RO33_", colnames(additional_genotypedata))
  
  defaults <- list(n_chains = 4,
                   rhat_threshold = 1.01,
                   ess_threshold = 400,
                   chunk_size = 10000,
                   max_iterations = 100000,
                   burn_in_frac = 0.25,
                   record_hidden_alleles = FALSE)
  config <- utils::modifyList(defaults, mcmc_config)
  
  if (data_type == "length_polymorphic") {
    data_marker_columns <- grep("_\\d+$", colnames(genotypedata_latefailures), value = TRUE)
    available_base_markers <- unique(gsub("_\\d+$", "", data_marker_columns))
    requested_base_markers <- marker_information$marker_id
    markers_to_use <- intersect(requested_base_markers, available_base_markers)
    
    if (length(markers_to_use) == 0) {
      stop("None of the markers specified in `marker_information` were found in the input data.")
    }
    missing_markers <- setdiff(requested_base_markers, available_base_markers)
    if (length(missing_markers) > 0 && verbose) {
      warning("The following requested markers were NOT found and will be ignored: ",
              paste(missing_markers, collapse = ", "), call. = FALSE)
    }
    
    extra_markers_in_data <- setdiff(available_base_markers, requested_base_markers)
    if (length(extra_markers_in_data) > 0 && verbose) {
      message("INFO: The following markers were found in your data but not requested: ",
              paste(extra_markers_in_data, collapse = ", "))
    }
    if (verbose) message("Proceeding with analysis for these markers: ", paste(markers_to_use, collapse = ", "))
    
  } else if (data_type == "ampseq") {
    data_marker_columns <- grep("_allele_\\d+$", colnames(genotypedata_latefailures), value = TRUE)
    available_base_markers <- unique(gsub("_allele_\\d+$", "", data_marker_columns))
    requested_base_markers <- marker_information$marker_id
    markers_to_use <- intersect(requested_base_markers, available_base_markers)
    
    if (length(markers_to_use) == 0) {
      stop("None of the markers from the ampseq data were found in the marker_information sheet.")
    }
    if (verbose) {
      message("Amplicon sequencing data detected. The following haplotype markers will be used for analysis: ",
              paste(markers_to_use, collapse = ", "))
    }
  } else {
    stop("Unknown data type detected: ", data_type)
  }
  
  final_marker_info <- marker_information[marker_information$marker_id %in% markers_to_use, ]
  all_marker_base_names <- if (data_type == "ampseq") {
    gsub("_allele_\\d+$", "", data_marker_columns)
  } else {
    gsub("_\\d+$", "", data_marker_columns)
  }
  final_marker_cols <- data_marker_columns[all_marker_base_names %in% markers_to_use]
  
  id_cols <- setdiff(colnames(genotypedata_latefailures), grep("_allele_\\d+|_\\d+$", colnames(genotypedata_latefailures), value = TRUE))
  latefailures_subset <- genotypedata_latefailures[, c(id_cols, final_marker_cols)]
  additional_subset <- if (nrow(additional_genotypedata) > 0) {
    additional_genotypedata[, c(id_cols, final_marker_cols)]
  } else {
    additional_genotypedata
  }
  
  if (verbose) message("The model is running MCMC analysis...")
  results <- run_all_sites(
    genotypedata_latefailures = latefailures_subset,
    additional_genotypedata = additional_subset,
    marker_info_subset = final_marker_info,
    mcmc_config = config,
    data_type = data_type,
    output_folder = output_folder,
    verbose = verbose
  )
  
  if (is.null(results) || length(results$ids) == 0) {
    warning("MCMC results are empty. No summary table can be generated.", call. = FALSE)
    return(
      list(
        summary = data.frame(Site=character(), Sample.ID=character(), Probability=numeric(), N_Available_D0=integer(), N_Available_DF=integer(), N_Comparable_Loci=integer()),
        marker_details = NULL
      )
    )
  }
  
  summary_list_probs <- list()
  for (site_name in names(results$ids)) {
    site_ids <- results$ids[[site_name]]
    site_probs_matrix <- results$classifications[[site_name]]
    
    if (is.null(site_ids) || is.null(site_probs_matrix) || length(site_ids) == 0 || nrow(site_probs_matrix) == 0) {
      next
    }
    
    probabilities <- rowMeans(site_probs_matrix, na.rm = TRUE)
    summary_list_probs[[site_name]] <- data.frame(
      Site = site_name,
      Sample.ID = site_ids,
      Probability = probabilities,
      stringsAsFactors = FALSE
    )
  }
  
  if (length(summary_list_probs) == 0) {
    warning("No probability data was successfully summarized across all sites.", call. = FALSE)
    return(
      list(
        summary = data.frame(Site=character(), Sample.ID=character(), Probability=numeric(), N_Available_D0=integer(), N_Available_DF=integer(), N_Comparable_Loci=integer()),
        marker_details = NULL
      )
    )
  }
  
  probabilities_df <- dplyr::bind_rows(summary_list_probs)
  locus_summary_df <- dplyr::bind_rows(results$locus_summary, .id = "Site")
  locus_summary_df <- dplyr::rename(locus_summary_df,
                                    Sample.ID = "patient_id",
                                    N_Available_D0 = "n_available_d0",
                                    N_Available_DF = "n_available_df",
                                    N_Comparable_Loci = "n_comparable_loci")
  
  summary_df <- dplyr::left_join(probabilities_df, locus_summary_df, by = c("Sample.ID", "Site"))
  
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  csv_path <- file.path(output_folder, "probability_of_recrudescence_summary.csv")
  utils::write.csv(summary_df, csv_path, row.names = FALSE)
  if (verbose) message("Summary table saved to: ", csv_path)
  
  marker_details_list <- list()
  for (site_name in names(results$ids)) {
    if (is.null(results$locus_lrs[[site_name]])) next
    
    site_ids <- results$ids[[site_name]]
    site_lrs <- results$locus_lrs[[site_name]]
    site_dists <- results$locus_dists[[site_name]]
    marker_names <- results$locinames[[site_name]] %||% dimnames(site_lrs)[[2]] %||% paste0("Locus_", seq_len(dim(site_lrs)[2]))
    
    mean_lrs <- apply(site_lrs, c(1, 2), mean, na.rm = TRUE)
    mean_dists <- apply(site_dists, c(1, 2), mean, na.rm = TRUE)
    
    if (length(site_ids) == nrow(mean_lrs) && length(marker_names) == ncol(mean_lrs)) {
      dimnames(mean_lrs) <- list(site_ids, marker_names)
    }
    if (length(site_ids) == nrow(mean_dists) && length(marker_names) == ncol(mean_dists)) {
      dimnames(mean_dists) <- list(site_ids, marker_names)
    }
    
    lrs_df <- as.data.frame(as.table(mean_lrs), stringsAsFactors = FALSE)
    colnames(lrs_df) <- c("Sample.ID", "Marker", "Mean_Likelihood_Ratio")
    
    dists_df <- as.data.frame(as.table(mean_dists), stringsAsFactors = FALSE)
    colnames(dists_df) <- c("Sample.ID", "Marker", "Mean_Distance")
    
    site_marker_details <- dplyr::full_join(lrs_df, dists_df, by = c("Sample.ID", "Marker"))
    site_marker_details$Site <- site_name
    marker_details_list[[site_name]] <- site_marker_details
  }
  
  marker_details_df <- NULL
  if (length(marker_details_list) > 0) {
    marker_details_df <- dplyr::bind_rows(marker_details_list)
    marker_details_df$Interpretation <- ifelse(marker_details_df$Mean_Likelihood_Ratio > 1, "Supports Recrudescence", "Supports Reinfection/Error")
    marker_csv_path <- file.path(output_folder, "marker_level_summary.csv")
    utils::write.csv(marker_details_df, marker_csv_path, row.names = FALSE)
    if (verbose) message("Detailed marker-level summary saved to: ", marker_csv_path)
  }
  
  return(
    list(
      summary = summary_df,
      marker_details = marker_details_df
    )
  )
}

# Helper function for NULL coalescing
`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}