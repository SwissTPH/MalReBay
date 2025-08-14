#' @title Classify infections using a Bayesian MCMC Framework
#' @description
#'   This is the main user-facing function to run the Bayesian analysis. It
#'   manages the entire workflow, including data import via its internal helper
#'   `import_data`, MCMC simulation, convergence assessment, and summarization
#'   of results.
#'
#' @param input_filepath A string containing the path to the input Excel file.
#' @param mcmc_config A list containing MCMC configuration parameters. See Details.
#' @param output_folder The folder where the results will be saved.
#'
#' @details
#' The `mcmc_config` list can contain the following parameters:
#' \itemize{
#'   \item `n_chains`: Number of parallel chains (Default: 4).
#'   \item `chunk_size`: Iterations per convergence check (Default: 2000).
#'   \item `max_iterations`: Max iterations before stopping (Default: 20000).
#'   \item `record_hidden_alleles`: Record hidden alleles in output (Default: FALSE).
#'   \item ... and other parameters like `rhat_threshold`, `ess_threshold`, etc.
#' }
#'
#' @returns A data frame summarizing the posterior probability of recrudescence for
#'   each patient. Also writes a summary CSV file to the output folder.
#'
#' @export
#'
#' @examples
#' # First, get the path to the example data file included with the package
#' example_file <- system.file("extdata", "Haplotype_data.xlsx",
#'                             package = "MalReBay")
#'
#' # Create a temporary directory to save the results
#' temp_output_dir <- tempdir()
#'
#' # Define a minimal MCMC configuration for a quick example run
#' quick_mcmc_config <- list(
#'   n_chains = 2,
#'   max_iterations = 500,
#'   chunk_size = 500
#' )
#'
#' \dontrun{
#' # Run the full analysis
#' final_summary <- classify_infections(
#'   input_filepath = example_file,
#'   mcmc_config = quick_mcmc_config,
#'   output_folder = temp_output_dir
#' )
#'
#' # View the results
#' print(head(final_summary))
#' }
#'
classify_infections <- function(input_filepath,
                                mcmc_config,
                                output_folder) {
  
  # Data import function
  imported_data <- import_data(filepath = input_filepath)
  genotypedata_latefailures <- imported_data$late_failures
  additional_genotypedata <- imported_data$additional
  marker_information <- imported_data$marker_info
  data_type <- imported_data$data_type
  
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
    if (length(missing_markers) > 0) {
      warning("The following requested markers were NOT found and will be ignored: ",
              paste(missing_markers, collapse = ", "))
    }

    extra_markers_in_data <- setdiff(available_base_markers, requested_base_markers)
    if (length(extra_markers_in_data) > 0) {
      message("INFO: The following markers were found in your data but are not on the package's predefined list: ",
              paste(extra_markers_in_data, collapse = ", "))
    }
    
    message("Proceeding with analysis for these markers: ", paste(markers_to_use, collapse = ", "))

  } else if (data_type == "ampseq") {
    data_marker_columns <- grep("_allele_\\d+$", colnames(genotypedata_latefailures), value = TRUE)
    available_base_markers <- unique(gsub("_allele_\\d+$", "", data_marker_columns))
    requested_base_markers <- marker_information$marker_id
    markers_to_use <- intersect(requested_base_markers, available_base_markers)

    if (length(markers_to_use) == 0) {
      stop("None of the markers from the ampseq data were found in the marker_information sheet.")
    }
    message("Amplicon sequencing data detected. The following haplotype markers will be used for analysis: ",
            paste(available_base_markers, collapse = ", "))
  }else {
    stop("Unknown data type detected: ", data_type)
  }
  
  final_marker_info <- marker_information[marker_information$marker_id %in% markers_to_use, ]
  marker_pattern <- paste0("^(", paste(markers_to_use, collapse = "|"), ")(_allele_\\d+|_\\d+)$")
  final_marker_cols <- grep(marker_pattern, colnames(genotypedata_latefailures), value = TRUE)
  id_cols <- setdiff(colnames(genotypedata_latefailures), grep("_allele_\\d+|_\\d+$", colnames(genotypedata_latefailures), value = TRUE))
  latefailures_subset <- genotypedata_latefailures[, c(id_cols, final_marker_cols)]
  additional_subset <- additional_genotypedata[, c(id_cols, final_marker_cols)]
  
  message("The model is running MCMC analysis...")
  results <- run_all_sites(
    genotypedata_latefailures = latefailures_subset,
    additional_genotypedata = additional_subset,
    marker_info_subset = final_marker_info,
    mcmc_config = config,
    data_type = data_type, 
    output_folder = output_folder
  )
  
  
  cat("\n--- Debug: Structure of MCMC results ---\n")
  print(utils::str(results))
  cat("----------------------------------------\n\n")
  
  if (is.null(results) || length(results$ids) == 0) {
    warning("MCMC results are empty. No summary table can be generated.")
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
      cat(sprintf("INFO: No valid patient data for probabilities for site '%s'. Skipping.\n", site_name))
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
     warning("No probability data was successfully summarized across all sites.")
     return(data.frame(Site=character(), Sample.ID=character(), Probability=numeric(), 
                       N_Available_D0=integer(), N_Available_DF=integer(), N_Comparable_Loci=integer()))
  }
  
  probabilities_df <- dplyr::bind_rows(summary_list_probs)
  locus_summary_df <- dplyr::bind_rows(results$locus_summary, .id = "Site")
  locus_summary_df <- dplyr::rename(locus_summary_df, 
                                    Sample.ID = patient_id,
                                    N_Available_D0 = n_available_d0,
                                    N_Available_DF = n_available_df,
                                    N_Comparable_Loci = n_comparable_loci)

  # Use `summary_df` as the final variable name for the main summary
  summary_df <- dplyr::left_join(probabilities_df, 
                                 locus_summary_df, 
                                 by = c("Sample.ID", "Site"))
  
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  csv_path <- file.path(output_folder, "probability_of_recrudescence_summary.csv")
  utils::write.csv(summary_df, csv_path, row.names = FALSE)
  message("Summary table saved to: ", csv_path)

  # --- MARKER DETAILS DATAFRAME CREATION ---
  marker_details_list <- list()
  for (site_name in names(results$ids)) {
    if (is.null(results$locus_lrs[[site_name]])) next 
    
    site_ids <- results$ids[[site_name]]
    site_lrs <- results$locus_lrs[[site_name]]
    site_dists <- results$locus_dists[[site_name]]
    
    if (!is.null(results$locinames) && !is.null(results$locinames[[site_name]])) {
      marker_names <- results$locinames[[site_name]]
    } else {
      marker_names <- dimnames(site_lrs)[[2]]
      if (is.null(marker_names)) {
        marker_names <- paste0("Locus_", 1:ncol(site_lrs))
      }
    }

    mean_lrs <- apply(site_lrs, c(1, 2), mean, na.rm = TRUE)
    mean_dists <- apply(site_dists, c(1, 2), mean, na.rm = TRUE)
    rownames(mean_lrs) <- rownames(mean_dists) <- site_ids
    colnames(mean_lrs) <- colnames(mean_dists) <- marker_names
    lrs_df <- as.data.frame(as.table(mean_lrs), stringsAsFactors = FALSE)
    colnames(lrs_df) <- c("Sample.ID", "Marker", "Mean_Likelihood_Ratio")
    dists_df <- as.data.frame(as.table(mean_dists), stringsAsFactors = FALSE)
    colnames(dists_df) <- c("Sample.ID", "Marker", "Mean_Distance")
    site_marker_details <- dplyr::full_join(lrs_df, dists_df, by = c("Sample.ID", "Marker"))
    site_marker_details$Site <- site_name
    marker_details_list[[site_name]] <- site_marker_details
  }

  if (length(marker_details_list) > 0) {
    marker_details_df <- dplyr::bind_rows(marker_details_list)
    marker_details_df$Interpretation <- ifelse(marker_details_df$Mean_Likelihood_Ratio > 1, "Supports Recrudescence", "Supports Reinfection/Error")
    marker_csv_path <- file.path(output_folder, "marker_level_summary.csv")
    utils::write.csv(marker_details_df, marker_csv_path, row.names = FALSE)
    message("Detailed marker-level summary saved to: ", marker_csv_path)
  }
  
return(
    list(
        summary = summary_df,
        marker_details = if (exists("marker_details_df")) marker_details_df else NULL
    )
)
}