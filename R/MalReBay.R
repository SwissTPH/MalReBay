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
  
  defaults <- list(n_chains = 4, 
                   rhat_threshold = 1.01, 
                   ess_threshold = 400,    
                   chunk_size = 2000, 
                   max_iterations = 20000, 
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
    message("Proceeding with analysis for these markers: ", paste(markers_to_use, collapse = ", "))
  } else if (data_type == "ampseq") {
    markers_to_use <- marker_information$marker_id
    message("Amplicon sequencing data detected. The following haplotype markers will be used for analysis: ",
            paste(available_base_markers, collapse = ", "))
    markers_to_use <- available_base_markers
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
    # Return an empty dataframe with the correct column names
    return(data.frame(Site=character(), ID=character(), Probability=numeric()))
  }
  
  summary_list <- list() 
  for (i in 1:length(results$ids)) {
    site_name <- names(results$ids)[i]
    if (is.null(site_name)) {
      site_name <- paste0("Site_", i)
    }
    site_ids <- results$ids[[i]]
    site_probs_matrix <- results$classifications[[i]]
    if (is.null(site_ids) || is.null(site_probs_matrix) || length(site_ids) == 0 || nrow(site_probs_matrix) == 0) {
      cat(sprintf("INFO: No valid patient data to summarize for site '%s'. Skipping.\n", site_name))
      next
    }
    if (length(site_ids) != nrow(site_probs_matrix)) {
      warning(sprintf("For site '%s', the number of IDs (%d) does not match the number of rows in the probability matrix (%d). Skipping site.", 
                      site_name, length(site_ids), nrow(site_probs_matrix)))
      next
    }
    
    probabilities <- rowMeans(site_probs_matrix, na.rm = TRUE)
    site_df <- data.frame(
      Site = site_name, 
      Sample.ID = site_ids, 
      Probability = probabilities,
      stringsAsFactors = FALSE
    )
    summary_list[[i]] <- site_df
  }
  
  if (length(summary_list) == 0) {
    warning("No data was successfully summarized across all sites.")
    summary_df <- data.frame(Site=character(), Sample.ID=character(), Probability=numeric())
  } else {
    summary_df <- do.call(rbind, summary_list)
  }
  
  
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  csv_path <- file.path(output_folder, "probability_of_recrudescence_summary.csv")
  utils::write.csv(summary_df, csv_path, row.names = FALSE)
  message("Summary table saved to: ", csv_path)
  
  return(summary_df)
}