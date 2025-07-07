#' @title Classify infections using a Bayesian MCMC Framework
#' @description
#'  This is the main analysis function for the BayesMalariaTES package. It analyzes  
#'  microsatellite data from paired day 0 and day of failure samples to differentiate 
#'  between recrudescent and new infections.The function processes the entire 
#'  workflow, including data import, MCMC simulation, convergence assessment, 
#'  and summarization of results.
#' 
#' 
#' @param input_data_path the excel file path containing the genotyping data.
#' @param marker_information a data frame with marker names and their repeat lengths.
#' @param mcmc_config a list containing MCMC configuration parameters.
#' @param output_folder the folder where the results will be saved.
#'
#' @returns A data frame summarizing the posterior probability of recrudescence for
#'  each sample. Also writes the summary to `probability_of_recrudescence_summary.csv` 
#'  in the output folder. 
#'  
#' @export
#'
#' @examples
#' # First, get the path to the example data file included with the package
#' example_file <- system.file("extdata", "Angola_2021_TES_7NMS.xlsx",
#'                             package = "BayesMalariaTES")
#'
#' # Create a temporary directory to save the results
#' temp_output_dir <- tempdir()
#'
#' # Define the marker information for the analysis
#' marker_info <- data.frame(
#'   Marker = c("313", "383", "TA1", "POLYA", "PFPK2", "2490"),
#'   RepeatLength = c(2, 2, 3, 3, 3, 3)
#' )
#'
#' # Define a minimal MCMC configuration for a quick example run
#' quick_mcmc_config <- list(
#'   n_chains = 2,
#'   max_iterations = 500,
#'   burn_in_frac = 0.2
#' )
#'
#' \dontrun{
#' # Run the full analysis
#' final_summary <- classify_infections(
#'   input_data_path = example_file,
#'   marker_information = marker_info,
#'   mcmc_config = quick_mcmc_config,
#'   output_folder = temp_output_dir
#' )
#'
#' # View the results
#' print(head(final_summary))
#' }



classify_infections <- function(input_data_path, marker_information, mcmc_config, output_folder) {
  
  # Set up MCMC Configuration with Defaults
  defaults <- list(
    n_chains = 4, rhat_threshold = 1.01, ess_threshold = 400,
    chunk_size = 2000, max_iterations = 10000, burn_in_frac = 0.25
  )
  config <- utils::modifyList(defaults, mcmc_config)
  
  # Importing data
  imported_data <- import_and_clean_data(filepath = input_data_path)
  genotypedata_latefailures <- imported_data$late_failures
  additional_genotypedata <- imported_data$additional
  
  # validating the markers available in the data
  message("Validating and filtering markers...")
  data_marker_columns <- grep("_\\d+$", colnames(genotypedata_latefailures), value = TRUE)
  available_base_markers <- unique(gsub("_\\d+$", "", data_marker_columns))
  requested_base_markers <- marker_information$Marker
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
  
  final_marker_info <- marker_information[marker_information$Marker %in% markers_to_use, ]
  final_locirepeats <- final_marker_info$RepeatLength
  names(final_locirepeats) <- final_marker_info$Marker
  
  id_cols <- setdiff(colnames(genotypedata_latefailures), data_marker_columns)
  cols_to_keep_pattern <- paste(markers_to_use, collapse = "|")
  final_marker_cols <- grep(cols_to_keep_pattern, data_marker_columns, value = TRUE)
  
  latefailures_subset <- genotypedata_latefailures[, c(id_cols, final_marker_cols)]
  additional_subset <- additional_genotypedata[, c(id_cols, final_marker_cols)]
  
  # Running MCMC Analysis
  message("Running MCMC analysis...")
  results <- run_all_sites(
    genotypedata_latefailures = latefailures_subset,
    additional_genotypedata = additional_subset,
    locirepeats = final_locirepeats,
    n_chains = config$n_chains,
    R_hat_threshold = config$rhat_threshold,
    ESS_threshold = config$ess_threshold,
    chunk_size = config$chunk_size,
    max_iterations = config$max_iterations,
    burn_in_frac = config$burn_in_frac
  )
  
  if (length(results$ids) == 0) stop("No MCMC results were generated.")
  
  # save results
  message("Creating summary of the results")
  summary_list <- lapply(names(results$ids), function(site) {
    site_ids <- results$ids[[site]]
    site_probs <- results$classifications[[site]]
    if (is.null(site_ids) || is.null(site_probs)) return(NULL)
    data.frame(Site = site, ID = site_ids, Probability = rowMeans(site_probs))
  })
  summary_df <- do.call(rbind, Filter(Negate(is.null), summary_list))
  
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  csv_path <- file.path(output_folder, "probability_of_recrudescence_summary.csv")
  utils::write.csv(summary_df, csv_path, row.names = FALSE)
  message("Summary table saved to: ", csv_path)
  
  return(summary_df)
}