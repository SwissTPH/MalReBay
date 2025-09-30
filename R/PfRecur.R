#' Run PfRecur analysis with specific data preparation
#'
#' @description Replicates the data preparation for a specific dataset structure
#' (e.g., Angola TES) and then runs the PfRecur analysis.
#'
#' @param raw_data_df The raw data frame from the input Excel file.
#' @param output_csv_path Path for the final results CSV file.
#' @param epsilon The error rate parameter for PfRecur.
#' @param omega_vals The omega parameter for PfRecur.
#' @param beta The beta parameter for PfRecur.
#' @return A data frame with the PfRecur analysis results.
#' @export
run_pfrecur_analysis_original_prep <- function(raw_data_df,
                                               output_csv_path = "pfrecur_analysis_results.csv",
                                               epsilon = 0.05,
                                               omega_vals= 0.75,
                                               beta = 0.25) {
  
  # Preparing data
  
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr package is required.")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr package is required.")
  if (!requireNamespace("stringr", quietly = TRUE)) stop("stringr package is required.")
  
  ANGOLA_MARKERS <- c("TA109", "M313", "M383", "TA1", "POLYA", "PFPK2", "M2490")
  df_long <- raw_data_df %>%
    dplyr::rename(SampleID_raw = `Sample.ID`) %>%
    tidyr::extract(SampleID_raw, c("PatientID", "Timepoint"), "([A-Z]{2}[0-9-]{6,})(D[0-9]+)", remove = FALSE) %>%
    tidyr::unite(TA109, dplyr::all_of(colnames(raw_data_df)[grepl("TA109", colnames(raw_data_df))]), sep="/", na.rm=TRUE) %>%
    tidyr::unite(M313,  dplyr::all_of(colnames(raw_data_df)[grepl("313",   colnames(raw_data_df))]), sep="/", na.rm=TRUE) %>%
    tidyr::unite(M383,  dplyr::all_of(colnames(raw_data_df)[grepl("383",   colnames(raw_data_df))]), sep="/", na.rm=TRUE) %>%
    tidyr::unite(TA1,   dplyr::all_of(colnames(raw_data_df)[grepl("TA1_",  colnames(raw_data_df))]), sep="/", na.rm=TRUE) %>%
    tidyr::unite(POLYA, dplyr::all_of(colnames(raw_data_df)[grepl("POLYA", colnames(raw_data_df))]), sep="/", na.rm=TRUE) %>%
    tidyr::unite(PFPK2, dplyr::all_of(colnames(raw_data_df)[grepl("PFPK2", colnames(raw_data_df))]), sep="/", na.rm=TRUE) %>%
    tidyr::unite(M2490, dplyr::all_of(colnames(raw_data_df)[grepl("2490",  colnames(raw_data_df))]), sep="/", na.rm=TRUE)
  
  df_long <- as.data.frame(df_long)
  rownames(df_long) <- df_long$SampleID_raw
  
  my_marker_set <- lapply(ANGOLA_MARKERS, function(marker_name) {
    all_alleles <- df_long[[marker_name]][df_long[[marker_name]] != ""]
    all_split_alleles <- unlist(strsplit(all_alleles, "/"))
    clean_split_alleles <- all_split_alleles[all_split_alleles != "NA"]
    unique_alleles <- as.character(sort(as.numeric(unique(clean_split_alleles))))
    return(unique_alleles)
  })
  names(my_marker_set) <- ANGOLA_MARKERS
  
  my_genotype_matrix <- lapply(ANGOLA_MARKERS, function(marker_name) {
    alleles_for_this_marker <- my_marker_set[[marker_name]]
    if (length(alleles_for_this_marker) == 0) return(NULL)
    genotype_mat <- matrix(0,
                           nrow = nrow(df_long),
                           ncol = length(alleles_for_this_marker),
                           dimnames = list(df_long$SampleID_raw, alleles_for_this_marker))
    for (i in 1:nrow(df_long)) {
      sample_alleles_str <- df_long[[marker_name]][i]
      if (sample_alleles_str != "" && !is.na(sample_alleles_str)) {
        alleles_detected <- unlist(strsplit(sample_alleles_str, "/"))
        alleles_detected <- alleles_detected[!is.na(alleles_detected)]
        valid_alleles_to_fill <- intersect(alleles_detected, colnames(genotype_mat))
        if (length(valid_alleles_to_fill) > 0) {
          genotype_mat[df_long$SampleID_raw[i], valid_alleles_to_fill] <- 1
        }
      }
    }
    return(genotype_mat)
  })
  my_genotype_matrix <- my_genotype_matrix[!sapply(my_genotype_matrix, is.null)]
  names(my_genotype_matrix) <- names(my_marker_set)[sapply(my_marker_set, function(x) length(x) > 0)]
  
  day0_samples <- df_long %>% dplyr::filter(Timepoint == "D0")
  recurrent_samples <- df_long %>% dplyr::filter(Timepoint != "D0")
  paired_data <- dplyr::inner_join(
    day0_samples %>% dplyr::select(PatientID, Site, D0_SampleID = SampleID_raw),
    recurrent_samples %>% dplyr::select(PatientID, Recurrent_SampleID = SampleID_raw),
    by = "PatientID"
  )
  day0_isolates_by_site <- split(day0_samples$SampleID_raw, day0_samples$Site)
  my_isolates <- list()
  for (i in 1:nrow(paired_data)) {
    current_pair <- paired_data[i, ]
    current_site <- current_pair$Site
    if (!is.null(current_site) && !is.na(current_site) && current_site %in% names(day0_isolates_by_site)) {
      my_isolates[[i]] <- list(
        recurrent = current_pair$Recurrent_SampleID,
        ref_C = current_pair$D0_SampleID,
        ref_I = setdiff(day0_isolates_by_site[[current_site]], current_pair$D0_SampleID)
      )
    }
  }
  my_isolates <- my_isolates[!sapply(my_isolates, is.null)]
  
  if (!requireNamespace("PfRecur", quietly = TRUE)) stop("PfRecur package is required.")
  if (!requireNamespace("parallel", quietly = TRUE)) stop("parallel package is required.")
  if (!requireNamespace("purrr", quietly = TRUE)) stop("purrr package is required.")
  genotype_matrix_clean <- lapply(my_genotype_matrix, function(marker_matrix) {
    row_sums <- rowSums(marker_matrix)
    clean_matrix <- marker_matrix[row_sums > 0, , drop = FALSE]
    return(clean_matrix)
  })
  
  my_error_matrix <- PfRecur::geometric_error_matrix(my_marker_set[names(genotype_matrix_clean)], EPSILON = epsilon)
  
  full_results <- parallel::mclapply(my_isolates, function(x)
    PfRecur::evaluate_posterior(
      recurrent = x[["recurrent"]],
      ref_C = x[["ref_C"]],
      ref_I = x[["ref_I"]],
      genotype_matrix = genotype_matrix_clean,
      error_matrix = my_error_matrix,
      omega_vals = omega_vals,
      beta = beta,
      keep_markers = "all"
    )
  )
  
  final_table <- purrr::map_dfr(full_results, ~ .x$metrics, .id = "pair_index")
  
  write.csv(final_table, output_csv_path, row.names = FALSE)
  message(paste("PfRecur analysis complete. Results saved to:", output_csv_path))
  
  return(final_table)
}