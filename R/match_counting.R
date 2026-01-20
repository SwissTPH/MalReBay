#' Perform Allele Match-Counting Analysis
#'
#' @description
#' Implements a traditional match-counting algorithm to classify infections.
#' This provides a simpler, deterministic alternative to the Bayesian MCMC
#' framework for assessing whether a follow-up infection is a recrudescence or
#' a new infection.
#'
#' @details
#' This function iterates through each patient with a paired "Day 0" and "Day of
#' Failure" sample. For each genetic marker with data at both time points, it
#' compares the set of alleles and determines if they "match".
#'
#' The definition of a match depends on the marker's `binning_method` as
#' specified in `marker_info`:
#' \itemize{
#'   \item **`microsatellite`**: A match occurs if for every Day 0 allele, there
#'     is a Day of Failure allele within the specified `repeatlength` tolerance,
#'     OR vice-versa.
#'   \item **`cluster`**: Alleles from both time points are first grouped into
#'     clusters based on the `cluster_gap_threshold`. A match occurs if the set
#'     of Day 0 clusters is a subset of the Day of Failure clusters, OR vice-versa.
#'   \item **`exact`** (e.g., for ampseq data): A match occurs if the set of
#'     Day 0 alleles is a subset of the Day of Failure alleles, OR vice-versa.
#' }
#' The function returns a summary data frame with the total number of matches
#' for each patient and a detailed breakdown of the result for each locus.
#'
#' @param genotypedata_latefailures The processed late failures data frame,
#'   containing paired samples for all patients.
#' @param marker_info The marker information data frame, which provides the
#'   `binning_method` and other parameters for each locus.
#' @return A data frame summarizing the match-counting results. The key columns
#'   are:
#'   \item{Sample.ID}{The unique patient identifier.}
#'   \item{Number_Matches}{The total count of loci classified as a match (R).}
#'   \item{Number_Loci_Compared}{The total number of loci with data at both Day 0
#'     and Day of Failure for that patient.}
#'   It also includes one column for each locus, containing a code for the
#'   comparison result:
#'   \itemize{
#'     \item \strong{`R`}: Recrudescence (a match).
#'     \item \strong{`NI`}: New Infection (not a match).
#'     \item \strong{`IND`}: Indeterminate (data was missing for either Day 0 or
#'       Day of Failure).
#'     \item \strong{`ERR`}: An error occurred (e.g., could not find a unique
#'       paired sample for the patient).
#'   }
#'
#' @keywords internal
#' @noRd
#'

assign_clusters <- function(alleles, threshold) {
  if (length(alleles) == 0) return(character(0))
  unique_sorted_alleles <- sort(unique(alleles))
  if (length(unique_sorted_alleles) <= 1) {
    cluster_ids <- rep(1, length(unique_sorted_alleles))
  } else {
    gaps <- diff(unique_sorted_alleles)
    cluster_ids <- cumsum(c(1, gaps > threshold))
  }
  return(stats::setNames(cluster_ids, unique_sorted_alleles))
}

#' @importFrom dplyr mutate filter
#' @importFrom stats setNames

perform_match_counting <- function(genotypedata_latefailures, marker_info) {
  
  # 1. Universal ID parsing (Kept from your current version)
  id_pattern <- "(_?D[0-9]+| Day 0| Day Failure)$"
  
  paired_data <- genotypedata_latefailures %>%
    dplyr::mutate(
      Patient.ID = trimws(gsub(id_pattern, "", .data$Sample.ID)),
      Day = ifelse(grepl("D0|Day 0", .data$Sample.ID), "Day 0", "Day X")
    )
  
  patient_ids <- unique(paired_data$Patient.ID)
  nids <- length(patient_ids)
  
  # Identify markers
  marker_cols <- grep("(_allele_\\d+|_\\d+)$", colnames(paired_data), value = TRUE)
  locinames <- unique(gsub("(_allele_\\d+|_\\d+)$", "", marker_cols))
  
  # Ensure we only use loci present in marker_info
  locinames <- intersect(locinames, marker_info$marker_id)
  nloci <- length(locinames)
  
  # 2. Lookups (Restored from your original version)
  bin_method_lookup <- stats::setNames(marker_info$binning_method, marker_info$marker_id)
  bin_lookup <- stats::setNames(marker_info$repeatlength, marker_info$marker_id)
  cluster_thresh_lookup <- stats::setNames(marker_info$cluster_gap_threshold, marker_info$marker_id)
  
  # Initialize output
  number_matches <- rep(NA, nids)
  match_output <- matrix("", nids, nloci) 
  colnames(match_output) <- locinames
  number_loci_compared <- rep(NA, nids)
  
  for (i in 1:nids) {
    current_patient_id <- patient_ids[i]
    
    day0_row <- paired_data %>% dplyr::filter(.data$Patient.ID == current_patient_id, .data$Day == "Day 0")
    dayf_row <- paired_data %>% dplyr::filter(.data$Patient.ID == current_patient_id, .data$Day == "Day X")
    
    # Check for valid pair
    if (nrow(day0_row) != 1 || nrow(dayf_row) != 1) {
      match_output[i, ] <- "ERR"
      next
    }
    
    nmatches_temp <- 0
    nloci_temp <- 0
    
    for (j in 1:nloci) {
      current_locus <- locinames[j]
      locus_cols <- grep(paste0("^", current_locus, "(_allele_\\d+|_\\d+)$"), colnames(paired_data), value = TRUE)
      
      method <- bin_method_lookup[current_locus]
      
      # Extract alleles
      day0_alleles <- unlist(day0_row[, locus_cols])
      dayf_alleles <- unlist(dayf_row[, locus_cols])
      day0_alleles <- day0_alleles[!is.na(day0_alleles) & day0_alleles != ""]
      dayf_alleles <- dayf_alleles[!is.na(dayf_alleles) & dayf_alleles != ""]
      
      if (length(day0_alleles) > 0 && length(dayf_alleles) > 0) {
        nloci_temp <- nloci_temp + 1
        is_match <- FALSE 
        
        # 3. Logic Branching (Restored and fixed)
        if (is.na(method) || method == "amplicon" || method == "exact") {
          # Standard exact string matching
          is_match <- any(dayf_alleles %in% day0_alleles)
          
        } else if (method == "microsatellite") {
          # Numeric matching with repeat length tolerance
          d0_num <- as.numeric(day0_alleles)
          df_num <- as.numeric(dayf_alleles)
          bin_size <- bin_lookup[current_locus]
          
          # Check if any allele in Day Failure has a match in Day 0 within the bin size
          is_match <- any(sapply(df_num, function(df) any(abs(df - d0_num) <= bin_size)))
          
        } else if (method == "cluster") {
          # Clustering logic
          d0_num <- as.numeric(day0_alleles)
          df_num <- as.numeric(dayf_alleles)
          thresh <- cluster_thresh_lookup[current_locus]
          
          all_vals <- unique(c(d0_num, df_num))
          cluster_map <- assign_clusters(all_vals, thresh)
          
          d0_clusters <- unique(cluster_map[as.character(d0_num)])
          df_clusters <- unique(cluster_map[as.character(df_num)])
          
          # Check if all failure clusters were present at Day 0
          is_match <- any(df_clusters %in% d0_clusters)
        }
        
        if (is_match) {
          nmatches_temp <- nmatches_temp + 1
          match_output[i, j] <- "R"   # Recrudescence
        } else {
          match_output[i, j] <- "NI"  # New Infection
        }
      } else {
        match_output[i, j] <- "IND"   # Indeterminate
      }
    }
    number_matches[i] <- nmatches_temp
    number_loci_compared[i] <- nloci_temp
  }
  
  # Prepare final dataframe
  match_results_df <- as.data.frame(match_output)
  match_results_df$Sample.ID <- patient_ids
  match_results_df$Number_Matches <- number_matches
  match_results_df$Number_Loci_Compared <- number_loci_compared
  
  # Reorder columns
  match_results_df <- match_results_df[, c("Sample.ID", "Number_Matches", "Number_Loci_Compared", locinames)]
  
  return(match_results_df)
}