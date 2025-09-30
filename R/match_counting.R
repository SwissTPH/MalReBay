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
  
  # Universal ID parsing that handles both formats
  id_pattern <- "(D[0-9]+| Day 0| Day Failure)$"
  
  paired_data <- genotypedata_latefailures %>%
    dplyr::mutate(
      Patient.ID = trimws(gsub(id_pattern, "", .data$Sample.ID)),
      Day = ifelse(grepl("D0|Day 0", .data$Sample.ID), "Day 0", "Day X")
    )
  
  patient_ids <- unique(paired_data$Patient.ID)
  nids <- length(patient_ids)
  
  marker_cols <- grep("(_allele_\\d+|_\\d+)$", colnames(paired_data), value = TRUE)
  locinames <- unique(gsub("(_allele_\\d+|_\\d+)$", "", marker_cols))
  nloci <- length(locinames)
  
  number_matches <- rep(NA, nids)
  match_output <- matrix("", nids, nloci) 
  colnames(match_output) <- locinames
  number_loci_compared <- rep(NA, nids)
  bin_method_lookup <- stats::setNames(marker_info$binning_method, marker_info$marker_id)
  
  for (i in 1:nids) {
    current_patient_id <- patient_ids[i]
    
    day0_row <- paired_data %>% dplyr::filter(.data$Patient.ID == current_patient_id, .data$Day == "Day 0")
    dayf_row <- paired_data %>% dplyr::filter(.data$Patient.ID == current_patient_id, .data$Day == "Day X")
    
    if (nrow(day0_row) != 1 || nrow(dayf_row) != 1) {
      next
    }
    
    nmatches_temp <- 0
    nloci_temp <- 0
    
    for (j in 1:nloci) {
      current_locus <- locinames[j]
      locus_cols <- grep(paste0("^", current_locus, "(_allele_\\d+|_\\d+)$"), colnames(paired_data), value = TRUE)
      
      day0_alleles <- as.character(unlist(day0_row[, locus_cols]))
      dayf_alleles <- as.character(unlist(dayf_row[, locus_cols]))
      day0_alleles <- day0_alleles[!is.na(day0_alleles)]
      dayf_alleles <- dayf_alleles[!is.na(dayf_alleles)]
      
      if (length(day0_alleles) > 0 && length(dayf_alleles) > 0) {
        nloci_temp <- nloci_temp + 1
        is_match <- all(dayf_alleles %in% day0_alleles)
        
        if (is_match) {
          nmatches_temp <- nmatches_temp + 1
          match_output[i, j] <- "R"
        } else {
          match_output[i, j] <- "NI"
        }
      } else {
        match_output[i, j] <- "IND"
      }
    }
    number_matches[i] <- nmatches_temp
    number_loci_compared[i] <- nloci_temp
  }
  
  match_results_df <- as.data.frame(match_output)
  match_results_df$Sample.ID <- patient_ids
  match_results_df$Number_Matches <- number_matches
  match_results_df$Number_Loci_Compared <- number_loci_compared
  match_results_df <- match_results_df[, c("Sample.ID", "Number_Matches", "Number_Loci_Compared", locinames)]
  
  return(match_results_df)
}