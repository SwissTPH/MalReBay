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

perform_match_counting <- function(genotypedata_latefailures, marker_info) {
  
  # 1. Robustly generate patient IDs from the data
  ids <- unique(trimws(gsub(" .*", "", genotypedata_latefailures$Sample.ID)))
  
  # 2. Identify all loci from the column names
  marker_cols <- grep("(_allele_\\d+|_\\d+)$", colnames(genotypedata_latefailures), value = TRUE)
  locinames <- unique(gsub("(_allele_\\d+|_\\d+)$", "", marker_cols))
  
  nloci <- length(locinames)
  nids <- length(ids)
  
  # 3. Validate that all found loci are defined in the marker metadata
  if (!all(locinames %in% marker_info$marker_id)) {
    missing_loci <- locinames[!locinames %in% marker_info$marker_id]
    warning("Loci in data but not in marker_info (will be ignored): ", paste(missing_loci, collapse=", "))
    locinames <- intersect(locinames, marker_info$marker_id)
    nloci <- length(locinames)
  }
  
  # 4. Initialize output structures
  number_matches <- rep(NA, nids)
  match_output <- matrix("", nids, nloci) 
  colnames(match_output) <- locinames
  number_loci_compared <- rep(NA, nids)
  
  # Create lookup tables for efficiency
  type_lookup <- stats::setNames(marker_info$markertype, marker_info$marker_id)
  bin_lookup <- stats::setNames(marker_info$repeatlength, marker_info$marker_id)
  bin_method_lookup <- stats::setNames(marker_info$binning_method, marker_info$marker_id)
  cluster_thresh_lookup <- stats::setNames(marker_info$cluster_gap_threshold, marker_info$marker_id)
  
  # 5. Loop through each unique patient
  for (i in 1:nids) {
    nmatches_temp <- 0
    nloci_temp <- 0
    current_id <- ids[i]
    
    # Find the row indices for this patient's Day 0 and Day Failure entries
    day0_string_to_find <- paste(current_id, "Day 0")
    dayf_string_to_find <- paste(current_id, "Day Failure")
    
    day0_idx <- which(trimws(genotypedata_latefailures$Sample.ID) == day0_string_to_find)
    dayf_idx <- which(trimws(genotypedata_latefailures$Sample.ID) == dayf_string_to_find)
    
    # Critical Safety Check: Ensure a unique D0/DF pair exists
    if (length(day0_idx) != 1 || length(dayf_idx) != 1) {
      warning(paste("Could not find a unique Day 0/DF pair for patient", current_id, ". Skipping."))
      match_output[i, ] <- "ERR" 
      next
    }
    
    # Loop through each locus for this patient
    for (j in 1:nloci) {
      current_locus <- locinames[j]
      current_marker_type <- type_lookup[current_locus]
      locus_cols <- grep(paste0("^", current_locus, "(_allele_\\d+|_\\d+)$"), colnames(genotypedata_latefailures), value = TRUE)
      method <- bin_method_lookup[current_locus]
      
      # Extract alleles using the validated row indices
      if (!is.na(method) && (method == "microsatellite" || method == "cluster")) {
        day0_alleles <- as.numeric(unlist(genotypedata_latefailures[day0_idx, locus_cols]))
        dayf_alleles <- as.numeric(unlist(genotypedata_latefailures[dayf_idx, locus_cols]))
      } else {
        day0_alleles <- as.character(unlist(genotypedata_latefailures[day0_idx, locus_cols]))
        dayf_alleles <- as.character(unlist(genotypedata_latefailures[dayf_idx, locus_cols]))
      }
      
      # Clean out any NA values
      day0_alleles <- day0_alleles[!is.na(day0_alleles)]
      dayf_alleles <- dayf_alleles[!is.na(dayf_alleles)]
      
      # Perform comparison only if data exists for both timepoints
      if (length(day0_alleles) > 0 && length(dayf_alleles) > 0) {
        nloci_temp <- nloci_temp + 1
        is_match <- FALSE 
        
        
        if (is.na(method)) {
          warning("No binning method found for locus: ", current_locus, ". Skipping comparison.")
          match_output[i, j] <- "NO_INFO"
          next 
        }
        
        if (method == "microsatellite") {
          current_bin_size <- bin_lookup[current_locus]
          d0_in_df <- all(sapply(day0_alleles, function(d0) any(abs(d0 - dayf_alleles) <= current_bin_size)))
          df_in_d0 <- all(sapply(dayf_alleles, function(df) any(abs(df - day0_alleles) <= current_bin_size)))
          if (d0_in_df || df_in_d0) {
            is_match <- TRUE
          }
        } else if (method == "cluster") {
          current_thresh <- cluster_thresh_lookup[current_locus]
          all_locus_alleles <- unique(c(day0_alleles, dayf_alleles))
          cluster_map <- assign_clusters(all_locus_alleles, current_thresh)
          d0_clusters <- if(length(day0_alleles) > 0) unique(cluster_map[as.character(day0_alleles)]) else character(0)
          df_clusters <- if(length(dayf_alleles) > 0) unique(cluster_map[as.character(dayf_alleles)]) else character(0)
          
          if (all(d0_clusters %in% df_clusters) || all(df_clusters %in% d0_clusters)) {
            is_match <- TRUE
          }
        } else {
          if (all(day0_alleles %in% dayf_alleles) || all(dayf_alleles %in% dayf_alleles)) {
            is_match <- TRUE
          }
        }
        
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
  
  # 6. Assemble and return the final results data frame
  match_results_df <- as.data.frame(match_output)
  match_results_df$Sample.ID <- ids
  match_results_df$Number_Matches <- number_matches
  match_results_df$Number_Loci_Compared <- number_loci_compared
  match_results_df <- match_results_df[, c("Sample.ID", "Number_Matches", "Number_Loci_Compared", locinames)]
  
  return(match_results_df)
}