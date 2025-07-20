#' Perform match-counting algorithm
#'
#' @param genotypedata_latefailures The processed late failures data frame.
#' @param marker_info The marker information data frame.
#' @return A data frame summarizing match-counting results.
#' @noRd
#' 
#' 
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

      # Extract alleles using the validated row indices
      if (current_marker_type == "microsatellite") {
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

        if (current_marker_type == "microsatellite") {
          current_bin_size <- bin_lookup[current_locus]
          all_pairs <- expand.grid(day0_alleles, dayf_alleles)
          min_distance <- min(abs(all_pairs$Var1 - all_pairs$Var2), na.rm = TRUE)
          if (is.finite(min_distance) && min_distance <= current_bin_size) {
            is_match <- TRUE
          }
        } else { 
          if (any(day0_alleles %in% dayf_alleles)) {
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