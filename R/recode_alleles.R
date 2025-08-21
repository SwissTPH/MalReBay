#' Recode a Single Allele Value to its Bin Index
#'
#' @description
#' Categorizes a single, raw, numeric allele value (e.g., a fragment length)
#' by finding the allele bin with the closest center point.
#'
#' @details
#' This function calculates the center of each provided allele bin and identifies
#' the bin whose center is numerically closest to the `proposed` allele value. An
#' optional `max_distance_allowed` can be used to invalidate matches that are
#' too far from any bin center, which is useful for handling outliers or clear
#' genotyping errors.
#'
#' @param alleles_definitions_subset A two-column matrix (`lower`, `upper`)
#'   defining the allele bins for a *single* genetic marker.
#' @param proposed A single, raw numeric allele value to be categorized.
#' @param max_distance_allowed A numeric threshold. If the absolute distance
#'   between `proposed` and the closest bin center exceeds this value,
#'   `NA_integer_` is returned. Default is `Inf` (no limit).
#'
#' @return An integer representing the row index of the best-matching bin in
#'   `alleles_definitions_subset`. Returns `NA_integer_` if the input is `NA`
#'   or if the match is invalidated by `max_distance_allowed`.
#'
#' @keywords internal
#' @noRd
#'
#'
recodeallele <- function(alleles_definitions_subset, proposed, max_distance_allowed = Inf) {
  if (is.na(proposed) || is.null(alleles_definitions_subset) || nrow(alleles_definitions_subset) == 0) {
    return(NA_integer_)
  }

  bin_centers <- rowMeans(alleles_definitions_subset, na.rm = TRUE)
  closest_bin_index <- which.min(abs(proposed - bin_centers))
  min_distance <- abs(proposed - bin_centers[closest_bin_index])
  if (min_distance > max_distance_allowed) {
    return(NA_integer_)
  }

  return(closest_bin_index)
}


#' Recode and Summarize an Entire Genotyping Dataset
#'
#' @description
#' Processes a full genotyping dataset, recoding raw allele values into their
#' discrete representations based on the marker type, and then summarizes the
#' results for each patient into a standardized format.
#'
#' @details
#' This function serves as a major pre-processing step. Its workflow is as follows:
#' \enumerate{
#'   \item It first calculates the Multiplicity of Infection (MOI) for each
#'     patient's Day 0 and Day of Failure samples.
#'   \item It then iterates through each genetic locus, applying the correct
#'     recoding strategy based on the `binning_method` in `marker_info`:
#'     \itemize{
#'       \item **`microsatellite`** or **`cluster`**: Calls the `recodeallele`
#'         helper for each observed allele to convert its raw value into an
#'         integer bin index.
#'       \item **`exact`**: Uses the raw allele values (e.g., character sequences
#'         from ampseq data) directly without transformation.
#'     }
#'   \item Finally, for each patient at each locus, it creates a summary string
#'     of the unique, sorted alleles in the format "D0-allele1-D0-allele2/DF-allele1".
#' }
#'
#' @param genotypedata A data frame containing the combined genotyping data for all samples.
#' @param alleles_definitions A named list where each element corresponds to a
#'   locus and contains the allele bin definitions matrix for that locus.
#' @param marker_info A data frame containing marker metadata, used to determine
#'   the `binning_method` for each locus.
#'
#' @return A list containing two elements:
#' \item{observeddatamatrix}{A list where each element is a named locus. Each
#'   of these contains a character vector (one entry per patient) summarizing
#'   their recoded Day 0 and Day of Failure alleles in the "D0/DF" format.}
#' \item{MOI}{A matrix with patient IDs as rownames, containing the calculated
#'   MOI for Day 0 (`MOI0`) and Day of Failure (`MOIf`).}
#'
#' @keywords internal
#' @noRd
#'
recode_alleles <- function(genotypedata, alleles_definitions, marker_info) {
  ids <- unique(gsub(" Day 0", "", genotypedata$Sample.ID[grepl(" Day 0", genotypedata$Sample.ID)]))
  locinames <- names(alleles_definitions)
  nids <- length(ids)
  nloci <- length(locinames)
  
  MOI0 <- rep(0, nids)
  MOIf <- rep(0, nids)
  for (i in 1:nids) {
    for (j in 1:nloci) {
      locicolumns <- grepl(paste0("^", locinames[j], "_"), colnames(genotypedata))
      nalleles0 <- sum(!is.na(genotypedata[grepl(paste(ids[i], "Day 0"), genotypedata$Sample.ID), locicolumns]))
      nallelesf <- sum(!is.na(genotypedata[grepl(paste(ids[i], "Day Failure"), genotypedata$Sample.ID), locicolumns]))
      MOI0[i] <- max(MOI0[i], nalleles0)
      MOIf[i] <- max(MOIf[i], nallelesf)
    }
  }
  
  observeddatamatrix <- list()
  for (locus_name in locinames) {
    
    # Get the binning method for this specific locus
    locus_marker_info <- marker_info[marker_info$marker_id == locus_name, ]
    binning_method <- locus_marker_info$binning_method[1]
    
    # Get the raw allele data for this locus
    locus_cols_pattern <- paste0("^", locus_name, "_")
    locus_cols_logical <- grepl(locus_cols_pattern, colnames(genotypedata))
    raw_alleles_matrix <- as.matrix(genotypedata[, locus_cols_logical, drop = FALSE])
    
    recoded_matrix <- NULL
    if (binning_method %in% c("microsatellite", "cluster")) {
      current_definitions <- alleles_definitions[[locus_name]]
      recoded_matrix <- apply(raw_alleles_matrix, MARGIN = c(1, 2), FUN = function(cell_value) {
        recodeallele(alleles_definitions_subset = current_definitions, proposed = cell_value)
      })
      
    } else if (binning_method == "exact") {
      recoded_matrix <- raw_alleles_matrix
      
    } else {
      stop("Unknown binning_method '", binning_method, "' for locus: ", locus_name)
    }
    
    if (is.null(dim(recoded_matrix))) {
      recoded_matrix <- matrix(recoded_matrix, nrow = nrow(raw_alleles_matrix))
    }
    
    tempobservedata <- character(nids)
    for (i_patient_idx in 1:nids) {
      patient_id_str <- ids[i_patient_idx]
      day0_logical_idx <- grepl(paste(patient_id_str, "Day 0"), genotypedata$Sample.ID)
      dayf_logical_idx <- grepl(paste(patient_id_str, "Day Failure"), genotypedata$Sample.ID)
      
      day0_alleles <- recoded_matrix[day0_logical_idx, , drop = FALSE]
      day0_unique_sorted <- sort(unique(day0_alleles[!is.na(day0_alleles)]))
      
      dayf_alleles <- recoded_matrix[dayf_logical_idx, , drop = FALSE]
      dayf_unique_sorted <- sort(unique(dayf_alleles[!is.na(dayf_alleles)]))
      
      day0_str <- paste(day0_unique_sorted, collapse = "-")
      dayf_str <- paste(dayf_unique_sorted, collapse = "-")
      tempobservedata[i_patient_idx] <- paste(day0_str, dayf_str, sep = "/")
    }
    
    observeddatamatrix[[locus_name]] <- tempobservedata
  }
  
  MOItemp <- cbind(MOI0, MOIf)
  rownames(MOItemp) <- ids
  list(observeddatamatrix = observeddatamatrix, MOI = MOItemp)
}