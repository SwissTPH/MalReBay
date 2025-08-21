#' Define Allele Bins for Length-Polymorphic Markers
#'
#' @description
#' Groups raw, continuous allele measurements (e.g., fragment lengths) into
#' discrete, well-defined bins. This function is a critical pre-processing step
#' for length-polymorphic data, defining the set of possible alleles (A) for
#' each locus before frequency calculation.
#'
#' @details
#' The function iterates through each locus and applies a binning method as
#' specified in the `marker_info_subset` data frame.
#' \itemize{
#'   \item For the **`cluster`** method, it identifies significant gaps between
#'     the unique sorted allele sizes to define bin boundaries.
#'   \item For the **`microsatellite`** method, it assumes true alleles fall on
#'     a regular grid defined by the repeat length. It identifies the optimal
#'     grid phasing and "snaps" observed values to the nearest expected allele center.
#'   \item Allele definition is skipped for markers with the **`exact`** binning
#'     method (e.g., amplicon data), as their alleles are already discrete.
#' }
#' After initial binning, it can filter the alleles. If `maxk` is specified,
#' only the `maxk` most frequent allele bins (based on observation counts) are
#' retained for each locus.
#'
#' @param genotypedata A data frame containing the combined and cleaned
#'   genotyping data for all samples.
#' @param marker_info_subset A data frame containing the metadata for the
#'   markers to be processed, specifically the `binning_method` and
#'   `repeatlength` columns.
#' @param maxk An integer. The maximum number of allele bins to define per
#'   locus. Bins are prioritized by frequency. Default is `Inf` (no limit).
#'
#' @return A named list, where each element corresponds to a locus. Each element
#'   contains a two-column matrix (`lower`, `upper`) defining the boundaries for
#'   each allele bin for that marker.
#'
#' @keywords internal
#' @noRd
#'
define_alleles <- function(genotypedata, marker_info_subset, maxk = Inf) {
  
  locus_names <- marker_info_subset$marker_id
  final_alleles <- vector("list", length(locus_names))
  names(final_alleles) <- locus_names
  
  for (locus_name in locus_names) {
    marker_details <- marker_info_subset[marker_info_subset$marker_id == locus_name, ]
    current_binning_method <- marker_details$binning_method
    
    message(paste("    -> Binning method is:", current_binning_method))
    
    locus_cols <- grep(paste0("^", locus_name, "_"), colnames(genotypedata), value = TRUE)
    message(paste("    -> Found", length(locus_cols), "columns for this locus."))
    
    if (current_binning_method == "exact") {
      message(paste("INFO: Skipping allele definition for amplicon locus:", locus_name, "(alleles are exact codes)"))
      final_alleles[[locus_name]] <- matrix(NA, ncol = 2, nrow = 0, dimnames = list(NULL, c("lower", "upper")))
      next 
    }
    current_repeat_length <- marker_details$repeatlength
    
    locus_cols <- grep(paste0("^", locus_name, "_"), colnames(genotypedata), value = TRUE)
    if (length(locus_cols) == 0) {
      final_alleles[[locus_name]] <- matrix(NA, ncol = 2, nrow = 0, dimnames = list(NULL, c("lower", "upper")))
      next
    }
    
    raw_alleles <- stats::na.omit(as.vector(as.matrix(genotypedata[, locus_cols, drop = FALSE])))
    if (length(raw_alleles) == 0) {
      final_alleles[[locus_name]] <- matrix(NA, ncol = 2, nrow = 0, dimnames = list(NULL, c("lower", "upper")))
      next
    }
    
    binned_alleles <- matrix(NA, ncol = 3, nrow = 0, dimnames = list(NULL, c("lower", "upper", "count")))
    
    if (current_binning_method == "cluster") {
      current_gap_threshold <- marker_details$cluster_gap_threshold
      unique_observed <- sort(unique(raw_alleles))
      
      if(length(unique_observed) == 1) {
        binned_alleles <- matrix(
          c(unique_observed - 0.5, unique_observed + 0.5, length(raw_alleles)), 
          nrow = 1, dimnames = list(NULL, c("lower", "upper", "count"))
        )
        
      } else {
        gaps <- diff(unique_observed)
        break_points <- which(gaps >= current_gap_threshold)
        cluster_indices_start <- c(1, break_points + 1)
        cluster_indices_end <- c(break_points, length(unique_observed))
        lower_bounds <- unique_observed[cluster_indices_start] - 0.5
        upper_bounds <- unique_observed[cluster_indices_end] + 0.5
        counts <- sapply(seq_along(lower_bounds), function(i) sum(raw_alleles >= lower_bounds[i] & raw_alleles <= upper_bounds[i]))
        binned_alleles <- cbind(lower = lower_bounds, upper = upper_bounds, count = counts)
      }
      
    } else if (current_binning_method == "microsatellite") {
      
      if (diff(range(raw_alleles)) < current_repeat_length) {
        binned_alleles <- matrix(
          c(min(raw_alleles) - current_repeat_length / 2, max(raw_alleles) + current_repeat_length / 2, length(raw_alleles)),
          nrow = 1, dimnames = list(NULL, c("lower", "upper", "count"))
        )
      } else {
        breaks <- seq(floor(min(raw_alleles)) - 0.5, ceiling(max(raw_alleles)) + 0.5, by = 1)
        allele_values <- (breaks[-1] + breaks[-length(breaks)]) / 2
        hist_alleles <- graphics::hist(raw_alleles, breaks = breaks, plot = FALSE)
        counts_by_offset <- sapply(1:current_repeat_length, function(x) sum(hist_alleles$counts[seq(x, length(hist_alleles$counts), by = current_repeat_length)]))
        possible_alleles <- allele_values[seq(which.max(counts_by_offset), length(allele_values), by = current_repeat_length)]
        if (min(raw_alleles) <= (min(possible_alleles) - current_repeat_length / 2)) possible_alleles <- c(min(possible_alleles) - current_repeat_length, possible_alleles)
        if (max(raw_alleles) > (max(possible_alleles) + current_repeat_length / 2)) possible_alleles <- c(possible_alleles, max(possible_alleles) + current_repeat_length)
        clusters <- sapply(raw_alleles, function(x) which.min(abs(possible_alleles - x)))
        unique_clusters <- sort(unique(clusters))
        lower_bounds <- possible_alleles[unique_clusters] - current_repeat_length / 2
        upper_bounds <- possible_alleles[unique_clusters] + current_repeat_length / 2
        counts <- sapply(seq_along(lower_bounds), function(i) sum(raw_alleles > lower_bounds[i] & raw_alleles <= upper_bounds[i]))
        binned_alleles <- cbind(lower = lower_bounds, upper = upper_bounds, count = counts)
      }
    } else {
      warning("Unknown binning method for locus: ", locus_name,". Skipping.")
      next
    }
    
    if (nrow(binned_alleles) > 0) {
      current_maxk <- if (!is.null(names(maxk)) && locus_name %in% names(maxk)) {
        maxk[[locus_name]] 
      } else {
        maxk[1]
      }
      if (is.finite(current_maxk) && current_maxk < nrow(binned_alleles)) {
        message(paste("INFO: For locus", locus_name, "-> filtering to the top", current_maxk, "most frequent alleles."))
        sorted_indices <- order(binned_alleles[, "count"], decreasing = TRUE)[1:current_maxk]
        locus_alleles <- binned_alleles[sorted_indices, c("lower", "upper"), drop = FALSE]
      } else {
        locus_alleles <- binned_alleles[, c("lower", "upper"), drop = FALSE]
      }
      final_alleles[[locus_name]] <- locus_alleles[order(locus_alleles[, "lower"]), , drop = FALSE]
      
    } else {
      final_alleles[[locus_name]] <- matrix(NA, ncol = 2, nrow = 0, dimnames = list(NULL, c("lower", "upper")))
    }
  }
  
  return(final_alleles)
}