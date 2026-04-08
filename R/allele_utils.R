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
#' specified in the `marker_info` data frame.
#' \itemize{
#'   \item For the **`msp_glurp`** method, it identifies significant gaps between
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
#' @param marker_info A data frame containing the metadata for the
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
define_alleles <- function(genotypedata, marker_info, maxk = Inf) {
  
  locus_names <- marker_info$marker_id
  final_alleles <- vector("list", length(locus_names))
  names(final_alleles) <- locus_names
  
  for (locus_name in locus_names) {
    marker_details <- marker_info[marker_info$marker_id == locus_name, ]
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
    
    if (current_binning_method == "msp_glurp") {
      current_gap_threshold <- marker_details$repeatlength
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
#' @param locus_allele_bins A two-column matrix (`lower`, `upper`)
#'   defining the allele bins for a *single* genetic marker.
#' @param proposed A single, raw numeric allele value to be categorized.
#' @param max_distance_allowed A numeric threshold. If the absolute distance
#'   between `proposed` and the closest bin center exceeds this value,
#'   `NA_integer_` is returned. Default is `Inf` (no limit).
#'
#' @return An integer representing the row index of the best-matching bin in
#'   `locus_allele_bins`. Returns `NA_integer_` if the input is `NA`
#'   or if the match is invalidated by `max_distance_allowed`.
#'
#' @keywords internal
#' @noRd
#'
#'
recodeallele <- function(locus_allele_bins, proposed, max_distance_allowed = Inf) {
  if (is.na(proposed) || is.null(locus_allele_bins) || nrow(locus_allele_bins) == 0) {
    return(NA_integer_)
  }
  
  bin_centers <- rowMeans(locus_allele_bins, na.rm = TRUE)
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
#'     patient's Day 0 and Day of recurrence samples.
#'   \item It then iterates through each genetic locus, applying the correct
#'     recoding strategy based on the `binning_method` in `marker_info`:
#'     \itemize{
#'       \item **`microsatellite`** or **`msp_glurp`**: Calls the `recodeallele`
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
#' @param allele_definitions A named list where each element corresponds to a
#'   locus and contains the allele bin definitions matrix for that locus.
#' @param marker_info A data frame containing marker metadata, used to determine
#'   the `binning_method` for each locus.
#'
#' @return A list containing two elements:
#' \item{observeddatamatrix}{A list where each element is a named locus. Each
#'   of these contains a character vector (one entry per patient) summarizing
#'   their recoded Day 0 and Day of recurrence alleles in the "D0/DF" format.}
#' \item{MOI}{A matrix with patient IDs as rownames, containing the calculated
#'   MOI for Day 0 (`MOI0`) and Day of recurrence (`MOIf`).}
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
      nallelesf <- sum(!is.na(genotypedata[grepl(paste(ids[i], "recurrence"), genotypedata$Sample.ID), locicolumns]))
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
    if (binning_method %in% c("microsatellite", "msp_glurp")) {
      current_definitions <- alleles_definitions[[locus_name]]
      recoded_matrix <- apply(raw_alleles_matrix, MARGIN = c(1, 2), FUN = function(cell_value) {
        recodeallele(locus_allele_bins = current_definitions, proposed = cell_value)
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
      dayf_logical_idx <- grepl(paste(patient_id_str, "recurrence"), genotypedata$Sample.ID)
      
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


#' Calculate Allele Frequencies and Genetic Variability
#'
#' @description
#' This function computes the frequency of each defined allele for every locus
#' in the provided genotyping data. It also calculates a measure of genetic
#' variability and counts the number of unique alleles per locus.
#'
#' @details
#' The function processes each locus according to the `binning_method` specified
#' in `marker_info`.
#' \itemize{
#'   \item For **length polymorphism** data (`microsatellite`, `msp_glurp`), it
#'     uses the pre-defined allele bins to tabulate allele counts and calculate
#'     frequencies. Variability is calculated as the mean of the within-bin
#'     standard deviations.
#'   \item For **amplicon sequencing** data (`exact`), it calculates frequencies
#'     based on exact allele matching. Variability is calculated as the expected
#'     heterozygosity (He).
#' }
#' The results are compiled into a standardized list, including a frequency
#' matrix padded to accommodate the locus with the maximum number of alleles.
#'
#' @param genotypedata A data frame containing all raw allele observations from
#'   all samples (from both `late_failures` and `additional` data).
#' @param allele_definitionsA list where each element corresponds to a locus
#'   and contains the allele bin definitions for that locus, as created by
#'   the `define_alleles` function.
#' @param marker_info A data frame containing marker metadata, used to determine
#'   the appropriate calculation method for each locus.
#'
#' @return A list containing the following four elements:
#' \item{n_alleles}{A named integer vector with the total number of unique
#'   alleles identified for each locus.}
#' \item{freq_matrix}{A matrix where rows represent loci and columns represent
#'   alleles. Each cell contains the calculated frequency of an allele at a
#'   given locus.}
#' \item{variability}{A named numeric vector containing the calculated genetic
#'   variability (or error) for each locus.}
#' \item{allele_codes}{A named list that serves as a lookup table, mapping the
#'   integer-based column index of `freq_matrix` back to its true allele
#'   meaning (e.g., bin number or sequence).}
#'
#' @keywords internal
#' @noRd
#'
calculate_frequencies <- function(genotypedata, alleles_definitions, marker_info) {
  
  locus_names <- marker_info$marker_id
  n_loci <- length(locus_names)
  freq_list <- vector("list", n_loci)
  variability <- numeric(n_loci)
  n_alleles_per_locus <- integer(n_loci)
  allele_codes_list <- vector("list", n_loci)
  names(n_alleles_per_locus) <- locus_names
  names(variability) <- locus_names
  names(freq_list) <- locus_names
  names(allele_codes_list) <- locus_names
  
  for (locus_name in locus_names) {
    binning_method <- marker_info$binning_method[marker_info$marker_id == locus_name]
    locus_cols_pattern <- paste0("^", locus_name, "_")
    locus_cols <- grep(locus_cols_pattern, colnames(genotypedata), value = TRUE)
    
    if(length(locus_cols) == 0) { 
      raw_alleles <- numeric(0)
    } else {
      raw_alleles <- stats::na.omit(as.vector(as.matrix(genotypedata[, locus_cols, drop = FALSE])))
    }
    
    if (binning_method %in% c("microsatellite", "msp_glurp")) {
      current_definitions <- alleles_definitions[[locus_name]]
      
      num_bins <- if (!is.null(current_definitions)) nrow(current_definitions) else 0
      n_alleles_per_locus[locus_name] <- num_bins
      if (num_bins == 0 || length(raw_alleles) == 0) {
        freq_list[[locus_name]] <- if(num_bins > 0) rep(0, num_bins) else numeric(0)
        variability[locus_name] <- 0
        allele_codes_list[[locus_name]] <- if(num_bins > 0) 1:num_bins else integer(0) # <-- ADDED HERE
        next
      }
      recoded_bins <- findInterval(raw_alleles, current_definitions[, 1], rightmost.closed = TRUE)
      counts <- tabulate(recoded_bins[recoded_bins > 0], nbins = nrow(current_definitions))
      total_alleles <- sum(counts)
      freq_list[[locus_name]] <- if (total_alleles > 0) counts / total_alleles else rep(0, length(counts))
      sds_per_bin <- tapply(raw_alleles, recoded_bins, stats::sd, na.rm = TRUE)
      variability[locus_name] <- mean(sds_per_bin, na.rm = TRUE)
      if (is.nan(variability[locus_name])) { variability[locus_name] <- 0 }
      allele_codes_list[[locus_name]] <- 1:n_alleles_per_locus[locus_name]
      
    } else if (binning_method == "exact") {
      
      if (length(raw_alleles) == 0) {
        n_alleles_per_locus[locus_name] <- 0
        freq_list[[locus_name]] <- numeric(0)
        variability[locus_name] <- 0
        allele_codes_list[[locus_name]] <- numeric(0)
        next
      }
      allele_counts <- table(raw_alleles)
      total_alleles <- sum(allele_counts)
      current_freqs <- as.vector(allele_counts) / total_alleles
      freq_list[[locus_name]] <- current_freqs
      
      variability[locus_name] <- 1 - sum(current_freqs^2)
      n_alleles_per_locus[locus_name] <- length(allele_counts)
      allele_codes_list[[locus_name]] <- names(allele_counts)
    } else {
      stop("Unknown binning_method '", binning_method, "' for locus: ", locus_name)
    }
  }
  
  max_alleles <- if (length(n_alleles_per_locus) > 0) max(n_alleles_per_locus, na.rm = TRUE) else 0
  freq_matrix <- matrix(0, nrow = n_loci, ncol = max_alleles)
  if (n_loci > 0 && max_alleles > 0) {
    rownames(freq_matrix) <- locus_names
    colnames(freq_matrix) <- paste0("Allele_", 1:max_alleles)
    
    for (locus_name in locus_names) {
      num_alleles <- n_alleles_per_locus[locus_name]
      if (num_alleles > 0) {
        freq_matrix[locus_name, 1:num_alleles] <- freq_list[[locus_name]]
      }
    }
  }
  
  list(
    n_alleles = n_alleles_per_locus,
    freq_matrix = freq_matrix,
    variability = variability,
    allele_codes = allele_codes_list
  )
}

