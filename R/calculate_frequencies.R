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
#'   \item For **length polymorphism** data (`microsatellite`, `cluster`), it
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
#' @param alleles_definitions A list where each element corresponds to a locus
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
calculate_frequencies3 <- function(genotypedata, alleles_definitions, marker_info) {
  
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
    
    if (binning_method %in% c("microsatellite", "cluster")) {
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