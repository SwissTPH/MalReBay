define_alleles <- function(genotypedata, locirepeats, maxk) {
  
  locus_names <- names(locirepeats)
  n_loci <- length(locus_names)
  
  alleles_list <- vector("list", n_loci)
  names(alleles_list) <- locus_names
  
  # Loop over each locus by NAME, not by index
  for (locus_name in locus_names) {
    
    # Get the specific repeat length for this locus by name
    current_repeat_length <- locirepeats[[locus_name]]
    
    locus_cols_pattern <- paste0("^", locus_name, "_")
    locus_cols <- grep(locus_cols_pattern, colnames(genotypedata), value = TRUE)
    
    if (length(locus_cols) == 0) next
    
    raw_alleles <- stats::na.omit(as.vector(as.matrix(genotypedata[, locus_cols, drop = FALSE])))
    
    if (length(raw_alleles) == 0) {
      alleles_list[[locus_name]] <- matrix(NA, ncol = 3, nrow = 0, 
                                           dimnames = list(NULL, c("lower", "upper", "count")))
      next
    }
    
    if (diff(range(raw_alleles)) < current_repeat_length) {
      alleles_list[[locus_name]] <- matrix(
        c(min(raw_alleles) - current_repeat_length / 2, 
          max(raw_alleles) + current_repeat_length / 2, 
          length(raw_alleles)),
        nrow = 1,
        dimnames = list(NULL, c("lower", "upper", "count"))
      )
    } else {
      breaks <- seq(floor(min(raw_alleles)) - 0.5, ceiling(max(raw_alleles)) + 0.5, by = 1)
      allele_values <- (breaks[-1] + breaks[-length(breaks)]) / 2
      hist_alleles <- graphics::hist(raw_alleles, breaks = breaks, plot = FALSE)
      
      counts_by_offset <- sapply(1:current_repeat_length, function(x) {
        sum(hist_alleles$counts[seq(x, length(hist_alleles$counts), by = current_repeat_length)])
      })
      
      possible_alleles <- allele_values[seq(which.max(counts_by_offset), length(allele_values), by = current_repeat_length)]
      
      if (min(raw_alleles) <= (min(possible_alleles) - current_repeat_length / 2)) {
        possible_alleles <- c(min(possible_alleles) - current_repeat_length, possible_alleles)
      }
      if (max(raw_alleles) > (max(possible_alleles) + current_repeat_length / 2)) {
        possible_alleles <- c(possible_alleles, max(possible_alleles) + current_repeat_length)
      }
      
      clusters <- sapply(raw_alleles, function(x) which.min(abs(possible_alleles - x)))
      unique_clusters <- sort(unique(clusters))
      
      lower_bounds <- possible_alleles[unique_clusters] - current_repeat_length / 2
      upper_bounds <- possible_alleles[unique_clusters] + current_repeat_length / 2
      
      counts <- sapply(seq_along(lower_bounds), function(i) {
        sum(raw_alleles > lower_bounds[i] & raw_alleles <= upper_bounds[i])
      })
      
      alleles_list[[locus_name]] <- cbind(lower = lower_bounds, upper = upper_bounds, count = counts)
    }
  }
  
  final_alleles <- vector("list", n_loci)
  names(final_alleles) <- locus_names
  
  for (locus_name in locus_names) {
    locus_alleles <- alleles_list[[locus_name]]
    
    if (nrow(locus_alleles) == 0) {
      final_alleles[[locus_name]] <- matrix(NA, ncol = 2, nrow = 0)
      next
    }
    
    current_maxk <- ifelse(!is.null(names(maxk)), maxk[[locus_name]], maxk[1])
    num_alleles_to_keep <- min(current_maxk, nrow(locus_alleles))
    sorted_indices <- order(locus_alleles[, "count"], decreasing = TRUE)[1:num_alleles_to_keep]
    final_bins <- locus_alleles[sorted_indices, c("lower", "upper"), drop = FALSE]
    final_alleles[[locus_name]] <- final_bins[order(final_bins[, "lower"]), , drop = FALSE]
  }
  
  return(final_alleles)
}