
calculate_frequencies3 <- function(genotypedata, alleles_definitions) {
  
  locus_names <- names(alleles_definitions)
  n_loci <- length(locus_names)
  freq_list <- vector("list", n_loci)
  variability <- numeric(n_loci)
  n_alleles_per_locus <- integer(n_loci)
  
  for (i in seq_along(locus_names)) {
    locus_name <- locus_names[i]
    current_definitions <- alleles_definitions[[locus_name]]
    
    if (is.null(current_definitions) || nrow(current_definitions) == 0) {
      n_alleles_per_locus[i] <- 0
      freq_list[[i]] <- numeric(0)
      variability[i] <- 0
      next
    }
    
    locus_cols_pattern <- paste0("^", locus_name, "_")
    locus_cols <- grep(locus_cols_pattern, colnames(genotypedata), value = TRUE)
    
    if(length(locus_cols) == 0) { 
      raw_alleles <- numeric(0)
    } else {
      raw_alleles <- stats::na.omit(as.vector(as.matrix(genotypedata[, locus_cols, drop = FALSE])))
    }
    
    if (length(raw_alleles) == 0) {
      n_alleles_per_locus[i] <- nrow(current_definitions)
      freq_list[[i]] <- rep(0, nrow(current_definitions))
      variability[i] <- 0
      next
    }
    
    recoded_bins <- findInterval(raw_alleles, current_definitions[, 1], rightmost.closed = TRUE)
    counts <- tabulate(recoded_bins, nbins = nrow(current_definitions))
    sds_per_bin <- tapply(raw_alleles, recoded_bins, stats::sd, na.rm = TRUE)
    variability[i] <- mean(sds_per_bin, na.rm = TRUE)
    if (is.nan(variability[i])) { variability[i] <- 0 }
    n_alleles_per_locus[i] <- length(counts)
    total_alleles <- sum(counts)
    freq_list[[i]] <- if (total_alleles > 0) counts / total_alleles else rep(0, length(counts))
  }
  
  max_alleles <- if (length(n_alleles_per_locus) > 0) max(n_alleles_per_locus) else 0
  freq_matrix <- matrix(0, nrow = n_loci, ncol = max_alleles)
  
  if (n_loci > 0) {
    rownames(freq_matrix) <- locus_names
    for (j in 1:n_loci) {
      if (n_alleles_per_locus[j] > 0) {
        freq_matrix[j, 1:n_alleles_per_locus[j]] <- freq_list[[j]]
      }
    }
  }
  
  list(
    n_alleles_per_locus,
    freq_matrix,
    variability
  )
}