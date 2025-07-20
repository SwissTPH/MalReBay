# Function to find posterior frequencies for a given locus index

findposteriorfrequencies =  function(locus_index, tempdata, maxMOI, frequencies_RR) {
  data_cols <- (1:maxMOI) + (locus_index - 1) * maxMOI
  locus_data <- tempdata[, data_cols]
  n_alleles <- frequencies_RR$n_alleles[locus_index]
  
  if (is.na(n_alleles) || n_alleles == 0) {
    max_alleles_in_matrix <- ncol(frequencies_RR$freq_matrix)
    return(rep(NA, max_alleles_in_matrix))
  }
  allele_levels <- frequencies_RR$allele_codes[[locus_index]]
  freq_prior_alpha <- rep(1, n_alleles)
  observed_counts <- table(factor(c(locus_data), levels = allele_levels))
  freq_posterior_alpha <- freq_prior_alpha + observed_counts
  new_frequencies <- gtools::rdirichlet(1, freq_posterior_alpha)
  
  # frequencies_RR$freq_matrix[locus_index, 1:n_alleles] <- new_frequencies
  return(as.vector(new_frequencies))
}
