
switch_hidden_ampseq <- function(x, hidden0, hiddenf, recoded0, recodedf,
                                 classification, qq, frequencies_RR, nloci, maxMOI) { 
  
  z <- stats::runif(1)
  hidden_positions <- which(c(hidden0[x,], hiddenf[x,]) == 1)
  
  if (length(hidden_positions) == 0) {
    return(list(recoded0=recoded0, recodedf=recodedf))
  }
  
  chosen <- sample(hidden_positions, 1)
  is_day0_allele <- chosen <= ncol(hidden0)
  
  if (is_day0_allele) {
    chosenlocus <- floor((chosen - 1) / maxMOI) + 1
    old_id <- recoded0[x, chosen]
  } else {
    chosenlocus <- floor((chosen - 1 - ncol(hidden0)) / maxMOI) + 1
    old_id <- recodedf[x, chosen - ncol(hidden0)]
  }
  
  possible_alleles <- frequencies_RR$allele_codes[[chosenlocus]]
  if (length(possible_alleles) == 0) {
    return(list(recoded0=recoded0, recodedf=recodedf))
  }
  
  new_id <- sample(possible_alleles, 1)
  
  calculate_log_lr <- function(patient_recoded0, patient_recodedf) {
    log_lr_total <- 0
    for (locus_idx in 1:nloci) {
      d0_alleles <- patient_recoded0[(maxMOI*(locus_idx-1)+1):(maxMOI*locus_idx)]
      df_alleles <- patient_recodedf[(maxMOI*(locus_idx-1)+1):(maxMOI*locus_idx)]
      d0_alleles <- d0_alleles[!is.na(d0_alleles)]
      df_alleles <- df_alleles[!is.na(df_alleles)]
      
      if (length(df_alleles) == 0) next
      
      is_recrud <- all(df_alleles %in% d0_alleles)
      log_prob_recrud <- log(if(is_recrud) 1 - qq else qq + 1e-10)
      
      freqs <- frequencies_RR$freq_matrix[locus_idx, 1:frequencies_RR$n_alleles[locus_idx]]
      allele_codes <- frequencies_RR$allele_codes[[locus_idx]]
      matched_indices <- match(df_alleles, allele_codes)
      
      if (any(is.na(matched_indices))) {
        log_prob_reinfect <- -1e6
      } else {
        log_prob_reinfect <- sum(log(freqs[matched_indices] + 1e-10))
      }
      log_lr_total <- log_lr_total + (log_prob_recrud - log_prob_reinfect)
    }
    return(log_lr_total)
  }
  
  
  log_lr_old <- calculate_log_lr(recoded0[x,], recodedf[x,])
  recoded0_new <- recoded0[x,]
  recodedf_new <- recodedf[x,]
  if (is_day0_allele) {
    recoded0_new[chosen] <- new_id
  } else {
    recodedf_new[chosen - ncol(hidden0)] <- new_id
  }
  
  log_lr_new <- calculate_log_lr(recoded0_new, recodedf_new)
  old_freq <- if(length(which(possible_alleles == old_id)) > 0) frequencies_RR$freq_matrix[chosenlocus, which(possible_alleles == old_id)] else 0
  new_freq <- if(length(which(possible_alleles == new_id)) > 0) frequencies_RR$freq_matrix[chosenlocus, which(possible_alleles == new_id)] else 0
  acceptance_ratio <- exp(log_lr_new - log_lr_old) * (new_freq / (old_freq + 1e-10))
  
  alpha <- min(1, acceptance_ratio)
  
  if (z < alpha) {
    # If accepted, update the actual state matrices
    if (is_day0_allele) {
      recoded0[x, chosen] <- new_id
    } else {
      recodedf[x, chosen - ncol(hidden0)] <- new_id
    }
  }
  return(list(recoded0=recoded0, recodedf=recodedf))
}