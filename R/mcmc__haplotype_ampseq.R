run_one_chain_ampseq <- function(chain_id,
                                 nruns, burnin, record_interval,
                                 nids, ids, nloci, maxMOI, locinames,
                                 genotypedata_RR, additional_neutral,
                                 marker_info,
                                 record_hidden_alleles = FALSE)
{
  all_patient_ids_in_data <- gsub(" .*", "", genotypedata_RR$Sample.ID)
  ids <- unique(all_patient_ids_in_data)
  nids <- length(ids)
  MOI0 <- rep(0, nids)
  MOIf <- rep(0, nids)
  recoded0 <- matrix(NA_character_, nids, maxMOI * nloci)
  recodedf <- matrix(NA_character_, nids, maxMOI * nloci)
  hidden0 <- matrix(NA_integer_, nids, maxMOI * nloci)
  hiddenf <- matrix(NA_integer_, nids, maxMOI * nloci)

  for (i in 1:nids) {
    patient_id <- ids[i]
    row_idx_d0 <- which(genotypedata_RR$Sample.ID == paste(patient_id, "Day 0"))
    row_idx_df <- which(genotypedata_RR$Sample.ID == paste(patient_id, "Day Failure"))
    if (length(row_idx_d0) == 0 || length(row_idx_df) == 0) {
        next
    }
    row_idx_d0 <- row_idx_d0[1]
    row_idx_df <- row_idx_df[1]
    
    for (j in 1:nloci) {
        locus <- locinames[j]
        locicolumns <- grepl(paste0("^", locus, "_"), colnames(genotypedata_RR))
        
        # Extract alleles for this patient at D0
        d0_alleles <- as.character(genotypedata_RR[row_idx_d0, locicolumns])
        d0_alleles <- d0_alleles[!is.na(d0_alleles) & d0_alleles != "NA"]
        
        # Extract alleles for this patient at DF
        df_alleles <- as.character(genotypedata_RR[row_idx_df, locicolumns])
        df_alleles <- df_alleles[!is.na(df_alleles) & df_alleles != "NA"]
        
        # Update MOI for this patient
        MOI0[i] <- max(MOI0[i], length(d0_alleles))
        MOIf[i] <- max(MOIf[i], length(df_alleles))
        
        # Place the non-NA alleles into the correct row of the recoded matrices
        if (length(d0_alleles) > 0) {
            start_col <- maxMOI * (j - 1) + 1
            recoded0[i, start_col:(start_col + length(d0_alleles) - 1)] <- d0_alleles
        }
        if (length(df_alleles) > 0) {
            start_col <- maxMOI * (j - 1) + 1
            recodedf[i, start_col:(start_col + length(df_alleles) - 1)] <- df_alleles
        }
    }
}

  recoded_additional_neutral <- matrix(NA_character_, 0, maxMOI * nloci)
  if (length(additional_neutral) > 0 && nrow(additional_neutral) > 0) {
    # Initialize and fill additional neutral data here if needed, also using drop=FALSE
  }

  alleles_definitions_RR <- list() 
  frequencies_RR <- calculate_frequencies3(
    genotypedata = rbind(genotypedata_RR, additional_neutral),
    alleles_definitions = alleles_definitions_RR,
    marker_info = marker_info
  )

  # Assign random hidden alleles (initial imputation)
  for (i in 1:nids) {
    for (j in 1:nloci) {
      d0_cols <- (maxMOI*(j-1)+1):(maxMOI*j)
      n_observed_locus_d0 <- sum(!is.na(recoded0[i, d0_cols]))
      n_to_impute_d0 <- MOI0[i] - n_observed_locus_d0
      
      if (n_to_impute_d0 > 0) {
        actually_empty_indices_d0 <- d0_cols[is.na(recoded0[i, d0_cols])]
        num_to_actually_fill_d0 <- min(n_to_impute_d0, length(actually_empty_indices_d0))
        
        if (num_to_actually_fill_d0 > 0) {
          indices_to_fill_d0 <- actually_empty_indices_d0[1:num_to_actually_fill_d0]
          hidden0[i, indices_to_fill_d0] <- 1
          possible_alleles <- frequencies_RR$allele_codes[[j]]
          if (length(possible_alleles) > 0) {
             allele_freqs <- frequencies_RR$freq_matrix[j, 1:frequencies_RR$n_alleles[j]]
             new_ids <- sample(possible_alleles, num_to_actually_fill_d0, replace=TRUE, prob=allele_freqs)
             recoded0[i, indices_to_fill_d0] <- new_ids
          }
        }
      }
      
      hidden0[i, d0_cols[!is.na(recoded0[i, d0_cols])]] <- 0
      df_cols <- (maxMOI*(j-1)+1):(maxMOI*j)      
      n_observed_locus_df <- sum(!is.na(recodedf[i, df_cols]))
      n_to_impute_df <- MOIf[i] - n_observed_locus_df
      
      if (n_to_impute_df > 0) {
        actually_empty_indices_df <- df_cols[is.na(recodedf[i, df_cols])]
        num_to_actually_fill_df <- min(n_to_impute_df, length(actually_empty_indices_df))
        
        if (num_to_actually_fill_df > 0) {
          indices_to_fill_df <- actually_empty_indices_df[1:num_to_actually_fill_df]
          hiddenf[i, indices_to_fill_df] <- 1
          
          possible_alleles <- frequencies_RR$allele_codes[[j]]
          if (length(possible_alleles) > 0) {
             allele_freqs <- frequencies_RR$freq_matrix[j, 1:frequencies_RR$n_alleles[j]]
             new_ids <- sample(possible_alleles, num_to_actually_fill_df, replace=TRUE, prob=allele_freqs)
             recodedf[i, indices_to_fill_df] <- new_ids
          }
        }
      }
      hiddenf[i, df_cols[!is.na(recodedf[i, df_cols])]] <- 0
    }
  }

  qq <- mean(c(hidden0, hiddenf), na.rm = TRUE); if (is.na(qq)) qq <- 0.1 
  prob_recrud <- 0.5
  classification <- ifelse(stats::runif(nids) < prob_recrud, 1, 0)
  num_records <- floor((nruns - burnin) / record_interval)
  state_classification <- matrix(NA, nids, num_records)
  state_parameters <- matrix(NA, 2 + nloci, num_records)
  state_loglikelihood <- rep(NA_real_, num_records)

  if (record_hidden_alleles) {
    state_recoded0 <- array(NA, c(nids, maxMOI * nloci, num_records))
    state_recodedf <- array(NA, c(nids, maxMOI * nloci, num_records))
  } else {
    state_recoded0 <- NULL
    state_recodedf <- NULL
  } 
  
  runmcmc_step <- function(current_state) {
    # Unpack the current state
    classification <- current_state$classification
    recoded0 <- current_state$recoded0
    recodedf <- current_state$recodedf
    hidden0 <- current_state$hidden0
    hiddenf <- current_state$hiddenf
    qq <- current_state$qq
    frequencies_RR <- current_state$frequencies_RR 
    prob_recrud <- current_state$prob_recrud

    # a. Update classification
    likelihoodratio <- sapply(1:nids, function(patient_idx) {
        # Patient-level hard penalty check
        num_matching_loci <- sum(sapply(1:nloci, function(j) {
            d0 <- recoded0[patient_idx, (maxMOI*(j-1)+1):(maxMOI*j)]
            df <- recodedf[patient_idx, (maxMOI*(j-1)+1):(maxMOI*j)]
            d0_alleles <- d0[!is.na(d0)]
            df_alleles <- df[!is.na(df)]
            if (length(d0_alleles) > 0 && length(df_alleles) > 0) {
                return(any(df_alleles %in% d0_alleles))
            } else {
                return(FALSE) 
            }
        }))

        if (num_matching_loci == 0) {
            log_lr_total <- -1e6 
            log_prior_ratio <- log(prob_recrud / (1 - prob_recrud + 1e-10))
            return(exp(log_lr_total + log_prior_ratio))
        }

        # Locus-by-locus scoring
        log_lr_total <- 0
        for (locus_idx in 1:nloci) {
            d0_alleles <- recoded0[patient_idx, (maxMOI*(locus_idx-1)+1):(maxMOI*locus_idx)]
            df_alleles <- recodedf[patient_idx, (maxMOI*(locus_idx-1)+1):(maxMOI*locus_idx)]
            d0_alleles <- d0_alleles[!is.na(d0_alleles)]
            df_alleles <- df_alleles[!is.na(df_alleles)]
            
            if (length(df_alleles) == 0) next 

            has_data <- length(d0_alleles) > 0 && length(df_alleles) > 0
            is_mismatch <- all(!(df_alleles %in% d0_alleles))

            if (has_data && is_mismatch) {
                log_prob_recrud <- -1e6 
            } else {
                n_matched <- sum(df_alleles %in% d0_alleles)
                n_total <- length(unique(c(df_alleles, d0_alleles)))
                match_score <- if (n_total == 0) 1 else n_matched / n_total
                prob_recrud_locus <- (1 - qq) * match_score + qq * (1 - match_score)
                log_prob_recrud <- log(prob_recrud_locus + 1e-10)
            }
            
            matched_indices <- match(df_alleles, frequencies_RR$allele_codes[[locus_idx]])
            if (any(is.na(matched_indices))) { 
                log_prob_reinfect <- -1e6 
            } else { 
                freqs <- frequencies_RR$freq_matrix[locus_idx, 1:frequencies_RR$n_alleles[locus_idx]]
                log_prob_reinfect <- sum(log(freqs[matched_indices] + 1e-10)) 
            }
            log_lr_total <- log_lr_total + (log_prob_recrud - log_prob_reinfect)
        }
        log_prior_ratio <- log(prob_recrud / (1 - prob_recrud + 1e-10))
        ratio <- exp(log_lr_total + log_prior_ratio)
        return(pmin(pmax(ratio, 1e-100), 1e100))
    })
    
    z <- stats::runif(nids)
    newclassification <- classification
    newclassification[classification == 0 & z < likelihoodratio] <- 1
    newclassification[classification == 1 & z < 1/likelihoodratio] <- 0
    classification <- newclassification

    # b. Update hidden states
    for (i in 1:nids) {
        updated_states <- switch_hidden_ampseq(x = i, hidden0=hidden0, hiddenf=hiddenf, recoded0=recoded0, recodedf=recodedf,
                                               classification=classification, qq=qq, frequencies_RR=frequencies_RR, nloci=nloci, maxMOI=maxMOI)
        recoded0 <- updated_states$recoded0
        recodedf <- updated_states$recodedf
    }
    
    # c. Update q (error/dropout rate)
    q_posterior_alpha <- 1 + sum(c(hidden0, hiddenf) == 1, na.rm=TRUE)
    q_posterior_beta <- 1 + sum(c(hidden0, hiddenf) == 0, na.rm=TRUE)
    qq <- stats::rbeta(1, q_posterior_alpha, q_posterior_beta)
    
    # d. Update prior probability of recrudescence
    posterior_alpha_recrud <- 1 + sum(classification == 1)
    posterior_beta_recrud <- 1 + sum(classification == 0)
    prob_recrud <- stats::rbeta(1, posterior_alpha_recrud, posterior_beta_recrud)

    # e. Update population allele frequencies
    tempdata <- rbind(recoded0, recodedf, recoded_additional_neutral)
    new_freq_matrix <- frequencies_RR$freq_matrix
    for (locus_idx in 1:nloci) {
      new_row <- findposteriorfrequencies(locus_index=locus_idx, tempdata=tempdata, maxMOI=maxMOI, frequencies_RR=frequencies_RR)
      n_alleles_locus <- frequencies_RR$n_alleles[locus_idx]
      if (length(new_row) > 0 && !any(is.na(new_row)) && n_alleles_locus > 0) {
        new_freq_matrix[locus_idx, 1:n_alleles_locus] <- new_row
      }
    }
    frequencies_RR$freq_matrix <- new_freq_matrix
    
    # Return the updated state list
    return(list(classification = classification, recoded0 = recoded0, recodedf = recodedf, hidden0 = hidden0, hiddenf = hiddenf,
                qq = qq, prob_recrud = prob_recrud, frequencies_RR = frequencies_RR, likelihood_ratio = likelihoodratio))
  }
  
  mcmc_state <- list(classification = classification, recoded0 = recoded0, recodedf = recodedf, hidden0 = hidden0, hiddenf = hiddenf,
                     qq = qq, prob_recrud = prob_recrud, frequencies_RR = frequencies_RR)
                     
  for (iter in 1:nruns) {
    mcmc_state <- runmcmc_step(mcmc_state)
    if (iter > burnin & (iter - burnin) %% record_interval == 0) {
      record_idx <- (iter - burnin) / record_interval
      likelihoodratio_for_record <- mcmc_state$likelihood_ratio
      loglik_val <- sum(log(likelihoodratio_for_record[likelihoodratio_for_record > 0]), na.rm = TRUE)
      if (is.infinite(loglik_val)) loglik_val <- -1e6

      state_classification[, record_idx] <- mcmc_state$classification
      state_parameters[1, record_idx] <- mcmc_state$qq
      state_parameters[2, record_idx] <- mcmc_state$prob_recrud
      he_per_locus <- 1 - rowSums(mcmc_state$frequencies_RR$freq_matrix^2, na.rm=TRUE)
      state_parameters[3:(2+nloci), record_idx] <- he_per_locus
      state_loglikelihood[record_idx] <- loglik_val
    }
  }
  
  return(list(
    state_classification = state_classification,
    state_parameters = state_parameters,
    state_loglikelihood = state_loglikelihood,
    state_recoded0 = state_recoded0,
    state_recodedf = state_recodedf
  ))
}