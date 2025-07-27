
# Define the chain function for parallel processing
run_one_chain <- function(chain_id,
                          nruns, burnin, record_interval,
                          nids, ids, nloci, maxMOI, locinames,
                          genotypedata_RR, additional_neutral, alleles_definitions_RR,
                          marker_info,
                          record_hidden_alleles = FALSE
) 
{
  ##### calculate MOI
  MOI0 <- rep(0, nids)
  MOIf <- rep(0, nids)
  for (i in 1:nids) {
    for (j in 1:nloci) {
      locicolumns <- grepl(paste(locinames[j], "_", sep = ""), colnames(genotypedata_RR))
      nalleles0 <- sum(!is.na(genotypedata_RR[grepl(paste(ids[i], "Day 0"), genotypedata_RR$Sample.ID), locicolumns]))
      nallelesf <- sum(!is.na(genotypedata_RR[grepl(paste(ids[i], "Day Failure"), genotypedata_RR$Sample.ID), locicolumns]))
      MOI0[i] <- max(MOI0[i], nalleles0)
      MOIf[i] <- max(MOIf[i], nallelesf)
    }
  }
  
  ##### define state vector and create state 0
  
  alleles0 <- matrix(0, nids, maxMOI * nloci)
  recoded0 <- matrix(0, nids, maxMOI * nloci)
  hidden0 <- matrix(NA, nids, maxMOI * nloci)
  recr0 <- matrix(NA, nids, nloci)
  recr_repeats0 <- matrix(NA, nids, nloci)
  allelesf <- matrix(0, nids, maxMOI * nloci)
  recodedf <- matrix(0, nids, maxMOI * nloci)
  hiddenf <- matrix(NA, nids, maxMOI * nloci)
  hidden_crossfamily0 <- matrix(0, nids, maxMOI * nloci)
  hidden_crossfamilyf <- matrix(0, nids, maxMOI * nloci)
  recrf <- matrix(NA, nids, nloci)
  recr_repeatsf <- matrix(NA, nids, nloci) 
  if (length(additional_neutral) > 0 && nrow(additional_neutral) > 0) {
    recoded_additional_neutral <- matrix(0, nrow(additional_neutral), maxMOI * nloci)
  }
  recoded_additional_neutral <- matrix(0, nrow = 0, ncol = maxMOI * nloci)
  mindistance <- matrix(0, nids, nloci)
  alldistance <- array(NA, c(nids, nloci, maxMOI * maxMOI))
  allrecrf <- array(NA, c(nids, nloci, maxMOI * maxMOI))
  classification <- rep(0, nids)
  
  for (j in 1:nloci) {
    locus = locinames[j]
    locicolumns = grepl(paste0(locus, "_"), colnames(genotypedata_RR))
    oldalleles = as.matrix(genotypedata_RR[, locicolumns])
    if (is.null(dim(oldalleles))) { oldalleles = matrix(oldalleles, ncol=1) }
    ncolumns = ncol(oldalleles)
    newalleles = matrix(NA, nrow = nrow(oldalleles), ncol = ncolumns)
    
    for (i in 1:ncolumns) {
      temp_recode_col <- numeric(nrow(oldalleles))
      for (row_index in 1:nrow(oldalleles)) {
        # Assuming recodeallele is defined elsewhere
        recode_val <- recodeallele(alleles_definitions_RR[[j]], oldalleles[row_index, i])
        
        if (length(recode_val) != 1 || is.na(recode_val)) {
          temp_recode_col[row_index] <- NA 
        } else {
          temp_recode_col[row_index] <- recode_val
        }
      }
      newalleles[, i] <- temp_recode_col
    }
    
    newalleles = matrix(as.numeric(newalleles), nrow = nrow(oldalleles), ncol = ncolumns)
    newalleles[is.na(newalleles)] = 0
    oldalleles = matrix(as.numeric(oldalleles), nrow = nrow(oldalleles), ncol = ncolumns)
    oldalleles[is.na(oldalleles)] = 0
    oldalleles[newalleles == 0] = 0
    day0_rows = grepl("Day 0", genotypedata_RR$Sample.ID)
    dayf_rows = grepl("Day Failure", genotypedata_RR$Sample.ID)
    alleles0[, (maxMOI*(j-1)+1):(maxMOI*(j-1)+ncolumns)] = oldalleles[day0_rows, , drop=FALSE]
    allelesf[, (maxMOI*(j-1)+1):(maxMOI*(j-1)+ncolumns)] = oldalleles[dayf_rows, , drop=FALSE]
    recoded0[, (maxMOI*(j-1)+1):(maxMOI*(j-1)+ncolumns)] = newalleles[day0_rows, , drop=FALSE]
    recodedf[, (maxMOI*(j-1)+1):(maxMOI*(j-1)+ncolumns)] = newalleles[dayf_rows, , drop=FALSE]
  }
  
  if (length(additional_neutral) > 0 && nrow(additional_neutral) > 0) {
    recoded_additional_neutral = matrix(0, nrow = nrow(additional_neutral), ncol = maxMOI * nloci)
    for (j in 1:nloci) {
      locus = locinames[j]
      locicolumns = grepl(paste0(locus, "_"), colnames(additional_neutral))
      oldalleles = as.matrix(additional_neutral[, locicolumns])
      if (is.null(dim(oldalleles))) { oldalleles = matrix(oldalleles, ncol = 1) }
      ncolumns = ncol(oldalleles)
      newalleles = matrix(NA, nrow = nrow(oldalleles), ncol = ncolumns)
      
      # Robust recoding loop (replaces the second sapply)
      for (i in 1:ncolumns) {
        temp_recode_col <- numeric(nrow(oldalleles))
        for (row_index in 1:nrow(oldalleles)) {
          recode_val <- recodeallele(alleles_definitions_RR[[j]], oldalleles[row_index, i])
          
          if (length(recode_val) != 1 || is.na(recode_val)) {
            temp_recode_col[row_index] <- NA
          } else {
            temp_recode_col[row_index] <- recode_val
          }
        }
        newalleles[, i] <- temp_recode_col
      }
      
      newalleles = matrix(as.numeric(newalleles), nrow = nrow(oldalleles), ncol = ncolumns)
      newalleles[is.na(newalleles)] = 0
      oldalleles = matrix(as.numeric(oldalleles), nrow = nrow(oldalleles), ncol = ncolumns)
      oldalleles[is.na(oldalleles)] = 0
      oldalleles[newalleles == 0] = 0
      recoded_additional_neutral[, (maxMOI*(j-1)+1):(maxMOI*(j-1)+ncolumns)] = newalleles
    }
  } 
  
  frequencies_RR <- calculate_frequencies3(rbind(genotypedata_RR, additional_neutral), alleles_definitions_RR, marker_info)
  
  ## assign random hidden alleles and classifications
  for (i in 1:nids) {
    for (j in 1:nloci) {
      nalleles0 = sum(alleles0[i,(maxMOI*(j-1)+1) : (maxMOI*(j))] != 0); nmissing0 = MOI0[i] - nalleles0; whichnotmissing0 =
        ((maxMOI*(j-1)+1) : (maxMOI*(j)))[which(alleles0[i,(maxMOI*(j-1)+1) : (maxMOI*(j-1)+MOI0[i])] != 0)]; whichmissing0 =
          ((maxMOI*(j-1)+1) : (maxMOI*(j)))[which(alleles0[i,(maxMOI*(j-1)+1) : (maxMOI*(j-1)+MOI0[i])] == 0)];
        if (nalleles0 > 0) { hidden0[i,whichnotmissing0] = 0 }
        if (nmissing0 > 0) { 
          n_alleles_locus <- frequencies_RR$n_alleles[j]
          ### UPDATED ###
          allele_freqs = frequencies_RR$freq_matrix[j, 1:n_alleles_locus]
          if(is.na(n_alleles_locus) || n_alleles_locus == 0 || length(allele_freqs) == 0) next

          newhiddenalleles0 = sample(1:n_alleles_locus, nmissing0, replace=TRUE, prob=allele_freqs)
          recoded0[i,whichmissing0] = newhiddenalleles0
          alleles0[i,whichmissing0] = rowMeans(alleles_definitions_RR[[j]])[newhiddenalleles0]
          hidden0[i,whichmissing0] = 1 }

          nallelesf = sum(allelesf[i,(maxMOI*(j-1)+1) : (maxMOI*(j))] != 0); nmissingf = MOIf[i] - nallelesf; whichnotmissingf =
            ((maxMOI*(j-1)+1) : (maxMOI*(j)))[which(allelesf[i,(maxMOI*(j-1)+1) : (maxMOI*(j-1)+MOIf[i])] != 0)]; whichmissingf =
            ((maxMOI*(j-1)+1) : (maxMOI*(j)))[which(allelesf[i,(maxMOI*(j-1)+1) : (maxMOI*(j-1)+MOIf[i])] == 0)];
          if (nallelesf > 0) { hiddenf[i,whichnotmissingf] = 0 }
          if (nmissingf > 0) { 

            n_alleles_locus <- frequencies_RR$n_alleles[j]

            allele_freqs <- frequencies_RR$freq_matrix[j, 1:n_alleles_locus]
            if(is.na(n_alleles_locus) || n_alleles_locus == 0 || length(allele_freqs) == 0) next
            newhiddenallelesf <- sample(1:n_alleles_locus, nmissingf, replace=TRUE, prob=allele_freqs)
            recodedf[i,whichmissingf] <- newhiddenallelesf
            allelesf[i,whichmissingf] <- rowMeans(alleles_definitions_RR[[j]])[newhiddenallelesf]
            hiddenf[i,whichmissingf] <- 1
            }
      }
    }
  
  
  qq <- mean(c(hidden0, hiddenf), na.rm = TRUE)
  if (is.na(qq)) qq <- 0.1 
  dvect <- stats::dgeom(0:(round(max(sapply(1:nloci, function(x) diff(range(c(alleles_definitions_RR[[x]]))))))), 0.75)

  qq_crossfamily <- 10^-3 
  thresholddistance <- 15 

  getmode <- function(v) {
    v_clean <- v[!is.na(v)]
    if (length(v_clean) == 0) return(NA)
    uniqv <- unique(v_clean)
    uniqv[which.max(tabulate(match(v_clean, uniqv)))]
  }

  mode_allele_lengths = lapply(1:nloci, function (j) {
    sapply(1:nrow(alleles_definitions_RR[[j]]), function (y) {
      getmode(c(
        alleles0[,((j-1)*maxMOI+1):(j*maxMOI)][recoded0[,((j-1)*maxMOI+1):(j*maxMOI)]==y],
        allelesf[,((j-1)*maxMOI+1):(j*maxMOI)][recodedf[,((j-1)*maxMOI+1):(j*maxMOI)]==y]
      ))
    })
  })

  for(i in 1:nloci) {
    na_indices <- which(is.na(mode_allele_lengths[[i]]))
    if(length(na_indices) > 0) {
      mode_allele_lengths[[i]][na_indices] <- rowMeans(alleles_definitions_RR[[i]][na_indices, , drop=FALSE])
    }
  }
  ## randomly assign recrudescences/reinfections for this chain
  
  classification <- ifelse(stats::runif(nids) < 0.5, 1, 0)

  marker_info <- marker_info[match(locinames, marker_info$marker_id), ]
  
  for (i in 1:nids) {
    for (j in 1:nloci) {
      if (MOI0[i] > 0 && MOIf[i] > 0) {
        
        allpossiblerecrud <- expand.grid(1:MOI0[i], 1:MOIf[i])
        method <- marker_info$binning_method[j]

        if (method == "microsatellite" || method == "cluster") {
          distances <- sapply(1:nrow(allpossiblerecrud), function(x) {
            abs(alleles0[i, maxMOI * (j - 1) + allpossiblerecrud[x, 1]] - allelesf[i, maxMOI * (j - 1) + allpossiblerecrud[x, 2]])
          })
        } else { 
          distances <- sapply(1:nrow(allpossiblerecrud), function(x) {
            allele_day0_type <- recoded0[i, maxMOI * (j - 1) + allpossiblerecrud[x, 1]]
            allele_dayf_type <- recodedf[i, maxMOI * (j - 1) + allpossiblerecrud[x, 2]]
            if (is.na(allele_day0_type) || is.na(allele_dayf_type)) {
              return(Inf) 
            }
            return(ifelse(allele_day0_type == allele_dayf_type, 0, 1))
          })
        }

        if (all(is.na(distances))) {
          stop(paste("FATAL ERROR during initialization for Indiv", i, "Locus", j, ": All distances are NA. Halting."))
        }
        
        closestrecrud <- which.min(distances)
        mindistance[i, j] <- distances[closestrecrud]
        alldistance[i, j, 1:length(distances)] <- distances
        allrecrf[i, j, 1:nrow(allpossiblerecrud)] <- recodedf[i, maxMOI * (j - 1) + allpossiblerecrud[, 2]]
        recr0[i, j] <- maxMOI * (j - 1) + allpossiblerecrud[closestrecrud, 1]
        recrf[i, j] <- maxMOI * (j - 1) + allpossiblerecrud[closestrecrud, 2]
        recr_repeats0[i, j] <- sum(recoded0[i, (maxMOI * (j - 1) + 1):(maxMOI * j)] == recoded0[i, recr0[i, j]], na.rm = TRUE)
        recr_repeatsf[i, j] <- sum(recodedf[i, (maxMOI * (j - 1) + 1):(maxMOI * j)] == recodedf[i, recrf[i, j]], na.rm = TRUE)
      }
    }
  }
  correction_distance_matrix <- list()
  for (i in 1:nloci) { correction_distance_matrix[[i]] <- as.matrix(stats::dist(rowMeans(alleles_definitions_RR[[i]]))) }
  
  num_records <- floor((nruns - burnin) / record_interval)
  state_classification <- matrix(NA, nids, num_records)
  n_params_to_record <- 3 + (2 * nloci)
  state_parameters <- matrix(NA, n_params_to_record, num_records)
  state_loglikelihood <- rep(NA_real_, num_records)
  
  if (record_hidden_alleles) {
    state_alleles0 <- array(NA, c(nids, maxMOI * nloci, num_records))
    state_allelesf <- array(NA, c(nids, maxMOI * nloci, num_records))
  } else {
    state_alleles0 <- NULL
    state_allelesf <- NULL
  }
  
  # Define ALL variables needed by runmcmc() before it's defined.
  count <- 1
  dposterior <- 0.75

  # Define MCMC function
  runmcmc <- function() {
    current_marker_info <- marker_info 
    calculate_single_locus <- function(x, y) {
      should_print_debug <- (count < 5)
      distances <- round(alldistance[x, y, ])
      valid <- !is.na(distances) & distances >= 0 & distances < length(dvect)
      if (!any(valid)) {
        if (should_print_debug) {
          message(sprintf("\n--- DEBUG (Count %d, Indiv %d, Locus %d): EARLY EXIT ---", count, x, y))
          message("Reason: No valid distances found. All calculations for this locus will be skipped.")
          message(sprintf("Length of dvect: %d", length(dvect)))
          message("Distances found (includes NAs and out-of-bounds values):")
          print(utils::head(distances, 20)) 
        }
        return(1)
      }

      current_marker_type <- current_marker_info$markertype[y]
      method <- current_marker_info$binning_method[y]

      if (method == "microsatellite") {
          numerator <- dvect[distances[valid] + 1]
      } else if (method == "cluster") {
          threshold <- current_marker_info$cluster_gap_threshold[y]
          prob_match <- 1 - qq_crossfamily
          prob_mismatch <- qq_crossfamily
          numerator <- ifelse(distances[valid] <= threshold, prob_match, prob_mismatch)
      } else { 
          prob_match <- 1 - qq_crossfamily
          prob_mismatch <- qq_crossfamily
          numerator <- ifelse(distances[valid] == 0, prob_match, prob_mismatch)
      }
      
        denominators <- sapply(which(valid), function(z) {
        recr_allele <- allrecrf[x, y, z]
        if (is.na(recr_allele)) {
          return(NA) 
        }
      
      if (method == "microsatellite") {
          distances_to_pop <- correction_distance_matrix[[y]][, recr_allele]
          probs <- dvect[distances_to_pop + 1]
      } else if (method == "cluster") {
          threshold <- current_marker_info$cluster_gap_threshold[y]
          distances_to_pop <- correction_distance_matrix[[y]][, recr_allele]
          probs <- ifelse(distances_to_pop <= threshold, 1 - qq_crossfamily, qq_crossfamily)
      } else { # For "exact" markers
          all_allele_types <- 1:frequencies_RR$n_alleles[y]
          distances_to_pop <- ifelse(all_allele_types == recr_allele, 0, 1)
          probs <- ifelse(distances_to_pop == 0, 1 - qq_crossfamily, qq_crossfamily)
      }

      sum(frequencies_RR$freq_matrix[y, 1:frequencies_RR$n_alleles[y]] * probs, na.rm = TRUE)
    })
      epsilon <- 1e-10
      ratios <- numerator / (denominators + epsilon)
      
      final_ratios <- ratios[!is.na(ratios) & is.finite(ratios)]
      
      if (length(final_ratios) == 0) {
        if (should_print_debug) {
          message(sprintf("\n--- DEBUG (Count %d, Indiv %d, Locus %d): LATE EXIT ---", count, x, y))
          message("Reason: All calculated ratios were invalid (NA/Inf) and were filtered out.")
          message("Numerator(s) calculated:")
          print(utils::head(numerator, 20))
          message("Denominator(s) calculated (before adding epsilon):")
          print(utils::head(denominators, 20)) 
        }
        return(1) 
      }
      
      return(mean(final_ratios))
    }
    
    likelihoodratio <- sapply(1:nids, function(x) {
      loc_results <- sapply(1:nloci, function(y) calculate_single_locus(x, y))
      exp(sum(log(pmax(loc_results, 1e-10))))
    })
    z = stats::runif(nids); newclassification = classification; newclassification[classification == 0 & z < likelihoodratio] = 1; newclassification[classification == 1 & z < 1/likelihoodratio] = 0; classification <<- newclassification
    
    loglik_val <- sum(log(pmax(pmin(likelihoodratio, 1e100), 1e-100)))  # Clamp values
    # loglikelihood_chain[count] <<- ifelse(is.finite(loglik_val), loglik_val, NA)
    
    for (i in 1:nids) {
      updated_states <- switch_hidden_length(
        x = i,
        hidden0 = hidden0, hiddenf = hiddenf, recoded0 = recoded0, recodedf = recodedf, 
        alleles0 = alleles0, allelesf = allelesf, classification = classification,
        mindistance = mindistance, alldistance = alldistance, allrecrf = allrecrf,
        recr0 = recr0, recrf = recrf, recr_repeats0 = recr_repeats0, recr_repeatsf = recr_repeatsf,
        nloci = nloci, maxMOI = maxMOI, MOI0 = MOI0, MOIf = MOIf, qq = qq, dvect = dvect,
        alleles_definitions_RR = alleles_definitions_RR, frequencies_RR = frequencies_RR,
        correction_distance_matrix = correction_distance_matrix,
        marker_info = marker_info,
        mode_allele_lengths = mode_allele_lengths,
        thresholddistance = thresholddistance,
        qq_crossfamily = qq_crossfamily,
        hidden_crossfamily0 = hidden_crossfamily0,
        hidden_crossfamilyf = hidden_crossfamilyf
      )
      
      hidden0       <<- updated_states$hidden0
      hiddenf       <<- updated_states$hiddenf
      recoded0      <<- updated_states$recoded0
      recodedf      <<- updated_states$recodedf
      alleles0      <<- updated_states$alleles0
      allelesf      <<- updated_states$allelesf
      mindistance   <<- updated_states$mindistance
      alldistance   <<- updated_states$alldistance
      allrecrf      <<- updated_states$allrecrf
      recr0         <<- updated_states$recr0
      recrf         <<- updated_states$recrf
      recr_repeats0 <<- updated_states$recr_repeats0
      recr_repeatsf <<- updated_states$recr_repeatsf
      hidden_crossfamily0 <<- updated_states$hidden_crossfamily0
      hidden_crossfamilyf <<- updated_states$hidden_crossfamilyf
    }
    
    q_prior_alpha = 0; 
    q_prior_beta = 0; 
    q_posterior_alpha = q_prior_alpha + sum(c(hidden0,hiddenf) == 1,na.rm=TRUE); 
    q_posterior_beta = q_prior_beta + sum(c(hidden0,hiddenf)==0,na.rm=TRUE); 
    if (q_posterior_alpha == 0) { q_posterior_alpha =1 }; 
    qq <<- stats::rbeta(1, q_posterior_alpha , q_posterior_beta)

    q_cross_prior_alpha <- 1
    q_cross_prior_beta <- 999 # Prior belief: cross-family jumps are rare
    q_cross_posterior_alpha <- q_cross_prior_alpha + sum(c(hidden_crossfamily0, hidden_crossfamilyf) == 1, na.rm=TRUE)
    q_cross_posterior_beta <- q_cross_prior_beta + sum(c(hidden_crossfamily0, hidden_crossfamilyf) == 0, na.rm=TRUE)
    qq_crossfamily <<- stats::rbeta(1, q_cross_posterior_alpha, q_cross_posterior_beta)
    
    if (sum(classification == 1) >= 1) {
      d_prior_alpha <- 2
      d_prior_beta <- 5 
      valid_distances <- mindistance[classification == 1, ]
      valid_distances <- valid_distances[!is.na(valid_distances)]
      
      d_posterior_alpha <- d_prior_alpha + length(valid_distances)
      d_posterior_beta <- d_prior_beta + sum(valid_distances)
      if (d_posterior_beta <= 0) { d_posterior_beta <- 1 }
      dposterior <<- stats::rbeta(1, d_posterior_alpha, d_posterior_beta)
      # if (chain_id == 1 && count < 10) {
        # message(sprintf("Chain 1, Count %d: dposterior = %.4f", count, dposterior))
      # }
      dvect_new <- (1 - dposterior) ^ (seq_along(dvect) - 1) * dposterior
      dvect <<- dvect_new / sum(dvect_new)
    }
    
    tempdata <- recoded0
    recrudescence_indices <- which(classification == 1)
    if (length(recrudescence_indices) > 0) {
      for (idx in recrudescence_indices) { tempdata[idx, recr0[idx,]] <- 0 }
    }
    tempdata <- rbind(tempdata, recodedf)
    full_data <- rbind(tempdata, recoded_additional_neutral)
    
    new_freq_matrix <- frequencies_RR$freq_matrix
    
    for (locus_idx in 1:nloci) {
      new_row <- findposteriorfrequencies(
        locus_index = locus_idx, 
        tempdata = full_data, 
        maxMOI = maxMOI, 
        frequencies_RR = frequencies_RR
      )
      n_alleles_locus <- frequencies_RR$n_alleles[locus_idx]
      if (length(new_row) > 0 && !any(is.na(new_row)) && n_alleles_locus > 0) {
        new_freq_matrix[locus_idx, 1:n_alleles_locus] <- new_row
      }
    }
    
    
    if (count > burnin & count %% record_interval == 0) {
      record_idx <- (count - burnin) / record_interval
      loglik_val <- sum(log(pmax(pmin(likelihoodratio, 1e100), 1e-100)))
      state_classification[, record_idx] <<- classification
      state_parameters[1, record_idx] <<- qq
      state_parameters[2, record_idx] <<- qq_crossfamily
      state_parameters[3, record_idx] <<- dposterior
      state_parameters[4:(4+nloci-1), record_idx] <<- apply(frequencies_RR$freq_matrix, 1, max)
      state_parameters[(4+nloci):(4+2*nloci-1), record_idx] <<- sapply(1:nloci, function (x) sum(frequencies_RR$freq_matrix[x, 1:frequencies_RR$n_alleles[x]]^2))
      state_loglikelihood[record_idx] <<- ifelse(is.finite(loglik_val), loglik_val, NA)
      if (record_hidden_alleles) {
        state_alleles0[,, record_idx] <<- alleles0
        state_allelesf[,, record_idx] <<- allelesf
      }
    }
    count <<- count + 1
  }
  
  replicate(nruns, runmcmc())
  
  return(
    list(
      classification = state_classification,
      parameters     = state_parameters,
      alleles0       = state_alleles0,
      allelesf       = state_allelesf,
      loglikelihood  = state_loglikelihood
    )
  )
  
} 