
# Define the chain function for parallel processing
run_one_chain <- function(chain_id,
                          nruns, burnin, record_interval,
                          nids, ids, nloci, maxMOI, locinames,
                          genotypedata_RR, additional_neutral, alleles_definitions_RR
) {
  
  loglikelihood_chain <- rep(NA_real_, nruns)  # Track unset values
  attr(loglikelihood_chain, "bounds") <- c(-1e100, 1e100)  
  
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
  recrf <- matrix(NA, nids, nloci)
  recr_repeatsf <- matrix(NA, nids, nloci) 
  if (length(additional_neutral) > 0 && nrow(additional_neutral) > 0) {
    recoded_additional_neutral <- matrix(0, nrow(additional_neutral), maxMOI * nloci)
  }
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
    for (i in 1:ncolumns) { newalleles[, i] = sapply(1:nrow(oldalleles), function(x) recodeallele(alleles_definitions_RR[[j]],
                                                                                                  oldalleles[x, i])) }
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
      for (i in 1:ncolumns) { newalleles[, i] = sapply(1:nrow(oldalleles), function(x) recodeallele(alleles_definitions_RR[[j]],
                                                                                                    oldalleles[x, i])) }
      newalleles = matrix(as.numeric(newalleles), nrow = nrow(oldalleles), ncol = ncolumns)
      newalleles[is.na(newalleles)] = 0
      oldalleles = matrix(as.numeric(oldalleles), nrow = nrow(oldalleles), ncol = ncolumns)
      oldalleles[is.na(oldalleles)] = 0
      oldalleles[newalleles == 0] = 0
      recoded_additional_neutral[, (maxMOI*(j-1)+1):(maxMOI*(j-1)+ncolumns)] = newalleles
    }
  } else { recoded_additional_neutral = c() }
  
  
  frequencies_RR <- calculate_frequencies3(rbind(genotypedata_RR, additional_neutral), alleles_definitions_RR)
  
  
  ## assign random hidden alleles and classifications
  for (i in 1:nids) {
    for (j in 1:nloci) {
      nalleles0 = sum(alleles0[i,(maxMOI*(j-1)+1) : (maxMOI*(j))] != 0); nmissing0 = MOI0[i] - nalleles0; whichnotmissing0 =
        ((maxMOI*(j-1)+1) : (maxMOI*(j)))[which(alleles0[i,(maxMOI*(j-1)+1) : (maxMOI*(j-1)+MOI0[i])] != 0)]; whichmissing0 =
          ((maxMOI*(j-1)+1) : (maxMOI*(j)))[which(alleles0[i,(maxMOI*(j-1)+1) : (maxMOI*(j-1)+MOI0[i])] == 0)];
        if (nalleles0 > 0) { hidden0[i,whichnotmissing0] = 0 }
        if (nmissing0 > 0) { newhiddenalleles0 =
          sample(1:(frequencies_RR[[1]][j]),nmissing0,replace=TRUE,frequencies_RR[[2]][j,1:(frequencies_RR[[1]][j])]);
        recoded0[i,whichmissing0] = newhiddenalleles0; alleles0[i,whichmissing0] =
          rowMeans(alleles_definitions_RR[[j]])[newhiddenalleles0]; hidden0[i,whichmissing0] = 1; }
        nallelesf = sum(allelesf[i,(maxMOI*(j-1)+1) : (maxMOI*(j))] != 0); nmissingf = MOIf[i] - nallelesf; whichnotmissingf =
          ((maxMOI*(j-1)+1) : (maxMOI*(j)))[which(allelesf[i,(maxMOI*(j-1)+1) : (maxMOI*(j-1)+MOIf[i])] != 0)]; whichmissingf =
          ((maxMOI*(j-1)+1) : (maxMOI*(j)))[which(allelesf[i,(maxMOI*(j-1)+1) : (maxMOI*(j-1)+MOIf[i])] == 0)];
        if (nallelesf > 0) { hiddenf[i,whichnotmissingf] = 0 }
        if (nmissingf > 0) { newhiddenallelesf =
          sample(1:(frequencies_RR[[1]][j]),nmissingf,replace=TRUE,frequencies_RR[[2]][j,1:(frequencies_RR[[1]][j])]);
        recodedf[i,whichmissingf] = newhiddenallelesf; allelesf[i,whichmissingf] =
          rowMeans(alleles_definitions_RR[[j]])[newhiddenallelesf]; hiddenf[i,whichmissingf] = 1; }
    }
  }
  
  
  qq <- mean(c(hidden0, hiddenf), na.rm = TRUE)
  dvect <- stats::dgeom(0:(round(max(sapply(1:nloci, function(x) diff(range(c(alleles_definitions_RR[[x]]))))))), 0.75)
  
  ## randomly assign recrudescences/reinfections for this chain
  
  classification <- ifelse(stats::runif(nids) < 0.5, 1, 0)
  for (i in 1:nids) {
    for (j in 1:nloci) {
      allpossiblerecrud = expand.grid(1:MOI0[i],1:MOIf[i]); closestrecrud = which.min(sapply(1:dim(allpossiblerecrud)[1], function (x)
        abs(alleles0[i,maxMOI*(j-1)+allpossiblerecrud[x,1]] - allelesf[i,maxMOI*(j-1)+allpossiblerecrud[x,2]]))); mindistance[i,j] =
          abs(alleles0[i,maxMOI*(j-1)+allpossiblerecrud[closestrecrud,1]] - allelesf[i,maxMOI*(j-1)+allpossiblerecrud[closestrecrud,2]]);
        alldistance[i,j,1:dim(allpossiblerecrud)[1]] = sapply(1:dim(allpossiblerecrud)[1], function (x)
          abs(alleles0[i,maxMOI*(j-1)+allpossiblerecrud[x,1]] - allelesf[i,maxMOI*(j-1)+allpossiblerecrud[x,2]]));
        allrecrf[i,j,1:dim(allpossiblerecrud)[1]] = recodedf[i,maxMOI*(j-1)+allpossiblerecrud[,2]]; recr0[i,j] =
          maxMOI*(j-1)+allpossiblerecrud[closestrecrud,1]; recrf[i,j] = maxMOI*(j-1)+allpossiblerecrud[closestrecrud,2];
          recr_repeats0[i,j] = sum(recoded0[i,(maxMOI*(j-1)+1) : (maxMOI*(j))] == recoded0[i,recr0[i,j]])
          recr_repeatsf[i,j] = sum(recodedf[i,(maxMOI*(j-1)+1) : (maxMOI*(j))] == recodedf[i,recrf[i,j]])
    }
  }
  
  correction_distance_matrix <- list()
  for (i in 1:nloci) { correction_distance_matrix[[i]] <- as.matrix(stats::dist(rowMeans(alleles_definitions_RR[[i]]))) }
  
  
  
  
  state_classification <- matrix(NA, nids, (nruns - burnin) / record_interval)
  state_parameters <- matrix(NA, 2 + 2 * nloci, (nruns - burnin) / record_interval)
  state_alleles0 <- array(NA, c(nids, maxMOI * nloci, (nruns - burnin) / record_interval))
  state_allelesf <- array(NA, c(nids, maxMOI * nloci, (nruns - burnin) / record_interval))
  state_loglikelihood <- rep(NA_real_, (nruns - burnin) / record_interval)
  
  
  # Define ALL variables needed by runmcmc() before it's defined.
  count <- 1
  dposterior <- 0.75
  
  
  # Define MCMC function
  
  runmcmc <- function() {
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
      numerator <- dvect[distances[valid] + 1]
      
      denominators <- sapply(which(valid), function(z) {
        recr_allele <- allrecrf[x, y, z]
        if (is.na(recr_allele)) {
          return(NA) 
        }
        sum(
          frequencies_RR[[2]][y, 1:frequencies_RR[[1]][y]] *
            dvect[correction_distance_matrix[[y]][, recr_allele] + 1]
        )
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
    loglikelihood_chain[count] <<- ifelse(is.finite(loglik_val), loglik_val, NA)
    
    for (i in 1:nids) {
      updated_states <- switch_hidden(
        x = i,
        hidden0 = hidden0, hiddenf = hiddenf, recoded0 = recoded0, recodedf = recodedf, 
        alleles0 = alleles0, allelesf = allelesf, classification = classification,
        mindistance = mindistance, alldistance = alldistance, allrecrf = allrecrf,
        recr0 = recr0, recrf = recrf, recr_repeats0 = recr_repeats0, recr_repeatsf = recr_repeatsf,
        nloci = nloci, maxMOI = maxMOI, MOI0 = MOI0, MOIf = MOIf, qq = qq, dvect = dvect,
        alleles_definitions_RR = alleles_definitions_RR, frequencies_RR = frequencies_RR,
        correction_distance_matrix = correction_distance_matrix
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
    }
    
    q_prior_alpha = 0; 
    q_prior_beta = 0; 
    q_posterior_alpha = q_prior_alpha + sum(c(hidden0,hiddenf) == 1,na.rm=TRUE); 
    q_posterior_beta = q_prior_beta + sum(c(hidden0,hiddenf)==0,na.rm=TRUE); 
    if (q_posterior_alpha == 0) { q_posterior_alpha =1 }; 
    qq <<- stats::rbeta(1, q_posterior_alpha , q_posterior_beta)
    
    if (sum(classification==1) >= 1) {
      d_prior_alpha = 0; d_prior_beta = 0; 
      d_posterior_alpha = d_prior_alpha + length(c(mindistance[classification==1,])); 
      d_posterior_beta = d_prior_beta + sum(c(round(mindistance[classification==1,])));
      if (d_posterior_beta == 0) { d_posterior_beta = sum(c((mindistance[classification==1,]))) }
      if (d_posterior_beta == 0) { d_posterior_beta = 1 }
      dposterior <<- stats::rbeta(1, d_posterior_alpha , d_posterior_beta); 
      dvect = (1-dposterior) ^ (1:length(dvect)-1) * dposterior; 
      dvect <<- dvect / (sum(dvect))
    }
    
    
    tempdata <- recoded0
    recrudescence_indices <- which(classification == 1)
    if (length(recrudescence_indices) > 0) {
      for (idx in recrudescence_indices) { tempdata[idx, recr0[idx,]] <- 0 }
    }
    tempdata <- rbind(tempdata, recodedf)
    full_data <- rbind(tempdata, recoded_additional_neutral)
    
    for (locus_idx in 1:nloci) {
      updated_frequencies <- findposteriorfrequencies(
        locus_index = locus_idx,
        tempdata = full_data,
        maxMOI = maxMOI,
        frequencies_RR = frequencies_RR
      )
      frequencies_RR <<- updated_frequencies
    }
    
    
    if (count > burnin & count %% record_interval == 0) {
      record_idx <- (count - burnin) / record_interval
      state_classification[, record_idx] <<- classification
      state_alleles0[,, record_idx] <<- alleles0
      state_allelesf[,, record_idx] <<- allelesf
      state_parameters[1, record_idx] <<- qq
      state_parameters[2, record_idx] <<- dposterior
      state_parameters[3:(3+nloci-1), record_idx] <<- apply(frequencies_RR[[2]],1,max)
      state_parameters[(3+nloci):(3+2*nloci-1), record_idx] <<- sapply(1:nloci,function (x) sum(frequencies_RR[[2]][x,]^2))
      state_loglikelihood[record_idx] <<- loglikelihood_chain[count]
    }
    
    
    
    for (locus_idx in 1:nloci) {
      updated_frequencies <- findposteriorfrequencies(
        locus_index = locus_idx,
        tempdata = full_data,
        maxMOI = maxMOI,
        frequencies_RR = frequencies_RR
      )
      
      
      frequencies_RR <<- updated_frequencies
    }
    count <<- count + 1
  }
  
  # Execute the MCMC for the current chain - This remains unchanged
  replicate(nruns, runmcmc())
  
  
  return(
    list(
      classification = state_classification[, !is.na(colSums(state_classification)), drop = FALSE],
      parameters     = state_parameters[, !is.na(colSums(state_parameters)), drop = FALSE],
      alleles0       = state_alleles0[,, !is.na(colSums(state_parameters)), drop = FALSE],
      allelesf       = state_allelesf[,, !is.na(colSums(state_parameters)), drop = FALSE],
      loglikelihood  = state_loglikelihood
    )
  )
  
} 