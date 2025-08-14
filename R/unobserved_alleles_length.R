# Unobserved alleles 

switch_hidden_length = function(x, hidden0, hiddenf, recoded0, recodedf, alleles0, allelesf,
                                classification, mindistance, alldistance, allrecrf,
                                recr0, recrf, recr_repeats0, recr_repeatsf,
                                nloci, maxMOI, MOI0, MOIf, qq, dvect,
                                alleles_definitions_RR, frequencies_RR, correction_distance_matrix,
                                marker_info, mode_allele_lengths, thresholddistance, qq_crossfamily,
                                hidden_crossfamily0, hidden_crossfamilyf, is_locus_comparable) {
  
  if (MOI0[x] == 0 || MOIf[x] == 0) {
    return(list(
      hidden0=hidden0, hiddenf=hiddenf, recoded0=recoded0, recodedf=recodedf, 
      alleles0=alleles0, allelesf=allelesf, mindistance=mindistance, alldistance=alldistance, 
      allrecrf=allrecrf, recr0=recr0, recrf=recrf, 
      recr_repeats0=recr_repeats0, recr_repeatsf=recr_repeatsf,
      hidden_crossfamily0=hidden_crossfamily0, hidden_crossfamilyf=hidden_crossfamilyf
    ))
  }
  
  z = stats::runif(1)
  if (sum(hidden0[x,], hiddenf[x,],na.rm=TRUE) > 0) { # if hidden alleles exist 
    if (length(which(c(hidden0[x,], hiddenf[x,])==1)) > 1) {
      chosen <- sample(which(c(hidden0[x,], hiddenf[x,])==1), 1)
    } else {
      chosen <- which(c(hidden0[x,], hiddenf[x,])==1)
    }
    # ================================== REINFECTION ==================================
    if (classification[x] == 0) { # reinfection
      if (chosen <= nloci*maxMOI) { # day 0 hidden allele
        chosenlocus = ceiling(chosen/maxMOI)
        old = recoded0[x,chosen]
        n_alleles_locus <- frequencies_RR$n_alleles[chosenlocus]
        if(is.na(n_alleles_locus) || n_alleles_locus == 0) 
        return(list(hidden0=hidden0, hiddenf=hiddenf, recoded0=recoded0, recodedf=recodedf, 
                    alleles0=alleles0, allelesf=allelesf, mindistance=mindistance, alldistance=alldistance, 
                    allrecrf=allrecrf, recr0=recr0, recrf=recrf, 
                    recr_repeats0=recr_repeats0, recr_repeatsf=recr_repeatsf,
                    hidden_crossfamily0=hidden_crossfamily0, hidden_crossfamilyf=hidden_crossfamilyf))
        new = sample(1:n_alleles_locus, 1)
        oldalleles = recoded0[x,intersect(((chosenlocus-1)*maxMOI+1):((chosenlocus)*maxMOI),which(hidden0[x,] == 0))]
        repeatedold = qq
        repeatednew = qq
        
        if (sum(oldalleles == old) >= 1) { # if old allele is a repeat, don't penalize with missing probability
          repeatedold = 1;
        }
        if (sum(oldalleles == new) >= 1) { # if new allele is a repeat, don't penalize with missing probability
          repeatednew = 1;
        }
        
        crossfamily_new_flag <- 0
        penalty_numerator <- 1
        penalty_denominator <- qq_crossfamily^hidden_crossfamily0[x, chosen] 

        method <- marker_info$binning_method[chosenlocus]
        if (method == "cluster") {
          current_threshold <- marker_info$cluster_gap_threshold[chosenlocus]
          obs_lengths <- allelesf[x, which(hiddenf[x,] == 0 & floor((0:(ncol(allelesf)-1))/maxMOI)+1 == chosenlocus)]
          if (sum(!is.na(obs_lengths)) > 0) {
            if (min(abs(mode_allele_lengths[[chosenlocus]][new] - obs_lengths), na.rm=TRUE) > current_threshold) {
              penalty_numerator <- qq_crossfamily
              crossfamily_new_flag <- 1
            }
          }
        }
        
        method <- marker_info$binning_method[chosenlocus]
        
        if (method == "microsatellite") {
            numerator <- (sum(frequencies_RR$freq_matrix[chosenlocus, 1:n_alleles_locus]*dvect[correction_distance_matrix[[chosenlocus]][,new]+1], na.rm=TRUE) * repeatednew) * penalty_numerator
            denominator <- (sum(frequencies_RR$freq_matrix[chosenlocus, 1:n_alleles_locus]*dvect[correction_distance_matrix[[chosenlocus]][,old]+1], na.rm=TRUE) * repeatedold) * penalty_denominator
        } else if (method == "cluster") {
            threshold <- marker_info$cluster_gap_threshold[chosenlocus]
            probs_new <- ifelse(correction_distance_matrix[[chosenlocus]][,new] <= threshold, 1 - qq_crossfamily, qq_crossfamily)
            numerator <- (sum(frequencies_RR$freq_matrix[chosenlocus, 1:n_alleles_locus] * probs_new, na.rm=TRUE) * repeatednew) * penalty_numerator
            probs_old <- ifelse(correction_distance_matrix[[chosenlocus]][,old] <= threshold, 1 - qq_crossfamily, qq_crossfamily)
            denominator <- (sum(frequencies_RR$freq_matrix[chosenlocus, 1:n_alleles_locus] * probs_old, na.rm=TRUE) * repeatedold) * penalty_denominator
        } else { 
            probs_new <- ifelse(correction_distance_matrix[[chosenlocus]][,new] == 0, 1 - qq_crossfamily, qq_crossfamily)
            numerator <- (sum(frequencies_RR$freq_matrix[chosenlocus, 1:n_alleles_locus] * probs_new, na.rm=TRUE) * repeatednew) * penalty_numerator
            probs_old <- ifelse(correction_distance_matrix[[chosenlocus]][,old] == 0, 1 - qq_crossfamily, qq_crossfamily)
            denominator <- (sum(frequencies_RR$freq_matrix[chosenlocus, 1:n_alleles_locus] * probs_old, na.rm=TRUE) * repeatedold) * penalty_denominator
        }

        alpha <- 0
        if (is.na(numerator) || is.na(denominator)) {
          alpha <- 0
        } else if (denominator == 0 && numerator > 0) {
          alpha <- 1.00001 
        } else if (denominator > 0) {
          alpha <- numerator / denominator
        }
        
        if (z < alpha) { # switch made
          recoded0[x,chosen] = new
          sd_val <- frequencies_RR$variability[chosenlocus]
          newallele_length = mean(alleles_definitions_RR[[chosenlocus]][new,]) + stats::rnorm(1,mean=0,sd=sd_val)
          hidden_crossfamily0[x, chosen] <- crossfamily_new_flag
          alleles0[x,chosen] = newallele_length        
          allpossiblerecrud = expand.grid(1:MOI0[x],1:MOIf[x])
          closestrecrud = which.min(sapply(1:dim(allpossiblerecrud)[1], function (y) abs(alleles0[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,1]] - allelesf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,2]])))
          mindistance[x,chosenlocus] = abs(alleles0[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[closestrecrud,1]] - allelesf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[closestrecrud,2]])
          alldistance[x,chosenlocus,1:dim(allpossiblerecrud)[1]] = sapply(1:dim(allpossiblerecrud)[1], function (y) abs(alleles0[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,1]] - allelesf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,2]]))
          allrecrf[x,chosenlocus,1:dim(allpossiblerecrud)[1]] = recodedf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[,2]]
          recr0[x,chosenlocus] = maxMOI*(chosenlocus-1)+allpossiblerecrud[closestrecrud,1]
          recrf[x,chosenlocus] = maxMOI*(chosenlocus-1)+allpossiblerecrud[closestrecrud,2]
          recr_repeats0[x,chosenlocus] = sum(recoded0[x,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == recoded0[x,recr0[x,chosenlocus]])
          recr_repeatsf[x,chosenlocus] = sum(recodedf[x,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == recodedf[x,recrf[x,chosenlocus]])
          
        }
      } else { # day f hidden allele
        chosen = chosen - nloci*maxMOI
        chosenlocus = ceiling(chosen/maxMOI)
        old = recodedf[x,chosen]
        n_alleles_locus <- frequencies_RR$n_alleles[chosenlocus]
        if(is.na(n_alleles_locus) || n_alleles_locus == 0) return(list(hidden0=hidden0, hiddenf=hiddenf, recoded0=recoded0, recodedf=recodedf, 
                                                                       alleles0=alleles0, allelesf=allelesf, mindistance=mindistance, alldistance=alldistance, 
                                                                       allrecrf=allrecrf, recr0=recr0, recrf=recrf, 
                                                                       recr_repeats0=recr_repeats0, recr_repeatsf=recr_repeatsf,
                                                                       hidden_crossfamily0=hidden_crossfamily0, hidden_crossfamilyf=hidden_crossfamilyf))
        new = sample(1:n_alleles_locus, 1)
        oldalleles = recodedf[x,intersect(((chosenlocus-1)*maxMOI+1):((chosenlocus)*maxMOI),which(hiddenf[x,] == 0))]
        repeatedold = qq
        repeatednew = qq
        if (sum(oldalleles == old) >= 1) { # if old allele is a repeat, don't penalize with missing probability
          repeatedold = 1;
        }
        if (sum(oldalleles == new) >= 1) { # if new allele is a repeat, don't penalize with missing probability
          repeatednew = 1;
        }
        
        crossfamily_new_flag <- 0
        penalty_numerator <- 1
        penalty_denominator <- qq_crossfamily^hidden_crossfamilyf[x, chosen]

        method <- marker_info$binning_method[chosenlocus]
        if (method == "cluster") {
          current_threshold <- marker_info$cluster_gap_threshold[chosenlocus]
          obs_lengths <- allelesf[x, which(hiddenf[x,] == 0 & floor((0:(ncol(allelesf)-1))/maxMOI)+1 == chosenlocus)]
          if (sum(!is.na(obs_lengths)) > 0) {
            if (min(abs(mode_allele_lengths[[chosenlocus]][new] - obs_lengths), na.rm=TRUE) > current_threshold) {
              penalty_numerator <- qq_crossfamily
              crossfamily_new_flag <- 1
            }
          }
        }
        
        method <- marker_info$binning_method[chosenlocus]
        
        if (method == "microsatellite") {
            numerator <- (sum(frequencies_RR$freq_matrix[chosenlocus, 1:n_alleles_locus]*dvect[correction_distance_matrix[[chosenlocus]][,new]+1], na.rm=TRUE) * repeatednew) * penalty_numerator
            denominator <- (sum(frequencies_RR$freq_matrix[chosenlocus, 1:n_alleles_locus]*dvect[correction_distance_matrix[[chosenlocus]][,old]+1], na.rm=TRUE) * repeatedold) * penalty_denominator
        } else if (method == "cluster") {
            threshold <- marker_info$cluster_gap_threshold[chosenlocus]
            probs_new <- ifelse(correction_distance_matrix[[chosenlocus]][,new] <= threshold, 1 - qq_crossfamily, qq_crossfamily)
            numerator <- (sum(frequencies_RR$freq_matrix[chosenlocus, 1:n_alleles_locus] * probs_new, na.rm=TRUE) * repeatednew) * penalty_numerator
            probs_old <- ifelse(correction_distance_matrix[[chosenlocus]][,old] <= threshold, 1 - qq_crossfamily, qq_crossfamily)
            denominator <- (sum(frequencies_RR$freq_matrix[chosenlocus, 1:n_alleles_locus] * probs_old, na.rm=TRUE) * repeatedold) * penalty_denominator
        } else { 
            probs_new <- ifelse(correction_distance_matrix[[chosenlocus]][,new] == 0, 1 - qq_crossfamily, qq_crossfamily)
            numerator <- (sum(frequencies_RR$freq_matrix[chosenlocus, 1:n_alleles_locus] * probs_new, na.rm=TRUE) * repeatednew) * penalty_numerator
            probs_old <- ifelse(correction_distance_matrix[[chosenlocus]][,old] == 0, 1 - qq_crossfamily, qq_crossfamily)
            denominator <- (sum(frequencies_RR$freq_matrix[chosenlocus, 1:n_alleles_locus] * probs_old, na.rm=TRUE) * repeatedold) * penalty_denominator
        }

        alpha <- 0
        if (is.na(numerator) || is.na(denominator)) {
          alpha <- 0
        } else if (denominator == 0 && numerator > 0) {
          alpha <- 1.00001 
        } else if (denominator > 0) {
          alpha <- numerator / denominator
        }
        
        
        if (z < alpha) { # switch made
          recodedf[x,chosen] = new
          sd_val <- frequencies_RR$variability[chosenlocus]
          newallele_length = mean(alleles_definitions_RR[[chosenlocus]][new,]) + stats::rnorm(1,mean=0,sd=sd_val)
          hidden_crossfamily0[x, chosen] <- crossfamily_new_flag
          allelesf[x,chosen] = newallele_length
          allpossiblerecrud = expand.grid(1:MOI0[x],1:MOIf[x])
          closestrecrud = which.min(sapply(1:dim(allpossiblerecrud)[1], function (y) abs(alleles0[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,1]] - allelesf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,2]])))
          mindistance[x,chosenlocus] = abs(alleles0[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[closestrecrud,1]] - allelesf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[closestrecrud,2]])
          alldistance[x,chosenlocus,1:dim(allpossiblerecrud)[1]] = sapply(1:dim(allpossiblerecrud)[1], function (y) abs(alleles0[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,1]] - allelesf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,2]]))
          allrecrf[x,chosenlocus,1:dim(allpossiblerecrud)[1]] = recodedf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[,2]]
          recr0[x,chosenlocus] = maxMOI*(chosenlocus-1)+allpossiblerecrud[closestrecrud,1]
          recrf[x,chosenlocus] = maxMOI*(chosenlocus-1)+allpossiblerecrud[closestrecrud,2]
          recr_repeats0[x,chosenlocus] = sum(recoded0[x,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == recoded0[x,recr0[x,chosenlocus]])
          recr_repeatsf[x,chosenlocus] = sum(recodedf[x,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == recodedf[x,recrf[x,chosenlocus]])
          
        }
      }
      
      # ================================== RECRUDESCENCE ==================================
    } else { # recrudescence
      if (chosen <= nloci*maxMOI) { # day 0 hidden allele
        chosenlocus = ceiling(chosen/maxMOI)
        
        if (is_locus_comparable[x, chosenlocus]) {
          old = recoded0[x,chosen]
          n_alleles_locus <- frequencies_RR$n_alleles[chosenlocus]
          if(is.na(n_alleles_locus) || n_alleles_locus == 0) 
          return(list(hidden0=hidden0, hiddenf=hiddenf, recoded0=recoded0, recodedf=recodedf, alleles0=alleles0, 
                      allelesf=allelesf, mindistance=mindistance, alldistance=alldistance, allrecrf=allrecrf, 
                      recr0=recr0, recrf=recrf, recr_repeats0=recr_repeats0, recr_repeatsf=recr_repeatsf,
                      hidden_crossfamily0=hidden_crossfamily0, hidden_crossfamilyf=hidden_crossfamilyf))
          new = sample(1:n_alleles_locus, 1)
          sd_val <- frequencies_RR$variability[chosenlocus]
          newallele_length = mean(alleles_definitions_RR[[chosenlocus]][new,]) + stats::rnorm(1,mean=0,sd=sd_val)
          oldalleles = recoded0[x,intersect(((chosenlocus-1)*maxMOI+1):((chosenlocus)*maxMOI),which(hidden0[x,] == 0))]
          repeatedold = qq; repeatednew = qq
          if (sum(oldalleles == old) >= 1) { repeatedold = 1 }
          if (sum(oldalleles == new) >= 1) { repeatednew = 1 }
          
          crossfamily_new_flag <- 0; penalty_numerator <- 1
          penalty_denominator <- qq_crossfamily^hidden_crossfamily0[x, chosen]
          
          allpossiblerecrud = expand.grid(1:MOI0[x],1:MOIf[x])
          tempalleles = alleles0[x,maxMOI*(chosenlocus-1)+1:maxMOI]; tempalleles[chosen-(chosenlocus-1)*maxMOI] = newallele_length 
          temprecoded = recoded0[x,maxMOI*(chosenlocus-1)+1:maxMOI]; temprecoded[chosen-(chosenlocus-1)*maxMOI] = new
          
          newclosestrecrud = which.min(sapply(1:dim(allpossiblerecrud)[1], function (y) abs(tempalleles[allpossiblerecrud[y,1]] - allelesf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,2]])))
          newmindistance = abs(tempalleles[allpossiblerecrud[newclosestrecrud,1]] - allelesf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[newclosestrecrud,2]])
          newalldistance = sapply(1:dim(allpossiblerecrud)[1], function (y) abs(tempalleles[allpossiblerecrud[y,1]] - allelesf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,2]]))
          newallrecrf = recodedf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[,2]]
          epsilon <- 1e-9; method <- marker_info$binning_method[chosenlocus]
          calculate_likelihood <- function(distances_vec, recrf_vec, freq_mat, correction_mat, method, threshold) {
              if (method == "microsatellite") { num_probs <- dvect[round(distances_vec) + 1] } else if (method == "cluster") { num_probs <- ifelse(distances_vec <= threshold, 1 - qq_crossfamily, qq_crossfamily) } else { num_probs <- ifelse(distances_vec == 0, 1 - qq_crossfamily, qq_crossfamily) }
              den_probs <- sapply(recrf_vec, function(recr_allele) { if (is.na(recr_allele)) return(NA); if (method == "microsatellite") { pop_probs <- dvect[correction_mat[, recr_allele] + 1] } else if (method == "cluster") { pop_probs <- ifelse(correction_mat[, recr_allele] <= threshold, 1 - qq_crossfamily, qq_crossfamily) } else { pop_probs <- ifelse(correction_mat[, recr_allele] == 0, 1 - qq_crossfamily, qq_crossfamily) }; sum(freq_mat * pop_probs, na.rm = TRUE) })
              mean(num_probs / (den_probs + epsilon), na.rm = TRUE)
          }
          threshold <- if(method == "cluster") marker_info$cluster_gap_threshold[chosenlocus] else NA
          likelihoodnew <- calculate_likelihood(newalldistance, newallrecrf, frequencies_RR$freq_matrix[chosenlocus, 1:n_alleles_locus], correction_distance_matrix[[chosenlocus]], method, threshold) * repeatednew * penalty_numerator
          likelihoodold <- calculate_likelihood(alldistance[x, chosenlocus, ], allrecrf[x, chosenlocus, ], frequencies_RR$freq_matrix[chosenlocus, 1:n_alleles_locus], correction_distance_matrix[[chosenlocus]], method, threshold) * repeatedold * penalty_denominator
          
          alpha <- 0
          if (!is.na(likelihoodold) && is.finite(likelihoodold) && likelihoodold > 0) { alpha <- likelihoodnew / likelihoodold } else if (!is.na(likelihoodnew) && likelihoodnew > 0) { alpha <- 2.0 }
          if (is.na(alpha) || !is.finite(alpha)) { alpha <- 0 }
          
          if (z < alpha) { # switch made
            recoded0[x,chosen] = new; alleles0[x,chosen] = newallele_length; hidden_crossfamily0[x, chosen] <- crossfamily_new_flag
            mindistance[x,chosenlocus] = newmindistance; alldistance[x,chosenlocus,1:dim(allpossiblerecrud)[1]] = newalldistance
            allrecrf[x,chosenlocus,1:dim(allpossiblerecrud)[1]] = newallrecrf
            recr0[x,chosenlocus] = maxMOI*(chosenlocus-1)+allpossiblerecrud[newclosestrecrud,1]; 
            recrf[x,chosenlocus] = maxMOI*(chosenlocus-1)+allpossiblerecrud[newclosestrecrud, 2]
            recr_repeats0[x,chosenlocus] = sum(recoded0[x,(maxMOI*(chosenlocus-1)+1):(maxMOI*chosenlocus)] == recoded0[x,recr0[x,chosenlocus]])
            recr_repeatsf[x,chosenlocus] = sum(recodedf[x,(maxMOI*(chosenlocus-1)+1):(maxMOI*chosenlocus)] == recodedf[x,recrf[x,chosenlocus]])
          }
        } else {
          old = recoded0[x,chosen]
          n_alleles_locus <- frequencies_RR$n_alleles[chosenlocus]
          if(is.na(n_alleles_locus) || n_alleles_locus == 0) 
          return(list(hidden0=hidden0, hiddenf=hiddenf, recoded0=recoded0, recodedf=recodedf, 
                      alleles0=alleles0, allelesf=allelesf, mindistance=mindistance, alldistance=alldistance, 
                      allrecrf=allrecrf, recr0=recr0, recrf=recrf, recr_repeats0=recr_repeats0, recr_repeatsf=recr_repeatsf,
                      hidden_crossfamily0=hidden_crossfamily0, hidden_crossfamilyf=hidden_crossfamilyf))
          new = sample(1:n_alleles_locus, 1)
          numerator <- frequencies_RR$freq_matrix[chosenlocus, new]
          denominator <- frequencies_RR$freq_matrix[chosenlocus, old]
          
          alpha <- if (!is.na(denominator) && denominator > 0) { numerator / denominator } else { 0 }
          
          if (z < alpha) { # switch made
            recoded0[x,chosen] = new
            sd_val <- frequencies_RR$variability[chosenlocus]
            newallele_length = mean(alleles_definitions_RR[[chosenlocus]][new,]) + stats::rnorm(1,mean=0,sd=sd_val)
            alleles0[x,chosen] = newallele_length        
          }
        }
      } else { # day f hidden allele
        chosen = chosen - nloci*maxMOI
        chosenlocus = ceiling(chosen/maxMOI)
        if (is_locus_comparable[x, chosenlocus]) {
           old = recodedf[x,chosen]; 
           n_alleles_locus <- frequencies_RR$n_alleles[chosenlocus]
           if(is.na(n_alleles_locus) || n_alleles_locus == 0) 
           return(list(hidden0=hidden0, hiddenf=hiddenf, recoded0=recoded0, recodedf=recodedf, 
                      alleles0=alleles0, allelesf=allelesf, mindistance=mindistance, alldistance=alldistance, 
                      allrecrf=allrecrf, recr0=recr0, recrf=recrf, recr_repeats0=recr_repeats0, 
                      recr_repeatsf=recr_repeatsf,hidden_crossfamily0=hidden_crossfamily0, hidden_crossfamilyf=hidden_crossfamilyf))
           new = sample(1:n_alleles_locus, 1); sd_val <- frequencies_RR$variability[chosenlocus]; 
           newallele_length = mean(alleles_definitions_RR[[chosenlocus]][new,]) + stats::rnorm(1,mean=0,sd=sd_val); 
           oldalleles = recodedf[x,intersect(((chosenlocus-1)*maxMOI+1):((chosenlocus)*maxMOI),which(hiddenf[x,] == 0))]; 
           repeatedold = qq; 
           repeatednew = qq
           if (sum(oldalleles == old) >= 1) { repeatedold = 1; }; if (sum(oldalleles == new) >= 1) { repeatednew = 1; }
           crossfamily_new_flag <- 0; penalty_numerator <- 1; 
           penalty_denominator <- qq_crossfamily^hidden_crossfamilyf[x, chosen]
           allpossiblerecrud = expand.grid(1:MOI0[x],1:MOIf[x]); 
           tempalleles = allelesf[x,maxMOI*(chosenlocus-1)+1:maxMOI]; 
           tempalleles[chosen-(chosenlocus-1)*maxMOI] = newallele_length; 
           temprecoded = recodedf[x,maxMOI*(chosenlocus-1)+1:maxMOI]; 
           temprecoded[chosen-(chosenlocus-1)*maxMOI] = new
           newclosestrecrud = which.min(sapply(1:dim(allpossiblerecrud)[1], function (y) abs(tempalleles[allpossiblerecrud[y,2]] - alleles0[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,1]])))
           newmindistance = abs(tempalleles[allpossiblerecrud[newclosestrecrud,2]] - alleles0[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[newclosestrecrud,1]])
           newalldistance = sapply(1:dim(allpossiblerecrud)[1], function (y) abs(tempalleles[allpossiblerecrud[y,2]] - alleles0[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,1]])); 
           newallrecrf = temprecoded[allpossiblerecrud[,2]]
           epsilon <- 1e-9; 
           method <- marker_info$binning_method[chosenlocus]
           calculate_likelihood <- function(distances_vec, recrf_vec, freq_mat, correction_mat, method, threshold) {
               if (method == "microsatellite") { num_probs <- dvect[round(distances_vec) + 1] } else if (method == "cluster") { num_probs <- ifelse(distances_vec <= threshold, 1 - qq_crossfamily, qq_crossfamily) } else { num_probs <- ifelse(distances_vec == 0, 1 - qq_crossfamily, qq_crossfamily) }
               den_probs <- sapply(recrf_vec, function(recr_allele) { if (is.na(recr_allele)) return(NA); 
               if (method == "microsatellite") { 
                pop_probs <- dvect[correction_mat[, recr_allele] + 1] 
                } else if (method == "cluster") { 
                  pop_probs <- ifelse(correction_mat[, recr_allele] <= threshold, 1 - qq_crossfamily, qq_crossfamily) 
                  } else { 
                    pop_probs <- ifelse(correction_mat[, recr_allele] == 0, 1 - qq_crossfamily, qq_crossfamily) }; 
                    sum(freq_mat * pop_probs, na.rm = TRUE) })
               mean(num_probs / (den_probs + epsilon), na.rm = TRUE)
           }
           threshold <- if(method == "cluster") marker_info$cluster_gap_threshold[chosenlocus] else NA
           likelihoodnew <- calculate_likelihood(newalldistance, newallrecrf, frequencies_RR$freq_matrix[chosenlocus, 1:n_alleles_locus], correction_distance_matrix[[chosenlocus]], method, threshold) * repeatednew * penalty_numerator
           likelihoodold <- calculate_likelihood(alldistance[x, chosenlocus, ], allrecrf[x, chosenlocus, ], frequencies_RR$freq_matrix[chosenlocus, 1:n_alleles_locus], correction_distance_matrix[[chosenlocus]], method, threshold) * repeatedold * penalty_denominator
           alpha <- 0
           if (!is.na(likelihoodold) && is.finite(likelihoodold) && likelihoodold > 0) { 
            alpha <- likelihoodnew / likelihoodold 
            } else if (!is.na(likelihoodnew) && likelihoodnew > 0) 
            { alpha <- 2.0 }
           if (is.na(alpha) || !is.finite(alpha)) { alpha <- 0 }
           if (z < alpha) {
            recodedf[x,chosen] = new; allelesf[x,chosen] = newallele_length; hidden_crossfamilyf[x, chosen] <- crossfamily_new_flag
            mindistance[x,chosenlocus] = newmindistance; alldistance[x,chosenlocus,1:dim(allpossiblerecrud)[1]] = newalldistance
            allrecrf[x,chosenlocus,1:dim(allpossiblerecrud)[1]] = newallrecrf
            recr0[x,chosenlocus] = maxMOI*(chosenlocus-1)+allpossiblerecrud[newclosestrecrud,1]; 
            recrf[x,chosenlocus] = maxMOI*(chosenlocus-1)+allpossiblerecrud[newclosestrecrud,2]
             recr_repeats0[x,chosenlocus] = sum(recoded0[x,(maxMOI*(chosenlocus-1)+1):(maxMOI*chosenlocus)] == recoded0[x,recr0[x,chosenlocus]])
            recr_repeatsf[x,chosenlocus] = sum(recodedf[x,(maxMOI*(chosenlocus-1)+1):(maxMOI*chosenlocus)] == recodedf[x,recrf[x,chosenlocus]])
           }
        } else {
          old = recodedf[x,chosen]
          n_alleles_locus <- frequencies_RR$n_alleles[chosenlocus]
          if(is.na(n_alleles_locus) || n_alleles_locus == 0) return(list(hidden0=hidden0, hiddenf=hiddenf, recoded0=recoded0, recodedf=recodedf, alleles0=alleles0, allelesf=allelesf, mindistance=mindistance, alldistance=alldistance, allrecrf=allrecrf, recr0=recr0, recrf=recrf, recr_repeats0=recr_repeats0, recr_repeatsf=recr_repeatsf,hidden_crossfamily0=hidden_crossfamily0, hidden_crossfamilyf=hidden_crossfamilyf))
          new = sample(1:n_alleles_locus, 1)
          numerator <- frequencies_RR$freq_matrix[chosenlocus, new]
          denominator <- frequencies_RR$freq_matrix[chosenlocus, old]
          
          alpha <- if (!is.na(denominator) && denominator > 0) { numerator / denominator } else { 0 }
          
          if (z < alpha) { # switch made
            recodedf[x,chosen] = new
            sd_val <- frequencies_RR$variability[chosenlocus]
            newallele_length = mean(alleles_definitions_RR[[chosenlocus]][new,]) + stats::rnorm(1,mean=0,sd=sd_val)
            allelesf[x,chosen] = newallele_length
          }
        }
      }
    }
  }
  return(list(
    hidden0=hidden0, hiddenf=hiddenf, recoded0=recoded0, recodedf=recodedf, 
    alleles0=alleles0, allelesf=allelesf, mindistance=mindistance, alldistance=alldistance, 
    allrecrf=allrecrf, recr0=recr0, recrf=recrf, 
    recr_repeats0=recr_repeats0, recr_repeatsf=recr_repeatsf,
    hidden_crossfamily0=hidden_crossfamily0, hidden_crossfamilyf=hidden_crossfamilyf
  ))
}