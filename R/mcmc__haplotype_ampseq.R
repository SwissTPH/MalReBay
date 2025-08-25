#' MCMC Chain for Amplicon Sequencing Data
#'
#' @description
#' Implements the full Gibbs sampling algorithm for a single MCMC chain to
#' classify infections as recrudescence or reinfection, specifically tailored for
#' amplicon sequencing (ampseq) data.
#'
#' @details
#' This function serves as the engine for the Bayesian classification model.
#' It initializes the model's state and iteratively updates each component in a
#' Gibbs sampling framework. The key steps within each iteration are:
#' \enumerate{
#'   \item **Data Preparation:** Extracts and recodes observed alleles for Day 0
#'     and Day of Failure samples. The Multiplicity of Infection (MOI) for each
#'     sample is determined.
#'   \item **Initial Imputation:** Guesses the identity of unobserved ("hidden")
#'     alleles for infections where the MOI is greater than the number of observed
#'     alleles.
#'   \item **Classification Update:** Updates the classification (recrudescence vs.
#'     reinfection) for each patient by calculating a likelihood ratio based on
#'     the current state of observed and imputed alleles.
#'   \item **Hidden Allele Update:** Re-samples the identity of the hidden alleles
#'     based on the current classification and allele frequencies.
#'   \item **Parameter Update:** Updates model hyperparameters, including the
#'     error rate (`qq`) and population-level allele frequencies, using their
#'     respective posterior distributions.
#' }
#' The function records the state of key parameters after the burn-in period at
#' specified intervals. It is designed to be run in parallel with other identical,
#' independent chains.
#'
#' @param chain_id An integer identifying the chain, mainly for logging purposes.
#' @param nruns The total number of MCMC iterations to perform.
#' @param burnin The number of initial iterations to discard as burn-in.
#' @param record_interval The interval at which to record the chain's state.
#' @param nids,ids,nloci,maxMOI,locinames Core dimensions and metadata for the
#'   analysis: number of patient IDs, the IDs themselves, number of loci, max
#'   MOI, and locus names.
#' @param genotypedata_RR,additional_neutral Data frames containing the primary
#'   paired samples and additional baseline samples, respectively.
#' @param marker_info A data frame containing metadata for each genetic marker.
#' @param is_locus_comparable A logical matrix indicating for each patient and
#'   locus whether a comparison is possible (i.e., data exists for both Day 0
#'   and Day of Failure).
#' @param record_hidden_alleles A logical flag. If `TRUE`, the full state of
#'   imputed hidden alleles is saved at each recorded step.
#'
#' @return A list containing the full history of the MCMC chain after burn-in.
#'   The list includes the following components:
#'   \item{state_classification}{Matrix of classification states (1 for
#'     recrudescence, 0 for reinfection).}
#'   \item{state_parameters}{Matrix of key model hyperparameters over time.}
#'   \item{state_loglikelihood}{Vector of the overall log-likelihood at each
#'     recorded step.}
#'   \item{state_recoded0, state_recodedf}{Optional arrays of the imputed alleles
#'     for Day 0 and Day of Failure samples if `record_hidden_alleles` is `TRUE`.}
#'   \item{locus_lrs, locus_dists}{Arrays containing the per-locus likelihood
#'     ratios and allele distances for each patient at each recorded step.}
#'
#' @keywords internal
#' @noRd
#'

run_one_chain_ampseq <- function(chain_id,
                                 nruns, burnin, record_interval,
                                 nids, ids, nloci, maxMOI, locinames,
                                 genotypedata_RR, additional_neutral,
                                 marker_info,
                                 is_locus_comparable,
                                 record_hidden_alleles = FALSE)
{
 
  MOI0 <- rep(0, nids)
  MOIf <- rep(0, nids)
  recoded0 <- matrix(NA_character_, nids, maxMOI * nloci)
  recodedf <- matrix(NA_character_, nids, maxMOI * nloci)
  hidden0 <- matrix(NA_integer_, nids, maxMOI * nloci)
  hiddenf <- matrix(NA_integer_, nids, maxMOI * nloci)
  mindistance <- matrix(NA_integer_, nids, nloci)

  for (i in 1:nids) {
    patient_id <- ids[i]
    row_idx_d0 <- which(genotypedata_RR$Sample.ID == paste(patient_id, "Day 0"))
    row_idx_df <- which(genotypedata_RR$Sample.ID == paste(patient_id, "Day Failure"))
    
    if (length(row_idx_d0) == 0 || length(row_idx_df) == 0) { next }
    patient_max_moi0 <- 0
    patient_max_moif <- 0
  
    for (j in 1:nloci) {
        locus <- locinames[j]
        locicolumns <- grepl(paste0("^", locus, "_"), colnames(genotypedata_RR))

        d0_alleles <- as.character(genotypedata_RR[row_idx_d0[1], locicolumns])
        d0_alleles <- d0_alleles[!is.na(d0_alleles) & d0_alleles != "NA"]
        patient_max_moi0 <- max(patient_max_moi0, length(d0_alleles))
        if (length(d0_alleles) > 0) {
          recoded0[i, (maxMOI * (j - 1) + 1):(maxMOI * (j - 1) + length(d0_alleles))] <- d0_alleles
        }

        df_alleles <- as.character(genotypedata_RR[row_idx_df[1], locicolumns])
        df_alleles <- df_alleles[!is.na(df_alleles) & df_alleles != "NA"]
        patient_max_moif <- max(patient_max_moif, length(df_alleles))
        if (length(df_alleles) > 0) {
          recodedf[i, (maxMOI * (j - 1) + 1):(maxMOI * (j - 1) + length(df_alleles))] <- df_alleles
        }
    }
    MOI0[i] <- patient_max_moi0
    MOIf[i] <- patient_max_moif
  }

  is_observed0 <- !is.na(recoded0)
  is_observedf <- !is.na(recodedf)
  original_recoded0 <- recoded0
  original_recodedf <- recodedf


  recoded_additional_neutral <- matrix(NA_character_, 0, maxMOI * nloci)
  if (!is.null(additional_neutral) && nrow(additional_neutral) > 0) {
  }

  frequencies_RR <- calculate_frequencies3(
    genotypedata = rbind(genotypedata_RR, additional_neutral),
    alleles_definitions = NULL,
    marker_info = marker_info
  )

  hidden0 <- matrix(NA_integer_, nids, maxMOI * nloci)
  hiddenf <- matrix(NA_integer_, nids, maxMOI * nloci)

  # Initial imputation of hidden alleles
  for (i in 1:nids) {
    for (j in 1:nloci) {
        d0_cols <- (maxMOI*(j-1)+1):(maxMOI*j)
        n_observed_d0 <- sum(!is.na(recoded0[i, d0_cols]))
        n_to_impute_d0 <- MOI0[i] - n_observed_d0
        if (n_to_impute_d0 > 0) {
            indices_to_fill_d0 <- d0_cols[is.na(recoded0[i, d0_cols])][1:n_to_impute_d0]
            hidden0[i, indices_to_fill_d0] <- 1
            possible_alleles <- frequencies_RR$allele_codes[[j]]
            if (length(possible_alleles) > 0) {
                allele_freqs <- frequencies_RR$freq_matrix[j, 1:frequencies_RR$n_alleles[j]]
                new_alleles <- sample(possible_alleles, n_to_impute_d0, replace=TRUE, prob=allele_freqs)
                recoded0[i, indices_to_fill_d0] <- new_alleles
            }
        }
        hidden0[i, d0_cols[!is.na(recoded0[i, d0_cols])]] <- 0
        
        df_cols <- (maxMOI*(j-1)+1):(maxMOI*j)
        n_observed_df <- sum(!is.na(recodedf[i, df_cols]))
        n_to_impute_df <- MOIf[i] - n_observed_df
        if(n_to_impute_df > 0) {
            indices_to_fill_df <- df_cols[is.na(recodedf[i, df_cols])][1:n_to_impute_df]
            hiddenf[i, indices_to_fill_df] <- 1
            possible_alleles <- frequencies_RR$allele_codes[[j]]
            if (length(possible_alleles) > 0) {
                allele_freqs <- frequencies_RR$freq_matrix[j, 1:frequencies_RR$n_alleles[j]]
                new_alleles <- sample(possible_alleles, n_to_impute_df, replace=TRUE, prob=allele_freqs)
                recodedf[i, indices_to_fill_df] <- new_alleles
            }
        }
        hiddenf[i, df_cols[!is.na(recodedf[i, df_cols])]] <- 0
    }
  }

  # Initialize parameters and classification
  qq <- mean(c(hidden0, hiddenf), na.rm = TRUE); if (is.na(qq)) qq <- 0.1
  q_loss <- 0.1
  prob_recrud <- 0.5
  classification <- ifelse(stats::runif(nids) < prob_recrud, 1, 0)
  num_records <- floor((nruns - burnin) / record_interval)
  
  classification_history <- matrix(NA, nids, num_records)
  parameters_history <- matrix(NA, 4 + (2 * nloci), num_records)
  loglikelihood_history <- rep(NA_real_, num_records)
  locus_lrs_history <- array(NA, c(nids, nloci, num_records))
  locus_dists_history <- array(NA, c(nids, nloci, num_records))
  
  alleles0_history <- if (record_hidden_alleles) array(NA, c(nids, maxMOI * nloci, num_records)) else NULL
  allelesf_history <- if (record_hidden_alleles) array(NA, c(nids, maxMOI * nloci, num_records)) else NULL

  # Main MCMC loop
  for (iter in 1:nruns) {
    mindistance <- matrix(NA_integer_, nids, nloci)
    locus_lrs_this_step <- matrix(1.0, nrow = nids, ncol = nloci)

    for (i in 1:nids) {
      for (j in 1:nloci) {
        d0_alleles <- unique(recoded0[i, (maxMOI*(j-1)+1):(maxMOI*j)])
        d0_alleles <- d0_alleles[!is.na(d0_alleles)]
        df_alleles <- unique(recodedf[i, (maxMOI*(j-1)+1):(maxMOI*j)])
        df_alleles <- df_alleles[!is.na(df_alleles)]
        if (length(d0_alleles) == 0 || length(df_alleles) == 0 || !is_locus_comparable[i, j]) {
          next
        }
        probs_per_failure_allele <- sapply(df_alleles, function(allele) {
          ifelse(allele %in% d0_alleles, 1 - qq, qq)
        })
        prob_part1 <- prod(probs_per_failure_allele)
        lost_alleles <- d0_alleles[!d0_alleles %in% df_alleles]
        prob_part2 <- q_loss ^ length(lost_alleles)
        prob_recrud_locus <- prob_part1 * prob_part2

        matched_indices <- match(df_alleles, frequencies_RR$allele_codes[[j]])
        valid_indices <- !is.na(matched_indices)
        
        if (!any(valid_indices)) {
          prob_reinfect_locus <- 1e-10
        } else {
          freqs <- frequencies_RR$freq_matrix[j, matched_indices[valid_indices]]
          prob_reinfect_locus <- prod(freqs, na.rm = TRUE)
        }
        locus_lrs_this_step[i, j] <- prob_recrud_locus / (prob_reinfect_locus + 1e-10)
        mindistance[i, j] <- ifelse(any(df_alleles %in% d0_alleles), 0, 1)
      }
    }

    likelihoodratio <- apply(locus_lrs_this_step, 1, function(r) prod(pmax(r, 1e-10)))
    prior_odds <- prob_recrud / (1 - prob_recrud + 1e-10)
    posterior_odds <- likelihoodratio * prior_odds    
    
    prob_is_recrudescence <- posterior_odds / (1 + posterior_odds)
    classification <- ifelse(stats::runif(nids) < prob_is_recrudescence, 1, 0)
    
    for (i in 1:nids) {
        updated_states <- switch_hidden_ampseq(x = i, hidden0 = hidden0, hiddenf = hiddenf, 
                                               recoded0 = recoded0, recodedf = recodedf,
                                               classification = classification, qq = qq, q_loss = q_loss,
                                               frequencies_RR = frequencies_RR, nloci = nloci, maxMOI = maxMOI)
        temp_recoded0 <- updated_states$recoded0
        temp_recodedf <- updated_states$recodedf
        temp_recoded0[is_observed0] <- original_recoded0[is_observed0]
        temp_recodedf[is_observedf] <- original_recodedf[is_observedf]
        recoded0 <- temp_recoded0
        recodedf <- temp_recodedf
    }

    q_posterior_alpha <- 1 + sum(c(hidden0, hiddenf) == 1, na.rm = TRUE)
    q_posterior_beta  <- 1 + sum(c(hidden0, hiddenf) == 0, na.rm = TRUE)
    qq <- stats::rbeta(1, q_posterior_alpha, q_posterior_beta)
    
    n_lost <- 0
    n_retained <- 0
    recrud_indices <- which(classification == 1)
    if(length(recrud_indices) > 0){
        for (i in recrud_indices) {
            for (j in 1:nloci) {
                d0_alleles <- unique(recoded0[i, (maxMOI*(j-1)+1):(maxMOI*j)])
                d0_alleles <- d0_alleles[!is.na(d0_alleles)]
                df_alleles <- unique(recodedf[i, (maxMOI*(j-1)+1):(maxMOI*j)])
                df_alleles <- df_alleles[!is.na(df_alleles)]
                
                if (length(d0_alleles) == 0) next
                
                n_lost <- n_lost + sum(!d0_alleles %in% df_alleles)
                n_retained <- n_retained + sum(d0_alleles %in% df_alleles)
            }
        }
    }
    q_loss_posterior_alpha <- 1 + n_lost
    q_loss_posterior_beta <- 1 + n_retained
    q_loss <- stats::rbeta(1, q_loss_posterior_alpha, q_loss_posterior_beta)

    tempdata <- rbind(recoded0, recodedf)
    new_freq_matrix <- frequencies_RR$freq_matrix
    for (locus_idx in 1:nloci) {
      new_row <- findposteriorfrequencies(locus_index = locus_idx, tempdata = tempdata, 
                                          maxMOI = maxMOI, frequencies_RR = frequencies_RR)
      n_alleles_locus <- frequencies_RR$n_alleles[locus_idx]
      if (length(new_row) > 0 && !any(is.na(new_row)) && n_alleles_locus > 0) {
        new_freq_matrix[locus_idx, 1:n_alleles_locus] <- new_row
      }
    }
    frequencies_RR$freq_matrix <- new_freq_matrix
    if (iter > burnin && (iter - burnin) %% record_interval == 0) {
      record_idx <- (iter - burnin) / record_interval
      loglik_val <- sum(log(pmax(likelihoodratio, 1e-10)))
      loglikelihood_history[record_idx] <- ifelse(is.finite(loglik_val), loglik_val, NA)
      classification_history[, record_idx] <- classification
      parameters_history[1, record_idx] <- qq
      parameters_history[2, record_idx] <- prob_recrud
      parameters_history[3, record_idx] <- q_loss 
      he_per_locus <- 1 - rowSums(frequencies_RR$freq_matrix^2, na.rm = TRUE)
      parameters_history[(5 + nloci):(5 + 2 * nloci - 1), record_idx] <- he_per_locus
      
      locus_lrs_history[,, record_idx] <- locus_lrs_this_step
      locus_dists_history[,, record_idx] <- mindistance
      
      if (record_hidden_alleles) {
        alleles0_history[,, record_idx] <- recoded0
        allelesf_history[,, record_idx] <- recodedf
      }
    }
  }
  
  return(list(
    state_classification = classification_history,
    state_parameters     = parameters_history,
    state_loglikelihood  = loglikelihood_history,
    state_recoded0       = alleles0_history,
    state_recodedf       = allelesf_history,
    locus_lrs            = locus_lrs_history,
    locus_dists          = locus_dists_history
  ))
}