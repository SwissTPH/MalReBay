#' Update a Single Hidden Allele for Amplicon Data via Metropolis-Hastings
#'
#' @description
#' A core MCMC sampler step that updates the imputed identity of a single, randomly
#' chosen hidden allele for one patient. This function is specifically for
#' amplicon sequencing (ampseq) data and uses a Metropolis-Hastings algorithm
#' to sample from the complex posterior distribution of the hidden allele's state.
#'
#' @details
#' This function executes one update for a single patient (`x`) within a Gibbs
#' sampling iteration. Its workflow is as follows:
#' \enumerate{
#'   \item **Selection:** Randomly selects one "hidden" (unobserved) allele
#'     slot from either the Day 0 or Day of Failure sample of the specified patient.
#'   \item **Proposal:** Proposes a new allele identity (a new character sequence)
#'     for this slot by sampling from the population-level allele frequency distribution.
#'   \item **Likelihood Calculation:** It calculates the log-likelihood of the
#'     patient's entire genetic profile under both the current state and the
#'     proposed new state. The likelihood logic depends on the patient's current
#'     `classification`:
#'     \itemize{
#'       \item **If Recrudescence:** The likelihood is based on whether alleles
#'         at Day of Failure are present at Day 0, modulated by the genotyping
#'         error/mutation probability (`qq`).
#'       \item **If Reinfection:** The likelihood is the product of the
#'         population frequencies of the observed Day of Failure alleles.
#'     }
#'   \item **Decision & Update:** The Metropolis-Hastings acceptance ratio is
#'     computed using the change in log-likelihood and the proposal probabilities
#'     (the frequencies of the old vs. new alleles). The proposed change is accepted
#'     stochastically based on this ratio. If accepted, the relevant `recoded0`
#'     or `recodedf` matrix is updated with the new allele.
#' }
#'
#' @param x The integer index of the patient to be updated.
#' @param hidden0,hiddenf The current state matrices of boolean flags indicating
#'   which allele slots are imputed (hidden).
#' @param recoded0,recodedf The current state matrices of character allele IDs for
#'   all samples (both observed and imputed).
#' @param classification A vector of the current classification state for all patients.
#' @param qq The current estimate of the genotyping error/mutation probability.
#' @param frequencies_RR The list object containing population genetic statistics,
#'   including allele codes (`allele_codes`) and their frequencies (`freq_matrix`).
#' @param nloci,maxMOI Core dimensions of the analysis (number of loci and max MOI).
#'
#' @return A list containing the updated `recoded0` and `recodedf` state matrices.
#'   If the proposed switch was accepted for patient `x`, the returned matrices
#'   will reflect this update; otherwise, they are returned unchanged.
#'
#' @keywords internal
#' @noRd
#'

switch_hidden_ampseq <- function(x, hidden0, hiddenf, recoded0, recodedf,
                                 classification, qq, q_loss, frequencies_RR, nloci, maxMOI) { 
  
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
  
  # new_id <- sample(possible_alleles, 1)
  
  allele_freqs <- frequencies_RR$freq_matrix[chosenlocus, 1:length(possible_alleles)]
  if(any(is.na(allele_freqs)) || sum(allele_freqs) == 0) {
    new_id <- sample(possible_alleles, 1)
  } else {
    new_id <- sample(possible_alleles, 1, prob = allele_freqs)
  }
  
  # log-likelihood for a given state
  calculate_log_lik <- function(patient_recoded0, patient_recodedf) {
    log_lik_total <- 0
    for (locus_idx in 1:nloci) {
      d0_alleles <- unique(patient_recoded0[(maxMOI*(locus_idx-1)+1):(maxMOI*locus_idx)])
      d0_alleles <- d0_alleles[!is.na(d0_alleles)]
      
      df_alleles <- unique(patient_recodedf[(maxMOI*(locus_idx-1)+1):(maxMOI*locus_idx)])
      df_alleles <- df_alleles[!is.na(df_alleles)]
      
      if (length(d0_alleles) == 0 || length(df_alleles) == 0) next
      
      if (classification[x] == 1) { # RECRUDESCENCE
        log_probs_per_allele <- sapply(df_alleles, function(allele) {
          log(ifelse(allele %in% d0_alleles, 1 - qq, qq) + 1e-10)
        })
        log_lik_total <- log_lik_total + sum(log_probs_per_allele)
        
        lost_alleles <- d0_alleles[!d0_alleles %in% df_alleles]
        if (length(lost_alleles) > 0) {
          log_lik_total <- log_lik_total + (length(lost_alleles) * log(q_loss + 1e-10))
        }  
      } else { # REINFECTION
        freqs <- frequencies_RR$freq_matrix[locus_idx, ]
        allele_codes <- frequencies_RR$allele_codes[[locus_idx]]
        matched_indices <- match(df_alleles, allele_codes)
        
        if (any(is.na(matched_indices))) {
          log_lik_total <- log_lik_total - 1e6 
        } else {
          log_lik_total <- log_lik_total + sum(log(freqs[matched_indices] + 1e-10))
        }
      }
    }
    return(log_lik_total)
  }
  
  log_lik_old <- calculate_log_lik(recoded0[x,], recodedf[x,])
  recoded0_new <- recoded0[x,]
  recodedf_new <- recodedf[x,]
  if (is_day0_allele) {
    recoded0_new[chosen] <- new_id
  } else {
    recodedf_new[chosen - ncol(hidden0)] <- new_id
  }
  
  log_lik_new <- calculate_log_lik(recoded0_new, recodedf_new)
  old_freq <- frequencies_RR$freq_matrix[chosenlocus, which(possible_alleles == old_id)]
  if(length(old_freq)==0) old_freq <- 0
  new_freq <- frequencies_RR$freq_matrix[chosenlocus, which(possible_alleles == new_id)]
  if(length(new_freq)==0) new_freq <- 0
  log_acceptance_ratio <- log_lik_new - log_lik_old
  
  alpha <- min(1, exp(log_acceptance_ratio))
  
  if (z < alpha) {
    if (is_day0_allele) {
      recoded0[x, chosen] <- new_id
    } else {
      recodedf[x, chosen - ncol(hidden0)] <- new_id
    }
  }
  
  return(list(recoded0=recoded0, recodedf=recodedf))
}