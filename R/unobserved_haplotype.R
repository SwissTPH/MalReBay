#' Update a Single Hidden Allele for Amplicon Data via Metropolis-Hastings
#'
#' @description
#' Proposes and accepts/rejects a new identity for one randomly chosen hidden
#' (imputed) allele for a single patient, using a Metropolis-Hastings step.
#'
#' @details
#' Executes one hidden allele update for patient \code{x} within a Gibbs
#' sampling iteration:
#' \enumerate{
#'   \item \strong{Selection:} Randomly picks one hidden allele slot from either
#'     the Day 0 or recurrence sample of patient \code{x}.
#'   \item \strong{Proposal:} Samples a new allele identity from the
#'     population-level allele frequency distribution.
#'   \item \strong{Likelihood Calculation:} Computes the log-likelihood of the
#'     patient's full genetic profile under the current and proposed states:
#'     \itemize{
#'       \item \strong{Recrudescence:} Likelihood based on allele sharing between
#'         Day 0 and recurrence, modulated by \code{q_mismatch} and \code{q_loss}.
#'       \item \strong{Reinfection:} Likelihood based on the product of
#'         population frequencies of recurrence alleles.
#'     }
#'   \item \strong{Accept/Reject:} Computes the Metropolis-Hastings acceptance
#'     ratio (including a proposal correction for the frequency-proportional
#'     proposal) and stochastically accepts or rejects the proposed change.
#' }
#'
#' @param x Integer index of the patient to update.
#' @param hidden0,hiddenf Integer matrices flagging which allele slots are
#'   imputed (1) or observed (0) for Day 0 and recurrence samples.
#' @param recoded0,recodedf Character matrices of allele IDs (observed and
#'   imputed) for Day 0 and recurrence samples.
#' @param classification Integer vector of current classifications for all
#'   patients (1 = recrudescence, 0 = reinfection).
#' @param q_mismatch Current estimate of the sequencing error/mismatch probability
#'   (\code{q_mismatch} in the calling chain).
#' @param q_loss Current estimate of the allele loss probability between
#'   Day 0 and recurrence.
#' @param allele_frequencies List containing allele codes (\code{allele_codes})
#'   and their population frequencies (\code{freq_matrix}).
#' @param nloci,maxMOI Number of loci and maximum MOI.
#'
#' @return A list with updated \code{recoded0} and \code{recodedf} matrices.
#'   Only the allele slot of patient \code{x} may differ from the input if the
#'   proposed change was accepted; otherwise returned unchanged.
#'
#' @noRd
switch_hidden_ampseq <- function(x, hidden0, hiddenf, recoded0, recodedf,
                                 classification, q_mismatch, q_loss, allele_frequencies,
                                 nloci, maxMOI) {
  
  # Identify hidden (imputed) allele positions for this patient
  hidden_positions <- which(c(hidden0[x, ], hiddenf[x, ]) == 1)
  if (length(hidden_positions) == 0)
    return(list(recoded0 = recoded0, recodedf = recodedf))
  
  # Randomly select one hidden position to propose a new value for
  chosen         <- sample(hidden_positions, 1)
  is_day0_allele <- chosen <= ncol(hidden0)
  
  if (is_day0_allele) {
    chosenlocus <- floor((chosen - 1) / maxMOI) + 1
    old_id      <- recoded0[x, chosen]
  } else {
    chosenlocus <- floor((chosen - 1 - ncol(hidden0)) / maxMOI) + 1
    old_id      <- recodedf[x, chosen - ncol(hidden0)]
  }
  
  possible_alleles <- allele_frequencies$allele_codes[[chosenlocus]]
  if (length(possible_alleles) == 0)
    return(list(recoded0 = recoded0, recodedf = recodedf))
  
  # Proposal distribution: population frequencies, normalised
  allele_freqs <- allele_frequencies$freq_matrix[chosenlocus, 1:length(possible_alleles)]
  if (any(is.na(allele_freqs)) || sum(allele_freqs) == 0)
    allele_freqs <- rep(1 / length(possible_alleles), length(possible_alleles))
  
  new_id <- sample(possible_alleles, 1, prob = allele_freqs)
  
  calculate_log_lik <- function(patient_recoded0, patient_recodedf) {
    log_lik_total <- 0
    
    for (locus_idx in 1:nloci) {
      d0_alleles <- patient_recoded0[(maxMOI*(locus_idx-1)+1):(maxMOI*locus_idx)]
      d0_alleles <- d0_alleles[!is.na(d0_alleles)]
      df_alleles <- patient_recodedf[(maxMOI*(locus_idx-1)+1):(maxMOI*locus_idx)]
      df_alleles <- df_alleles[!is.na(df_alleles)]
      
      if (length(d0_alleles) == 0 || length(df_alleles) == 0) next
      
      allele_codes <- allele_frequencies$allele_codes[[locus_idx]]
      freq_row     <- allele_frequencies$freq_matrix[locus_idx, ]
      
      if (classification[x] == 1) {
        # RECRUDESCENCE: pair-averaged likelihood accounting for mismatch (q_mismatch) and allele loss (q_loss)
        n_pairs  <- length(d0_alleles) * length(df_alleles)
        n_lost   <- sum(!d0_alleles %in% df_alleles)
        
        pair_sum <- 0
        for (r0 in seq_along(d0_alleles)) {
          for (rf in seq_along(df_alleles)) {
            w   <- ifelse(d0_alleles[r0] == df_alleles[rf], 1 - q_mismatch, q_mismatch)
            idx <- match(df_alleles[rf], allele_codes)
            # Protect against zero frequencies
            lam <- if (!is.na(idx)) max(freq_row[idx], 1e-6) else 1e-6
            pair_sum <- pair_sum + w / lam
          }
        }
        
        log_lik_total <- log_lik_total +
          n_lost * log(q_loss + 1e-10) +
          log(pair_sum / n_pairs + 1e-10)
        
      } else {
        # REINFECTION: log product of recurrence allele frequencies
        df_matched <- match(df_alleles, allele_codes)
        df_valid   <- !is.na(df_matched)
        if (!any(df_valid)) {
          log_lik_total <- log_lik_total - 1e6
        } else {
          log_lik_total <- log_lik_total +
            sum(log(pmax(freq_row[df_matched[df_valid]], 1e-6)))
        }
      }
    }
    return(log_lik_total)
  }
  
  # Evaluate likelihood of current and proposed states
  log_lik_old  <- calculate_log_lik(recoded0[x, ], recodedf[x, ])
  
  recoded0_new <- recoded0[x, ]
  recodedf_new <- recodedf[x, ]
  if (is_day0_allele) {
    recoded0_new[chosen] <- new_id
  } else {
    recodedf_new[chosen - ncol(hidden0)] <- new_id
  }
  
  log_lik_new <- calculate_log_lik(recoded0_new, recodedf_new)
  
  # MH acceptance ratio with proposal correction.
  # Since proposals are drawn proportional to allele frequency (independent of current state):
  #   log q(old|new) - log q(new|old) = log(freq_old) - log(freq_new)
  # Note: proposing the same allele as current (new_id == old_id) is possible,
  # in which case the acceptance ratio is 1 and the step is a no-op.
  old_idx  <- match(old_id, possible_alleles)
  new_idx  <- match(new_id, possible_alleles)
  old_freq <- if (!is.na(old_idx)) allele_freqs[old_idx] else 1e-10
  new_freq <- if (!is.na(new_idx)) allele_freqs[new_idx] else 1e-10
  
  log_proposal_ratio   <- log(old_freq + 1e-10) - log(new_freq + 1e-10)
  log_acceptance_ratio <- log_lik_new - log_lik_old + log_proposal_ratio
  
  if (stats::runif(1) < min(1, exp(log_acceptance_ratio))) {
    if (is_day0_allele) {
      recoded0[x, chosen] <- new_id
    } else {
      recodedf[x, chosen - ncol(hidden0)] <- new_id
    }
  }
  
  return(list(recoded0 = recoded0, recodedf = recodedf))
}