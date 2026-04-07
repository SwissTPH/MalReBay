#' MCMC Chain for Amplicon Sequencing Data
#'
#' @description
#' Runs a single Gibbs sampling MCMC chain to classify paired amplicon
#' sequencing samples as recrudescence or reinfection.
#'
#' @details
#' Each iteration performs the following steps:
#' \enumerate{
#'   \item \strong{Data Preparation:} Extracts observed alleles and determines
#'     MOI per locus for Day 0 and recurrence samples.
#'   \item \strong{Initial Imputation:} Samples hidden alleles where MOI exceeds
#'     the number of observed alleles.
#'   \item \strong{Classification Update:} Updates recrudescence/reinfection
#'     classification via likelihood ratios over observed and imputed alleles.
#'   \item \strong{Hidden Allele Update:} Re-samples hidden alleles conditional
#'     on current classification and allele frequencies.
#'   \item \strong{Parameter Update:} Updates \code{q_mismatch} (sequencing
#'     error rate), \code{q_dropout} (allelic dropout rate), \code{q_loss}
#'     (allele loss rate), and population-level allele frequencies from their
#'     posterior distributions.
#' }
#' The prior probability of recrudescence is fixed at 0.5. States are recorded
#' after burn-in at every \code{record_interval} iterations.
#'
#' @param chain_id Integer identifier for the chain, used for logging.
#' @param nruns Total number of MCMC iterations.
#' @param burnin Number of initial iterations to discard.
#' @param record_interval Number of iterations between recorded states.
#' @param nids,ids,nloci,maxMOI,locinames Number of patients, patient IDs,
#'   number of loci, maximum MOI, and locus names.
#' @param late_failures_site Data frame of paired Day 0 and recurrence samples.
#' @param additional_site Data frame of additional samples used to inform
#'   allele frequency priors.
#' @param marker_info Data frame of metadata for each genetic marker.
#' @param is_locus_comparable Logical matrix (\code{nids x nloci}) indicating
#'   whether data exists for both timepoints at each patient-locus combination.
#' @param record_hidden_alleles If \code{TRUE}, saves the full allele state
#'   arrays at each recorded step.
#'
#' @return A list with the following components:
#'   \item{state_classification}{Matrix of classifications (1 = recrudescence,
#'     0 = reinfection) across recorded steps.}
#'   \item{state_parameters}{3 x \code{num_records} matrix of sampled parameter
#'     values: row 1 = \code{q_mismatch}, row 2 = \code{q_dropout},
#'     row 3 = \code{q_loss}.}
#'   \item{state_loglikelihood}{Vector of log-likelihoods at each recorded step.}
#'   \item{state_recoded0, state_recodedf}{Full allele state arrays for Day 0
#'     and recurrence if \code{record_hidden_alleles = TRUE}, else \code{NULL}.}
#'   \item{locus_lrs}{Array of per-locus likelihood ratios per patient per
#'     recorded step.}
#'   \item{locus_dists}{Array of binary match indicators (0 = at least one
#'     allele matches, 1 = no match) per patient, locus, and recorded step.}
#'
#' @noRd
run_one_chain_ampseq <- function(chain_id,
                                 nruns, burnin, record_interval,
                                 nids, ids, nloci, maxMOI, locinames,
                                 late_failures_site, additional_site,
                                 marker_info,
                                 is_locus_comparable,
                                 record_hidden_alleles = FALSE)
{
  # Data extraction per locus
  MOI0_locus  <- matrix(0, nids, nloci)
  MOIf_locus  <- matrix(0, nids, nloci)
  recoded0    <- matrix(NA_character_, nids, maxMOI * nloci)
  recodedf    <- matrix(NA_character_, nids, maxMOI * nloci)
  
  for (i in 1:nids) {
    patient_id <- ids[i]
    row_idx_d0 <- which(late_failures_site$Sample.ID == paste(patient_id, "Day 0"))
    row_idx_df <- which(late_failures_site$Sample.ID == paste(patient_id, "recurrence"))
    if (length(row_idx_d0) == 0 || length(row_idx_df) == 0) next
    
    for (j in 1:nloci) {
      locus       <- locinames[j]
      locicolumns <- grepl(paste0("^", locus, "_"), colnames(late_failures_site))
      
      d0_alleles <- as.character(late_failures_site[row_idx_d0[1], locicolumns])
      d0_alleles <- d0_alleles[!is.na(d0_alleles) & d0_alleles != "NA"]
      MOI0_locus[i, j] <- length(d0_alleles)
      if (length(d0_alleles) > 0)
        recoded0[i, (maxMOI*(j-1)+1):(maxMOI*(j-1)+length(d0_alleles))] <- d0_alleles
      
      df_alleles <- as.character(late_failures_site[row_idx_df[1], locicolumns])
      df_alleles <- df_alleles[!is.na(df_alleles) & df_alleles != "NA"]
      MOIf_locus[i, j] <- length(df_alleles)
      if (length(df_alleles) > 0)
        recodedf[i, (maxMOI*(j-1)+1):(maxMOI*(j-1)+length(df_alleles))] <- df_alleles
    }
  }
  
  is_observed0      <- !is.na(recoded0)
  is_observedf      <- !is.na(recodedf)
  original_recoded0 <- recoded0
  original_recodedf <- recodedf
  
  # processing the additional data
  recoded_additional_site <- matrix(NA_character_, 0, maxMOI * nloci)
  if (!is.null(additional_site) && nrow(additional_site) > 0) {
    recoded_additional_site <- matrix(NA_character_, nrow(additional_site), maxMOI * nloci)
    for (i in 1:nrow(additional_site)) {
      for (j in 1:nloci) {
        locus         <- locinames[j]
        locicolumns   <- grepl(paste0("^", locus, "_"), colnames(additional_site))
        extra_alleles <- as.character(additional_site[i, locicolumns])
        extra_alleles <- extra_alleles[!is.na(extra_alleles) & extra_alleles != "NA"]
        if (length(extra_alleles) > 0) {
          n_to_fill    <- min(length(extra_alleles), maxMOI)
          cols_to_fill <- (maxMOI*(j-1)+1):(maxMOI*(j-1)+n_to_fill)
          recoded_additional_site[i, cols_to_fill] <- extra_alleles[1:n_to_fill]
        }
      }
    }
  }
  
  allele_frequencies <- calculate_frequencies(
    genotypedata        = rbind(late_failures_site, additional_site),
    alleles_definitions = NULL,
    marker_info         = marker_info
  )
  
  # Initializing hidden alleles
  hidden0 <- matrix(NA_integer_, nids, maxMOI * nloci)
  hiddenf <- matrix(NA_integer_, nids, maxMOI * nloci)
  
  for (i in 1:nids) {
    for (j in 1:nloci) {
      if (!is_locus_comparable[i, j]) next
      
      # Day 0
      d0_cols <- (maxMOI*(j-1)+1):(maxMOI*j)
      n_obs0  <- sum(!is.na(recoded0[i, d0_cols]))
      n_imp0  <- MOI0_locus[i, j] - n_obs0
      if (n_imp0 > 0) {
        fill0           <- d0_cols[is.na(recoded0[i, d0_cols])][1:n_imp0]
        hidden0[i, fill0] <- 1L
        possible        <- allele_frequencies$allele_codes[[j]]
        freqs           <- allele_frequencies$freq_matrix[j, 1:allele_frequencies$n_alleles[j]]
        recoded0[i, fill0] <- sample(possible, n_imp0, replace = TRUE, prob = freqs)
      }
      hidden0[i, d0_cols[!is.na(recoded0[i, d0_cols])]] <- 0L
      
      # Recurrence
      df_cols <- (maxMOI*(j-1)+1):(maxMOI*j)
      n_obsf  <- sum(!is.na(recodedf[i, df_cols]))
      n_impf  <- MOIf_locus[i, j] - n_obsf
      if (n_impf > 0) {
        fillf           <- df_cols[is.na(recodedf[i, df_cols])][1:n_impf]
        hiddenf[i, fillf] <- 1L
        possible        <- allele_frequencies$allele_codes[[j]]
        freqs           <- allele_frequencies$freq_matrix[j, 1:allele_frequencies$n_alleles[j]]
        recodedf[i, fillf] <- sample(possible, n_impf, replace = TRUE, prob = freqs)
      }
      hiddenf[i, df_cols[!is.na(recodedf[i, df_cols])]] <- 0L
    }
  }
  
  # Initialise parameters: sequencing error, allelic dropout, allele loss, and prior recrudescence probability
  q_mismatch  <- 0.01
  q_dropout  <- 0.05
  q_loss      <- 0.10
  prob_recrud <- 0.50
  classification <- ifelse(stats::runif(nids) < prob_recrud, 1, 0)
  
  num_records            <- floor((nruns - burnin) / record_interval)
  classification_history <- matrix(NA, nids, num_records)
  parameters_history     <- matrix(NA, 3,    num_records)
  # Row 1 = q_mismatch (sequencing error)
  # Row 2 = q_dropout (allelic dropout)
  # Row 3 = q_loss     (allele loss between timepoints)
  loglikelihood_history  <- rep(NA_real_, num_records)
  locus_lrs_history      <- array(NA, c(nids, nloci, num_records))
  locus_dists_history    <- array(NA, c(nids, nloci, num_records))
  alleles0_history <- if (record_hidden_alleles) array(NA, c(nids, maxMOI * nloci, num_records)) else NULL
  allelesf_history <- if (record_hidden_alleles) array(NA, c(nids, maxMOI * nloci, num_records)) else NULL
  
  # mcmc
  for (iter in 1:nruns) {
    locus_lrs_this_step <- matrix(NA_real_, nids, nloci)
    mindistance         <- matrix(NA_integer_, nids, nloci)
    
    for (i in 1:nids) {
      for (j in 1:nloci) {
        d0 <- recoded0[i, (maxMOI*(j-1)+1):(maxMOI*j)]
        d0 <- d0[!is.na(d0)]
        df <- recodedf[i, (maxMOI*(j-1)+1):(maxMOI*j)]
        df <- df[!is.na(df)]
        
        if (length(d0) == 0 || length(df) == 0 || !is_locus_comparable[i, j]) next
        
        n_lost  <- sum(!d0 %in% df)
        dropout <- q_loss ^ n_lost
        pair_sum <- 0
        
        for (rf in seq_along(df)) {
          
          idx <- match(df[rf], allele_frequencies$allele_codes[[j]])
          lam <- if (!is.na(idx)) max(allele_frequencies$freq_matrix[j, idx], 1e-6) else 1e-6
          
          # compute best match over all Day 0 alleles
          best_contrib <- 0
          
          for (r0 in seq_along(d0)) {
            if (d0[r0] == df[rf]) {
              contrib <- (1 - q_mismatch) / lam
            } else {
              contrib <- q_mismatch
            }
            best_contrib <- max(best_contrib, contrib)
          }
          
          pair_sum <- pair_sum + best_contrib
        }
        
        # accounting for extra alleles from recurrence
        locus_cols <- (maxMOI * (j - 1) + 1):(maxMOI * j)
        m_i        <- sum(c(hidden0[i, locus_cols], hiddenf[i, locus_cols]) == 1L, na.rm = TRUE)
        locus_lrs_this_step[i, j] <- dropout * (pair_sum / length(df)) * (q_dropout ^ m_i)
        
        mindistance[i, j] <- ifelse(any(df %in% d0), 0L, 1L)
      }
    }
    
    # Classification update
    likelihoodratio <- apply(locus_lrs_this_step, 1, function(r) {
      v <- r[!is.na(r)]
      if (length(v) == 0) return(1.0)
      prod(pmax(v, 1e-15))
    })
    
    prior_odds <- log(prob_recrud + 1e-10) - log(1 - prob_recrud + 1e-10)
    posterior_odds  <- prior_odds + log(likelihoodratio + 1e-10)
    prob_is_recrudescence <- 1 / (1 + exp(-posterior_odds))
    classification <- rbinom(nids, 1, prob_is_recrudescence)
    
    # Updating the hidden alleles
    for (i in 1:nids) {
      updated <- switch_hidden_ampseq(
        x = i, 
        hidden0 = hidden0, 
        hiddenf = hiddenf,
        recoded0 = recoded0, 
        recodedf = recodedf,
        classification = classification,
        q_mismatch = q_mismatch,
        q_loss =q_loss,
        allele_frequencies = allele_frequencies,
        nloci = nloci, 
        maxMOI = maxMOI
      )
      temp0 <- updated$recoded0
      tempf <- updated$recodedf
      temp0[is_observed0] <- original_recoded0[is_observed0]
      tempf[is_observedf] <- original_recodedf[is_observedf]
      recoded0 <- temp0
      recodedf <- tempf
    }
    
    #  Update q_loss: posterior probability of allele loss between Day 0 and recurrence
    n_lost_total <- 0
    n_ret_total  <- 0
    for (i in 1:nids) {
      w <- prob_is_recrudescence[i]
      if (w < 0.01) next
      for (j in 1:nloci) {
        if (!is_locus_comparable[i, j]) next
        d0u <- unique(recoded0[i, (maxMOI*(j-1)+1):(maxMOI*j)])
        d0u <- d0u[!is.na(d0u)]
        dfu <- unique(recodedf[i, (maxMOI*(j-1)+1):(maxMOI*j)])
        dfu <- dfu[!is.na(dfu)]
        if (length(d0u) == 0 || length(dfu) == 0) next
        n_lost_total <- n_lost_total + w * sum(!d0u %in% dfu)
        n_ret_total  <- n_ret_total  + w * sum( d0u %in% dfu)
      }
    }
    q_loss <- stats::rbeta(1, 1 + n_lost_total, 1 + n_ret_total)
    
    # Update q_mismatch: posterior sequencing error rate from allele match/mismatch pairs
    n_match_pairs    <- 0
    n_mismatch_pairs <- 0
    for (i in 1:nids) {
      w <- prob_is_recrudescence[i]
      if (w < 0.01) next 
      for (j in 1:nloci) {
        if (!is_locus_comparable[i, j]) next
        d0u <- unique(recoded0[i, (maxMOI*(j-1)+1):(maxMOI*j)])
        d0u <- d0u[!is.na(d0u)]
        dfu <- unique(recodedf[i, (maxMOI*(j-1)+1):(maxMOI*j)])
        dfu <- dfu[!is.na(dfu)]
        if (length(d0u) == 0 || length(dfu) == 0) next
        for (r0 in seq_along(d0u)) {
          for (rf in seq_along(dfu)) {
            if (d0u[r0] == dfu[rf]) {
              n_match_pairs    <- n_match_pairs    + w
            } else {
              n_mismatch_pairs <- n_mismatch_pairs + w
            }
          }
        }
      }
    }
    q_mismatch <- stats::rbeta(1, 1 + n_mismatch_pairs, 50 + n_match_pairs)
    
    # Update q_dropout: posterior allelic dropout rate from hidden vs observed allele counts
    n_hid <- sum(c(hidden0, hiddenf) == 1L, na.rm = TRUE)
    n_obs <- sum(c(hidden0, hiddenf) == 0L, na.rm = TRUE)
    q_dropout <- stats::rbeta(1, 1 + n_hid, 50 + n_obs)
    
    # update the allele frquencies
    reinfection_indices <- which(classification == 0)
    pool <- rbind(
      recoded0[reinfection_indices, , drop = FALSE],
      recodedf[reinfection_indices, , drop = FALSE],
      recoded_additional_site
    )
    for (j in 1:nloci) {
      new_row <- findposteriorfrequencies(
        locus_index        = j,
        tempdata           = pool,
        maxMOI             = maxMOI,
        allele_frequencies = allele_frequencies
      )
      n_all <- allele_frequencies$n_alleles[j]
      if (length(new_row) > 0 && !any(is.na(new_row)) && n_all > 0)
        allele_frequencies$freq_matrix[j, 1:n_all] <- new_row
    }
    
    # record the state variables
    if (iter > burnin && (iter - burnin) %% record_interval == 0) {
      ridx   <- (iter - burnin) / record_interval
      loglik <- sum(log(pmax(likelihoodratio, 1e-10)))
      
      classification_history[, ridx] <- classification
      parameters_history[1, ridx]    <- q_mismatch   # Row 1 = q_mismatch (sequencing error)
      parameters_history[2, ridx]    <- q_dropout    # Row 2 = q_dropout (allelic dropout)
      parameters_history[3, ridx]    <- q_loss       # Row 3 = q_loss     (allele loss between timepoints)
      loglikelihood_history[ridx]    <- loglik
      locus_lrs_history[,,  ridx]    <- locus_lrs_this_step
      locus_dists_history[,, ridx]   <- mindistance
      
      if (record_hidden_alleles) {
        alleles0_history[,, ridx] <- recoded0
        allelesf_history[,, ridx] <- recodedf
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