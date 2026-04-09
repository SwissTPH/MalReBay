# R/stan_data_prep.R
# Step 1 of Stan transition — pure R, no Stan dependency
#
# Purpose:
#   Packages allele_utils.R output into a named list that rstan::stan()
#   accepts as its `data` argument.
#
# Depends on:  R/allele_utils.R  (unchanged)
# Called by:   R/stan_interface.R  (Step 3)
# Replaces:    Initialisation block of run_one_chain() in
#              R/mcmc_length_polymorphic.R  (retired after validation)
#
# Debugging workflow:
#   sd <- prepare_stan_data(...)
#   validate_stan_data(sd)   # must print ALL CHECKS PASSED before running Stan
# prepare_stan_data()
#
# Arguments:
#   late_failures_site   data.frame — site-filtered rows, Site column removed
#   additional_site      data.frame — site-filtered additional data
#   allele_definitions   list       — output of define_alleles()
#   marker_info          data.frame — marker metadata, pre-filtered to site loci
#   ids                  character  — patient ID vector (length nids)
#   locinames            character  — locus name vector (length nloci), ordered
#   maxMOI               integer    — maximum multiplicity of infection
#   is_locus_comparable  matrix     — logical [nids × nloci]
#
# Returns:
#   Named list ready for rstan::stan(data = ...).
#   Also carries $allele_frequencies for use in extract_results().
#' Prepare Stan data for the length-polymorphic model
#'
#' Converts site-filtered genotype data and allele definitions into a named
#' list accepted by \code{rstan::sampling(data = ...)}.
#'
#' @param late_failures_site   Data frame of late-failure rows for one site
#'   (Site column removed).
#' @param additional_site      Data frame of additional background rows for one
#'   site (Site column removed).
#' @param allele_definitions   List of allele bin definition data frames, one
#'   per locus (output of \code{define_alleles}).
#' @param marker_info          Data frame of marker metadata pre-filtered to
#'   the site's loci.
#' @param ids                  Character vector of patient IDs.
#' @param locinames            Character vector of locus names.
#' @param maxMOI               Integer maximum multiplicity of infection.
#' @param is_locus_comparable  Logical matrix \code{[nids x nloci]} indicating
#'   whether both day-0 and recurrence data are available for each
#'   patient/locus combination.
#'
#' @return A named list ready for \code{rstan::sampling}. Also carries
#'   \code{$allele_frequencies} for use in result extraction.
#' @noRd
prepare_stan_data <- function(late_failures_site,
                              additional_site,
                              allele_definitions,
                              marker_info,
                              ids,
                              locinames,
                              maxMOI,
                              is_locus_comparable) {
  
  nids  <- length(ids)
  nloci <- length(locinames)
  
  # ---- 1. MOI per patient ----
  MOI0 <- rep(1L, nids); MOIf <- rep(1L, nids)
  for (i in seq_len(nids)) {
    for (j in seq_len(nloci)) {
      loci_cols <- grepl(paste0(locinames[j], "_"), colnames(late_failures_site))
      n0 <- sum(!is.na(late_failures_site[grepl(paste0("\\b", ids[i], " Day 0\\b"), late_failures_site$Sample.ID), loci_cols]))
      nf <- sum(!is.na(late_failures_site[grepl(paste0("\\b", ids[i], " recurrence\\b"), late_failures_site$Sample.ID), loci_cols]))
      MOI0[i] <- max(MOI0[i], n0); MOIf[i] <- max(MOIf[i], nf)
    }
  }

  # ---- 2. Site-Specific Allele Pruning & Recoding ----
  recoded0 <- matrix(0L, nids, maxMOI * nloci)
  recodedf <- matrix(0L, nids, maxMOI * nloci)
  
  K_site <- integer(nloci)
  pruned_defs <- list()
  
  # Temporary storage for counts (will be trimmed to max_K later)
  tmp_counts <- matrix(0L, nloci, 200) 

  for (j in seq_len(nloci)) {
    locus <- locinames[j]
    cols  <- grepl(paste0("^", locus, "_"), colnames(late_failures_site))
    
    obs_trial <- as.matrix(late_failures_site[, cols, drop=FALSE])
    obs_add   <- if(nrow(additional_site) > 0) as.matrix(additional_site[, cols, drop=FALSE]) else matrix(nrow=0, ncol=0)
    
    all_obs_raw <- c(as.vector(obs_trial), as.vector(obs_add))
    all_obs_ids <- unique(sapply(all_obs_raw, function(x) recodeallele(allele_definitions[[j]], x)))
    all_obs_ids <- sort(all_obs_ids[!is.na(all_obs_ids) & all_obs_ids > 0])
    
    K_site[j] <- length(all_obs_ids)
    if(K_site[j] == 0) {
       K_site[j] <- 1
       all_obs_ids <- 1 
    }
    
    id_map <- setNames(seq_along(all_obs_ids), all_obs_ids)
    pruned_defs[[j]] <- allele_definitions[[j]][all_obs_ids, , drop=FALSE]
    
    day0_rows <- grepl("Day 0",      late_failures_site$Sample.ID)
    dayf_rows <- grepl("recurrence", late_failures_site$Sample.ID)

    # Build a named lookup once for O(1) per-query patient-index resolution.
    pid_idx <- setNames(seq_along(ids), ids)

    for(row_idx in seq_len(nrow(late_failures_site))){
      p_id  <- gsub(" Day 0| recurrence", "", late_failures_site$Sample.ID[row_idx])
      p_idx <- pid_idx[p_id]
      if(is.na(p_idx)) next

      for(c_idx in seq_len(ncol(obs_trial))){
        old_id <- recodeallele(allele_definitions[[j]], obs_trial[row_idx, c_idx])
        if(!is.na(old_id) && old_id > 0){
          new_id <- id_map[as.character(old_id)]
          slot <- (maxMOI * (j - 1) + c_idx)
          if(day0_rows[row_idx]) recoded0[p_idx, slot] <- as.integer(new_id)
          if(dayf_rows[row_idx]) recodedf[p_idx, slot] <- as.integer(new_id)
        }
      }
    }
    
    # Count all observed alleles (trial + additional) into tmp_counts so that 
    # the Dirichlet prior on allele frequencies is informed by the observed data.
    
    # Count trial data alleles
    if(nrow(obs_trial) > 0){
      trial_ids_raw <- sapply(as.vector(obs_trial), function(x) recodeallele(allele_definitions[[j]], x))
      trial_ids_raw <- trial_ids_raw[!is.na(trial_ids_raw) & trial_ids_raw > 0]
      if(length(trial_ids_raw) > 0){
        new_trial_ids <- id_map[as.character(trial_ids_raw)]
        cnts <- table(factor(new_trial_ids, levels = seq_along(all_obs_ids)))
        tmp_counts[j, seq_along(all_obs_ids)] <-
          tmp_counts[j, seq_along(all_obs_ids)] + as.integer(cnts)
      }
    }
    
    # Count additional_site alleles (as before)
    if(nrow(obs_add) > 0){
      add_ids_raw <- sapply(obs_add, function(x) recodeallele(allele_definitions[[j]], x))
      add_ids_raw <- add_ids_raw[!is.na(add_ids_raw) & add_ids_raw > 0]
      if(length(add_ids_raw) > 0){
        new_add_ids <- id_map[as.character(add_ids_raw)]
        cnts <- table(factor(new_add_ids, levels=seq_along(all_obs_ids)))
        tmp_counts[j, seq_along(all_obs_ids)] <-
          tmp_counts[j, seq_along(all_obs_ids)] + as.integer(cnts)
      }
    }
  }

  max_K <- max(K_site)
  additional_counts <- tmp_counts[, 1:max_K, drop=FALSE]

  # Pruned Distance Array
  dist_array <- array(0, dim = c(nloci, max_K, max_K))
  for (j in seq_len(nloci)) {
    if (K_site[j] > 0) {
      d_mat <- as.matrix(stats::dist(rowMeans(pruned_defs[[j]])))
      dist_array[j, 1:K_site[j], 1:K_site[j]] <- d_mat
    }
  }
  
  # Method & Threshold logic
  method_int <- sapply(locinames, function(ln) {
    m <- marker_info$binning_method[match(ln, marker_info$marker_id)]
    if(is.na(m)) return(3L)
    switch(as.character(m), "microsatellite" = 1L, "msp_glurp" = 2L, 3L)
  })
  
  threshold_vec <- sapply(seq_len(nloci), function(j) {
    if (method_int[j] == 2L) {
      val <- marker_info$repeatlength[match(locinames[j], marker_info$marker_id)]
      return(if(is.na(val)) 0 else val)
    } else return(0)
  })

  # Derive max_dist from the largest allele distance in the data, 
  # ensuring the distance vector is no larger than necessary.
  computed_max_dist <- max(1L, ceiling(max(dist_array)))

  return(list(
    N        = nids,
    J        = nloci,
    maxMOI   = as.integer(maxMOI),
    max_K    = as.integer(max_K),
    max_dist = as.integer(computed_max_dist),
    K        = as.array(K_site),
    recoded0 = recoded0,
    recodedf = recodedf,
    hidden0  = matrix(as.integer(recoded0 == 0), nids, nloci * maxMOI),
    hiddenf  = matrix(as.integer(recodedf == 0), nids, nloci * maxMOI),
    MOI0     = as.array(MOI0),
    MOIf     = as.array(MOIf),
    dist_array = dist_array,
    method_int = as.array(method_int),
    threshold  = as.array(threshold_vec),
    comparable = matrix(as.integer(is_locus_comparable), nids, nloci),
    additional_counts = additional_counts
  ))
}

#' Validate Stan data list for the length-polymorphic model
#'
#' Runs a series of checks on the list returned by \code{prepare_stan_data}
#' and prints \code{[PASS]} / \code{[FAIL]} for each. Returns \code{TRUE}
#' only if all checks pass.
#'
#' @param sd Named list produced by \code{prepare_stan_data}.
#' @return Logical scalar.
#' @noRd
validate_stan_data <- function(sd) {
  ok  <- TRUE
  chk <- function(label, test) {
    if (!isTRUE(test)) {
      message("  [FAIL] ", label)
      ok <<- FALSE
    } else {
      message("  [PASS] ", label)
    }
  }
  
  message("=== validate_stan_data ===")
  chk("N >= 1",       sd$N >= 1)
  chk("J >= 1",       sd$J >= 1)
  chk("max_K >= 1",   sd$max_K >= 1)
  chk("length(K) == J", length(sd$K) == sd$J)
  chk("dim(recoded0) == [N, J*maxMOI]", all(dim(sd$recoded0) == c(sd$N, sd$J * sd$maxMOI)))
  chk("additional_counts exists", !is.null(sd$additional_counts))
  chk("additional_counts has data", sum(sd$additional_counts) > 0)
  
  message(if (ok) "=== ALL CHECKS PASSED ===" else "=== VALIDATION FAILED ===")
  return(ok)
}

#' Strip internal bookkeeping fields from a Stan data list
#'
#' Removes any element whose name starts with \code{.} (used internally for
#' diagnostics) before passing the list to \code{rstan::sampling}.
#'
#' @param sd Named list produced by \code{prepare_stan_data} or
#'   \code{prepare_stan_data_ampseq}.
#' @return The same list with private fields removed.
#' @noRd
stan_data_only <- function(sd) {
  sd[!startsWith(names(sd), ".")]
}
