#' Run the Full MCMC Analysis Pipeline Across All Sites
#'
#' @description
#' Manages MCMC runs across all geographical sites and returns combined results.
#' @param late_failures Data frame of late treatment failure samples.
#' @param additional Additional neutral marker data.
#' @param marker_info Marker definitions used for allele matching.
#' @param mcmc_config List of MCMC configuration options.
#' @param data_type Either "ampseq" or "length".
#' @param output_folder Directory where results will be saved.
#' @param verbose Logical; whether to print progress messages.
#'   and diagnostics to the console.
#'
#' @keywords internal
#' @noRd
#'
run_all_sites <- function(
    late_failures,
    additional,
    marker_info,
    mcmc_config,
    data_type,
    output_folder = NULL,
    verbose = TRUE) {

  # Extract mcmc configuration parameters
  n_chains              <- as.integer(mcmc_config$n_chains)
  rhat_threshold        <- as.numeric(mcmc_config$rhat_threshold)
  ess_threshold         <- as.numeric(mcmc_config$ess_threshold)
  chunk_size            <- as.integer(mcmc_config$chunk_size)
  max_iterations        <- as.integer(mcmc_config$max_iterations)
  burn_in_frac          <- as.numeric(
    if (!is.null(mcmc_config$burn_in_frac)) mcmc_config$burn_in_frac else 0.5
  )
  record_hidden_alleles <- as.logical(
    if (!is.null(mcmc_config$record_hidden_alleles)) mcmc_config$record_hidden_alleles else FALSE
  )
  base_seed             <- as.integer(
    if (!is.null(mcmc_config$random_seed)) mcmc_config$random_seed else 42
  )
  ess_tail_threshold  <- as.numeric(
    if (!is.null(mcmc_config$ess_tail_threshold)) mcmc_config$ess_tail_threshold else 400
  )
  geweke_threshold    <- as.numeric(
    if (!is.null(mcmc_config$geweke_threshold)) mcmc_config$geweke_threshold else 1.96
  )

  # Validate config 
  if (is.na(n_chains)       || n_chains < 1)       stop("mcmc_config$n_chains must be a positive integer.",        call. = FALSE)
  if (is.na(rhat_threshold) || rhat_threshold <= 1) stop("mcmc_config$rhat_threshold must be a number above 1.",   call. = FALSE)
  if (is.na(ess_threshold)  || ess_threshold < 1)   stop("mcmc_config$ess_threshold must be a positive number.",   call. = FALSE)
  if (is.na(chunk_size)     || chunk_size < 1)       stop("mcmc_config$chunk_size must be a positive integer.",     call. = FALSE)
  if (is.na(max_iterations) || max_iterations < 1)   stop("mcmc_config$max_iterations must be a positive integer.", call. = FALSE)
  if (is.na(burn_in_frac)   || burn_in_frac <= 0 || burn_in_frac >= 1) stop("mcmc_config$burn_in_frac must be between 0 and 1 exclusive.", call. = FALSE)
  if (is.na(record_hidden_alleles)) stop("mcmc_config$record_hidden_alleles must be TRUE or FALSE.", call. = FALSE)
  if (is.na(base_seed)      || base_seed < 0) stop("mcmc_config$random_seed must be a non-negative integer.", call. = FALSE)
  if (is.na(ess_tail_threshold) || ess_tail_threshold < 1) stop("mcmc_config$ess_tail_threshold must be a positive number.", call. = FALSE)
  if (is.na(geweke_threshold)   || geweke_threshold <= 0)  stop("mcmc_config$geweke_threshold must be a positive number.",  call. = FALSE)

  # Initialise site-level output lists
  site_names                 <- as.character(unique(late_failures$Site))
  local_sites_classification <- list()
  local_sites_alleles0       <- list()
  local_sites_allelesf       <- list()
  local_sites_ids            <- list()
  local_sites_locus_summary  <- list()
  local_sites_locus_lrs      <- list()
  local_sites_locus_dists    <- list()
  local_sites_locinames      <- list()
  local_sites_loglikelihood  <- list()

  for (site in site_names) {
    late_failures_site <- late_failures[late_failures$Site == site, ]
    additional_site    <- additional[additional$Site == site, ]
    ids  <- unique(gsub(" Day 0", "", late_failures_site$Sample.ID[grepl("Day 0", late_failures_site$Sample.ID)]))
    nids <- length(ids)

    if (verbose) message("Number of patient pairs (nids) found for site '", site, "': ", nids)
    if (nids == 0) {
      if (verbose) message("INFO: No valid Day 0/Failure pairs found for site '", site, "'. Skipping.")
      next
    }

    mcmc_engine_function <- NULL
    engine_args          <- NULL

    if (data_type == "ampseq") {
      if (verbose) message("INFO: Preparing AMPSEQ MCMC engine for site: ", site)
      marker_cols <- grep("(_allele_\\d+|_\\d+)$", colnames(late_failures_site), value = TRUE)
      maxMOI      <- if (length(marker_cols) > 0) max(as.numeric(gsub(".*(_allele_|_)(\\d+)$", "\\2", marker_cols)), na.rm = TRUE) else 0
      locinames   <- unique(gsub("(_allele_\\d+|_\\d+)$", "", marker_cols))
      local_sites_locinames[[site]] <- locinames
      nloci <- length(locinames)
      mcmc_engine_function <- run_one_chain_ampseq
      engine_args <- list(
        nids               = nids, ids = ids, nloci = nloci, maxMOI = maxMOI, locinames = locinames,
        late_failures_site = late_failures_site[, -which(names(late_failures_site) == "Site")],
        additional_site    = additional_site[,    -which(names(additional_site)    == "Site")],
        marker_info        = marker_info,
        record_hidden_alleles = record_hidden_alleles
      )
    } else {
      if (verbose) message("INFO: Preparing LENGTH-POLYMORPHIC MCMC engine for site: ", site)
      allele_definitions <- suppressMessages(define_alleles(
        genotypedata = rbind(late_failures_site, additional_site),
        marker_info  = marker_info
      ))
      locinames <- names(allele_definitions)
      local_sites_locinames[[site]] <- locinames
      nloci       <- length(locinames)
      marker_cols <- grep("_\\d+$", colnames(late_failures_site), value = TRUE)
      maxMOI      <- if (length(marker_cols) > 0) max(as.numeric(gsub(".*_", "", marker_cols)), na.rm = TRUE) else 0
      mcmc_engine_function <- run_one_chain
      engine_args <- list(
        nids               = nids, ids = ids, nloci = nloci, maxMOI = maxMOI, locinames = locinames,
        late_failures_site = late_failures_site[, -which(names(late_failures_site) == "Site")],
        additional_site    = additional_site[,    -which(names(additional_site)    == "Site")],
        allele_definitions = allele_definitions,
        marker_info        = marker_info,
        record_hidden_alleles = record_hidden_alleles
      )
    }

    if (nloci == 0) {
      if (verbose) message("INFO: No valid loci with data found for site '", site, "'. Skipping.")
      next
    }

    locus_summary       <- data.frame(patient_id = ids, n_available_d0 = 0, n_available_df = 0, n_comparable_loci = 0)
    is_locus_comparable <- matrix(FALSE, nrow = nids, ncol = nloci, dimnames = list(ids, locinames))

    for (i in 1:nrow(locus_summary)) {
      patient_id <- locus_summary$patient_id[i]
      pattern_d0 <- paste0("\\b", patient_id, " Day 0\\b")
      pattern_df <- paste0("\\b", patient_id, " recurrence\\b")
      d0_row <- late_failures_site[grepl(pattern_d0, late_failures_site$Sample.ID), ]
      df_row <- late_failures_site[grepl(pattern_df, late_failures_site$Sample.ID), ]
      if (nrow(d0_row) == 0 || nrow(df_row) == 0) next
      for (locus_name in locinames) {
        locus_cols    <- grep(paste0("^", locus_name, "_"), colnames(late_failures_site), value = TRUE)
        has_d0_allele <- any(!is.na(d0_row[, locus_cols]))
        has_df_allele <- any(!is.na(df_row[, locus_cols]))
        if (has_d0_allele) locus_summary$n_available_d0[i] <- locus_summary$n_available_d0[i] + 1
        if (has_df_allele) locus_summary$n_available_df[i] <- locus_summary$n_available_df[i] + 1
        if (has_d0_allele && has_df_allele) {
          locus_summary$n_comparable_loci[i]         <- locus_summary$n_comparable_loci[i] + 1
          is_locus_comparable[patient_id, locus_name] <- TRUE
        }
      }
    }
    local_sites_locus_summary[[site]] <- locus_summary
    engine_args$is_locus_comparable   <- is_locus_comparable

    # MCMC loop with convergence stop rule and post-burn-in sample collection
    full_loglik_history <- list()
    full_chain_results  <- list()
    total_iterations    <- 0L
    converged           <- FALSE

    while (!converged && total_iterations <= max_iterations) {
      total_iterations <- total_iterations + chunk_size

      chunk_results <- foreach::foreach(
        id             = 1:n_chains,
        .packages      = c("abind", "coda", "MalReBay"),
        .errorhandling = "stop"
      ) %dopar% {
        set.seed(base_seed + id)
        full_args_for_chunk <- c(
          list(chain_id = id, nruns = chunk_size,
               burnin = 0, record_interval = 10),
          engine_args
        )
        do.call(mcmc_engine_function, full_args_for_chunk)
      }

      if (length(full_chain_results) == 0) {
        full_chain_results  <- chunk_results
        full_loglik_history <- lapply(chunk_results, `[[`, "state_loglikelihood")
      } else {
        for (i in 1:n_chains) {
          full_chain_results[[i]]$state_parameters     <- cbind(full_chain_results[[i]]$state_parameters,     chunk_results[[i]]$state_parameters)
          full_chain_results[[i]]$state_classification <- cbind(full_chain_results[[i]]$state_classification, chunk_results[[i]]$state_classification)
          full_loglik_history[[i]]                     <- c(full_loglik_history[[i]], chunk_results[[i]]$state_loglikelihood)
          full_chain_results[[i]]$locus_lrs            <- abind::abind(full_chain_results[[i]]$locus_lrs,   chunk_results[[i]]$locus_lrs,   along = 3)
          full_chain_results[[i]]$locus_dists          <- abind::abind(full_chain_results[[i]]$locus_dists, chunk_results[[i]]$locus_dists, along = 3)
          if (record_hidden_alleles) {
            full_chain_results[[i]]$state_recoded0 <- abind::abind(full_chain_results[[i]]$state_recoded0, chunk_results[[i]]$state_recoded0, along = 3)
            full_chain_results[[i]]$state_recodedf <- abind::abind(full_chain_results[[i]]$state_recodedf, chunk_results[[i]]$state_recodedf, along = 3)
          }
        }
      }

      if (length(full_loglik_history) == 0 || length(full_loglik_history[[1]]) == 0) {
        if (verbose) message("INFO: Log-likelihood history not yet populated. Running another chunk...")
        next
      }

      # Build post-burn-in draws matrix: rows = iterations, cols = chains
      # Used by both posterior and coda diagnostics below
      post_burn_draws <- tryCatch({
        do.call(cbind, lapply(full_loglik_history, function(x) {
          start_idx <- floor(burn_in_frac * length(x)) + 1
          x[start_idx:length(x)]
        }))
      }, error = function(e) NULL)

      if (is.null(post_burn_draws) || nrow(post_burn_draws) < 50) next

      # When chains stabilise completely (near-constant trace), classical
      # Gelman-Rubin is undefined (W (within chain) -> 0). We check the chain stability 
      # before applying variance-based diagnostics.
      chain_sds <- apply(post_burn_draws, 2, sd, na.rm = TRUE)
      if (all(!is.na(chain_sds) & chain_sds < 0.1)) {
        converged <- TRUE
        if (verbose) message("INFO: Convergence reached for site '", site,
                             "' -- chains are stable (SD < 0.1) at iteration ",
                             total_iterations, ".")
        next
      }

      # Rank-normalised R-hat + bulk and tail ESS
      rhat_rank <- tryCatch(
        posterior::rhat(post_burn_draws),
        error = function(e) NA_real_
      )
      ess_bulk <- tryCatch(
        posterior::ess_bulk(post_burn_draws),
        error = function(e) NA_real_
      )
      ess_tail <- tryCatch(
        posterior::ess_tail(post_burn_draws),
        error = function(e) NA_real_
      )

      rhat_ok     <- !is.na(rhat_rank) && rhat_rank  < rhat_threshold
      ess_b_ok    <- !is.na(ess_bulk)  && ess_bulk   > ess_threshold
      ess_t_ok    <- !is.na(ess_tail)  && ess_tail   > ess_tail_threshold

      if (rhat_ok && ess_b_ok && ess_t_ok) {
        converged <- TRUE
        if (verbose) message("INFO: Convergence reached for site '", site,
                             "' at iteration ", total_iterations,
                             " (Rhat_rank=", round(rhat_rank, 4),
                             ", ESS_bulk=",  round(ess_bulk, 0),
                             ", ESS_tail=",  round(ess_tail, 0), ").")
        next
      }

      # Classical Gelman-Rubin + ESS fallback 
      mcmc_list_loglik <- tryCatch(
        coda::mcmc.list(lapply(full_loglik_history, function(x) {
          start_idx <- floor(burn_in_frac * length(x)) + 1
          coda::mcmc(x[start_idx:length(x)])
        })),
        error = function(e) NULL
      )

      if (!is.null(mcmc_list_loglik)) {
        r_hat_classic <- try(
          coda::gelman.diag(mcmc_list_loglik, autoburnin = FALSE)$psrf[1, 1],
          silent = TRUE
        )
        ess_classic <- try(
          coda::effectiveSize(mcmc_list_loglik)[1],
          silent = TRUE
        )
        if (!inherits(r_hat_classic, "try-error") &&
            !inherits(ess_classic,   "try-error") &&
            !is.na(r_hat_classic) && r_hat_classic < rhat_threshold &&
            !is.na(ess_classic)   && ess_classic   > ess_threshold) {
          converged <- TRUE
          if (verbose) message("INFO: Convergence reached for site '", site,
                               "' at iteration ", total_iterations,
                               " (classical Rhat=", round(r_hat_classic, 4),
                               ", ESS=", round(ess_classic, 0), ").")
          next
        }
      }

      # Geweke diagnostic
      if (!is.null(mcmc_list_loglik)) {
        geweke_vals <- tryCatch(
          sapply(mcmc_list_loglik, function(ch) coda::geweke.diag(ch)$z),
          error = function(e) NULL
        )
        if (!is.null(geweke_vals) &&
            all(!is.na(geweke_vals)) &&
            all(abs(geweke_vals) < geweke_threshold)) {
          converged <- TRUE
          if (verbose) message("INFO: Convergence reached for site '", site,
                               "' via Geweke diagnostic at iteration ",
                               total_iterations,
                               " (max |Z|=", round(max(abs(geweke_vals)), 4), ").")
        }
      }
    }

    if (!converged && verbose) {
      message("WARNING: Site '", site, "' did not converge within ",
              max_iterations, " iterations. Results may be unreliable.")
    }

    # Collect post-burn-in samples
    if (length(full_chain_results) > 0) {
      num_total_samples_per_chain <- ncol(full_chain_results[[1]]$state_parameters)
      burn_in_samples_per_chain   <- floor(burn_in_frac * num_total_samples_per_chain)

      if (burn_in_samples_per_chain < num_total_samples_per_chain) {
        keep_indices <- (burn_in_samples_per_chain + 1):num_total_samples_per_chain

        local_sites_classification[[site]] <- do.call(cbind, lapply(full_chain_results, function(x) x$state_classification[, keep_indices, drop = FALSE]))
        local_sites_locus_lrs[[site]]      <- do.call(abind::abind, list(lapply(full_chain_results, function(x) x$locus_lrs[,,   keep_indices, drop = FALSE]), along = 3))
        local_sites_locus_dists[[site]]    <- do.call(abind::abind, list(lapply(full_chain_results, function(x) x$locus_dists[,, keep_indices, drop = FALSE]), along = 3))
        local_sites_loglikelihood[[site]]  <- lapply(full_loglik_history, function(x) x[keep_indices])
        local_sites_ids[[site]]            <- ids
        local_sites_locus_summary[[site]]  <- locus_summary
        local_sites_locinames[[site]]      <- locinames
      } else {
        if (verbose) message("INFO: Burn-in equals or exceeds total samples for site '",
                             site, "'. No post-burn-in samples to summarise.")
      }
    }
  }

  return(list(
    classifications          = local_sites_classification,
    all_chains_loglikelihood = local_sites_loglikelihood,
    ids                      = local_sites_ids,
    locus_summary            = local_sites_locus_summary,
    locus_lrs                = local_sites_locus_lrs,
    locus_dists              = local_sites_locus_dists,
    locinames                = local_sites_locinames
  ))
}