# =============================================================================
# R/stan_ampseq_interface.R
# Stan interface for amplicon sequencing data (ampseq)
# Parallel structure to stan_interface.R for length-polymorphic data
# =============================================================================

run_stan_sites_ampseq <- function(late_failures,
                                   additional,
                                   marker_info,
                                   mcmc_config,
                                   verbose = TRUE) {

  # --- 1. Extract config ---
  n_chains     <- as.integer(mcmc_config$n_chains)
  iter_total   <- as.integer(mcmc_config$iter)
  burn_in_frac <- as.numeric(mcmc_config$burn_in_frac)
  base_seed    <- as.integer(mcmc_config$random_seed)
  adapt_delta  <- as.numeric(
    if (!is.null(mcmc_config$adapt_delta)) mcmc_config$adapt_delta else 0.85
  )

  iter_warmup   <- max(200L, as.integer(floor(burn_in_frac * iter_total)))
  iter_sampling <- max(200L, iter_total - iter_warmup)

  # --- 2. Stan global setup ---
  # cmdstanr compiles to a native binary and caches it automatically.
  # parallel_chains runs each chain as a separate OS process — reliable on Windows.

  stan_file <- system.file("stan", "malrebay_ampseq_model.stan", package = "MalReBay")
  if (!file.exists(stan_file)) {
    stan_file <- file.path("inst", "stan", "malrebay_ampseq_model.stan")
  }
  if (!file.exists(stan_file)) stop("Cannot find malrebay_ampseq_model.stan", call. = FALSE)

  if (verbose) message("Compiling Stan ampseq model...")
  stan_model_obj <- cmdstanr::cmdstan_model(stan_file)

  # --- 3. Output containers ---
  site_names          <- as.character(unique(late_failures$Site))
  out_classifications <- list()
  out_loglikelihood   <- list()
  out_ids             <- list()
  out_locus_summary   <- list()
  out_locus_lrs       <- list()
  out_locus_dists     <- list()
  out_locinames       <- list()
  out_stan_fits       <- list()

  # --- 4. Site loop ---
  for (site in site_names) {
    if (verbose) message("\n--- Processing Site: ", site, " ---")

    late_site <- late_failures[late_failures$Site == site, ]
    add_site  <- additional[additional$Site == site, ]

    late_site <- late_site[, colnames(late_site) != "Site", drop = FALSE]
    add_site  <- add_site[,  colnames(add_site)  != "Site", drop = FALSE]

    ids <- unique(gsub(" Day 0", "",
                       late_site$Sample.ID[grepl("Day 0", late_site$Sample.ID)]))
    if (length(ids) == 0) {
      if (verbose) message("  No patient pairs found. Skipping.")
      next
    }

    # Derive loci and maxMOI from column names
    marker_cols <- grep("(_allele_\\d+|_\\d+)$", colnames(late_site), value = TRUE)
    if (length(marker_cols) == 0) {
      if (verbose) message("  No allele columns found. Skipping.")
      next
    }
    maxMOI    <- max(as.integer(gsub(".*(_allele_|_)(\\d+)$", "\\2", marker_cols)),
                     na.rm = TRUE)
    locinames <- unique(gsub("(_allele_\\d+|_\\d+)$", "", marker_cols))
    nloci     <- length(locinames)

    # Locus comparability matrix
    locus_summary <- data.frame(
      patient_id        = ids,
      n_available_d0    = 0L,
      n_available_df    = 0L,
      n_comparable_loci = 0L
    )
    is_locus_comparable <- matrix(
      FALSE, nrow = length(ids), ncol = nloci,
      dimnames = list(ids, locinames)
    )

    for (i in seq_along(ids)) {
      pid    <- ids[i]
      d0_row <- late_site[grepl(paste0("\\b", pid, " Day 0\\b"),      late_site$Sample.ID), ]
      df_row <- late_site[grepl(paste0("\\b", pid, " recurrence\\b"), late_site$Sample.ID), ]
      if (nrow(d0_row) == 0 || nrow(df_row) == 0) next

      for (ln in locinames) {
        lc <- grep(paste0("^", ln, "_"), colnames(late_site), value = TRUE)
        has_d0 <- any(!is.na(d0_row[, lc]) & d0_row[, lc] != "NA")
        has_df <- any(!is.na(df_row[, lc]) & df_row[, lc] != "NA")
        if (has_d0) locus_summary$n_available_d0[i] <- locus_summary$n_available_d0[i] + 1L
        if (has_df) locus_summary$n_available_df[i] <- locus_summary$n_available_df[i] + 1L
        if (has_d0 && has_df) {
          locus_summary$n_comparable_loci[i]    <- locus_summary$n_comparable_loci[i] + 1L
          is_locus_comparable[pid, ln]           <- TRUE
        }
      }
    }

    # Prepare Stan data
    if (verbose) message("  Preparing Stan ampseq data...")
    sd <- prepare_stan_data_ampseq(
      late_failures_site  = late_site,
      additional_site     = add_site,
      marker_info         = marker_info,
      ids                 = ids,
      locinames           = locinames,
      maxMOI              = maxMOI,
      is_locus_comparable = is_locus_comparable
    )

    if (!validate_stan_data_ampseq(sd)) next

    # Init function: flat frequencies, low error/loss rates.
    # Closure over sd so freq dimensions are correct.
    # Note: q_dropout included to match all three parameters declared in the model.
    init_fun <- local({
      sd_local <- sd
      function() {
        list(
          q_mismatch = 0.01,   # matches old R initial value
          q_loss     = 0.10,   # matches old R initial value
          q_dropout  = 0.05,   # matches old R initial value
          freq       = lapply(seq_len(sd_local$J), function(j) {
            x <- rep(0.1 / sd_local$max_K, sd_local$max_K)
            x[1:sd_local$K[j]] <- 1.0 / sd_local$K[j]
            x / sum(x)  # ensure valid simplex
          })
        )
      }
    })

    # Run Stan
    if (verbose) message("  Running Stan (", n_chains, " chains, ",
                         iter_sampling, " samples)...")

    fit <- tryCatch({
      stan_model_obj$sample(
        data            = stan_data_only(sd),
        chains          = n_chains,
        parallel_chains = min(n_chains, parallel::detectCores(logical = FALSE)),
        iter_warmup     = iter_warmup,
        iter_sampling   = iter_sampling,
        init            = init_fun,
        adapt_delta     = adapt_delta,
        max_treedepth   = 12L,
        seed            = base_seed,
        refresh         = if (verbose) 100L else 0L,
        output_dir      = tempdir()
      )
    }, error = function(e) {
      warning("Stan sampling failed for site '", site, "': ", e$message)
      NULL
    })

    if (is.null(fit)) next
    
    # check all chains actually completed
    if (all(fit$return_codes() != 0)) {
      message(paste(fit$output(), collapse = "\n"))
      warning("All chains failed for site '", site, "'. Skipping.")
      next
    }

    # Extract results
    extracted <- extract_stan_results_ampseq(fit, ids, locinames, nloci, length(ids))

    out_classifications[[site]] <- extracted$p_recrud_draws
    out_loglikelihood[[site]]   <- extracted$loglik_chains
    out_ids[[site]]             <- ids
    out_locus_summary[[site]]   <- locus_summary
    out_locinames[[site]]       <- locinames
    out_locus_lrs[[site]]       <- extracted$locus_lrs
    out_locus_dists[[site]]     <- extracted$locus_dists
    out_stan_fits[[site]]       <- fit
  }

  return(list(
    classifications          = out_classifications,
    all_chains_loglikelihood = out_loglikelihood,
    ids                      = out_ids,
    locus_summary            = out_locus_summary,
    locus_lrs                = out_locus_lrs,
    locus_dists              = out_locus_dists,
    locinames                = out_locinames,
    stan_fits                = out_stan_fits
  ))
}

extract_stan_results_ampseq <- function(fit, ids, locinames, nloci, nids) {

  # draws_matrix: rows = draws (all chains merged), cols = p_recrud[1..N]
  p_recrud_draws <- fit$draws("p_recrud", format = "draws_matrix")

  # draws_array: [iterations, chains, variables] — preserves chain structure
  lp_array  <- fit$draws("lp__", format = "draws_array")
  n_iter    <- dim(lp_array)[1]
  n_chains  <- dim(lp_array)[2]
  lp_matrix <- matrix(as.numeric(lp_array), nrow = n_iter, ncol = n_chains)
  loglik_chains <- lapply(seq_len(n_chains), function(ch) lp_matrix[, ch])

  n_draws     <- nrow(p_recrud_draws)
  # locus_lrs and locus_dists not directly available from Stan output for ampseq.
  # They would require adding generated quantities per locus — left as NA for now.
  locus_lrs   <- array(NA_real_, dim = c(nids, nloci, n_draws))
  locus_dists <- array(NA_real_, dim = c(nids, nloci, n_draws))

  return(list(
    p_recrud_draws = p_recrud_draws,
    loglik_chains  = loglik_chains,
    locus_lrs      = locus_lrs,
    locus_dists    = locus_dists
  ))
}
