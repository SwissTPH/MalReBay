# =============================================================================
# R/stan_interface.R
# Single Stan interface for both marker families (length-polymorphic and
# amplicon sequencing / ampseq). Per-locus dispatch on
# marker_info$binning_method happens inside prepare_stan_data(); this file
# just runs the (single) compiled model against whatever prepare_stan_data()
# produces.
# =============================================================================

run_stan_sites <- function(late_failures,
                           additional,
                           marker_info,
                           mcmc_config,
                           verbose = TRUE) {

  # Extract MCMC configuration
  n_chains     <- as.integer(mcmc_config$n_chains)
  iter_total   <- as.integer(mcmc_config$iter)
  burn_in_frac <- as.numeric(mcmc_config$burn_in_frac)
  base_seed    <- as.integer(mcmc_config$random_seed)
  adapt_delta  <- as.numeric(if(!is.null(mcmc_config$adapt_delta)) mcmc_config$adapt_delta else 0.85)
  iter_warmup   <- max(200L, as.integer(floor(burn_in_frac * iter_total)))
  iter_sampling <- max(200L, iter_total - iter_warmup)

  # 2. Stan Global Setup
  # The model is pre-compiled at package install time via {instantiate}
  # (src/stan/malrebay_model.stan,src/install.libs.R). 
  # This just loads the shipped executable.
  # parallel_chains uses separate OS processes -- reliable on operating system.
  stan_model_obj <- instantiate::stan_package_model(
    name    = "malrebay_model",
    package = "MalReBay"
  )

  # Initialize Containers
  site_names          <- as.character(unique(late_failures$Site))
  out_classifications <- list()
  out_loglikelihood   <- list()
  out_ids             <- list()
  out_locus_summary   <- list()
  out_locus_lrs       <- list()
  out_locus_dists     <- list()
  out_locinames       <- list()
  out_stan_fits       <- list()


  # Site Loop
  for (site in site_names) {
    if (verbose) message("\n--- Processing Site: ", site, " ---")

    late_site <- late_failures[late_failures$Site == site, ]
    add_site  <- additional[additional$Site == site, ]

    # Remove Site column for internal functions
    late_site <- late_site[, colnames(late_site) != "Site", drop = FALSE]
    add_site  <- add_site[,  colnames(add_site)  != "Site", drop = FALSE]

    # Extract sample IDs
    ids <- unique(gsub(" Day 0", "",
                       late_site$Sample.ID[grepl("Day 0",
                                                 late_site$Sample.ID)]))
    if (length(ids) == 0) {
      warning("WARNING: Site '", site, "' has no Day 0 samples. Skipping this site.", call. = FALSE)
      next
    }

    # A. Define Alleles (Site-specific). define_alleles() already dispatches
    # per locus on marker_info$binning_method, skipping numeric binning for
    # "exact" (ampseq) loci — see R/allele_utils.R.
    allele_definitions <- suppressMessages(define_alleles(rbind(late_site, add_site), marker_info))
    locinames <- names(allele_definitions)
    nloci     <- length(locinames)
   
    if (nloci == 0) {
      warning("WARNING: Site '", site, "' has no valid loci/alleles. Skipping this site.", call. = FALSE)
      next
    }

    # Calculate maximum MOI per site. Matches both length-polymorphic
    # ("locus_1") and ampseq ("locus_allele_1") column-naming conventions.
    marker_cols <- grep("(_allele_\\d+|_\\d+)$", colnames(late_site), value = TRUE)
    maxMOI      <- if (length(marker_cols) > 0) max(as.integer(gsub(".*(_allele_|_)(\\d+)$", "\\2", marker_cols)), na.rm = TRUE) else 1L

    # B. Locus Comparability Matrix
    comparability       <- compute_locus_comparability(late_site, ids, locinames)
    locus_summary       <- comparability$locus_summary
    is_locus_comparable <- comparability$is_locus_comparable

    # C. Prepare Stan Data (Now includes additional_counts)
    if (verbose) message("  Preparing Stan data...")
    stan_data <- prepare_stan_data(
      late_failures_site  = late_site,
      additional_site     = add_site,
      allele_definitions  = allele_definitions,
      marker_info         = marker_info,
      ids                 = ids,
      locinames           = locinames,
      maxMOI              = maxMOI,
      is_locus_comparable = is_locus_comparable
    )

    if (!validate_stan_data(stan_data)) {
      warning("WARNING: Site '", site, "' failed Stan data validation. Skipping this site.", call. = FALSE)
      next
    }

    # Define starting values for the stan sampler
    init_fun <- local({
      stan_data_local <- stan_data
      function() {
        list(
          qq             = 0.1,
          qq_crossfamily = 0.001,
          d_param        = 0.5,
          freq           = lapply(seq_len(stan_data_local$J), function(j) {
            x <- rep(0.1 / stan_data_local$max_K, stan_data_local$max_K)
            x[1:stan_data_local$K[j]] <- 1.0 / stan_data_local$K[j]
            x / sum(x)
          })
        )
      }
    })

    # D. Run Stan Sampling
    if (verbose) message("  Running Stan (",
                         n_chains, " chains, ",
                         iter_sampling, " samples)...")

    fit <- tryCatch({
      stan_model_obj$sample(
        data            = stan_data_only(stan_data),
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

    # E. Extract and Store Results
    extracted <- extract_stan_results(fit, ids, locinames, nloci, length(ids))

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

extract_stan_results <- function(fit, ids, locinames, nloci, nids) {

  # A. Pull the continuous probability of recrudescence
  # draws_matrix: rows = draws (all chains merged), cols = p_recrud[1..N]
  p_recrud_draws <- fit$draws("p_recrud", format = "draws_matrix")

  # B. Pull the log-likelihood (lp__) per chain for diagnostics
  # draws_array: [iterations, chains, variables] — preserves chain structure
  lp_array  <- fit$draws("lp__", format = "draws_array")
  n_iter    <- dim(lp_array)[1]
  n_chains  <- dim(lp_array)[2]
  lp_matrix <- matrix(as.numeric(lp_array), nrow = n_iter, ncol = n_chains)
  loglik_chains <- lapply(seq_len(n_chains), function(ch) lp_matrix[, ch])

  n_draws     <- nrow(p_recrud_draws)
  locus_lrs   <- array(NA_real_, dim = c(nids, nloci, n_draws))
  locus_dists <- array(NA_real_, dim = c(nids, nloci, n_draws))

  return(list(
    p_recrud_draws = p_recrud_draws,
    loglik_chains  = loglik_chains,
    locus_lrs      = locus_lrs,
    locus_dists    = locus_dists
  ))
}
