#' Run Stan MCMC for length-polymorphic marker sites
#'
#' Loops over sites in the data, prepares Stan input, runs HMC sampling, and
#' collects results. Only supports \code{data_type = "length_polymorphic"}.
#'
#' @param late_failures  Data frame of late-failure patient rows (all sites).
#' @param additional     Data frame of additional background allele data (all sites).
#' @param marker_info    Data frame of marker metadata.
#' @param mcmc_config    Named list with fields \code{n_chains}, \code{iter},
#'   \code{burn_in_frac}, \code{random_seed}, and optionally \code{adapt_delta}.
#' @param data_type      Must be \code{"length_polymorphic"}.
#' @param verbose        Logical; print progress messages.
#'
#' @return A named list with elements \code{classifications},
#'   \code{all_chains_loglikelihood}, \code{ids}, \code{locus_summary},
#'   \code{locus_lrs}, \code{locus_dists}, \code{locinames}, and
#'   \code{stan_fits}, each a named list over sites.
#' @noRd
run_stan_sites <- function(late_failures,
                           additional,
                           marker_info,
                           mcmc_config,
                           data_type,
                           verbose = TRUE) {
  
  if (data_type != "length_polymorphic") {
    stop("run_stan_sites() only supports data_type='length_polymorphic'.")
  }
  
  # Extract mcmc configuration parameters 
  n_chains     <- as.integer(mcmc_config$n_chains)
  iter_total   <- as.integer(mcmc_config$iter) 
  burn_in_frac <- as.numeric(mcmc_config$burn_in_frac)
  base_seed    <- as.integer(mcmc_config$random_seed)
  adapt_delta  <- as.numeric(if(!is.null(mcmc_config$adapt_delta)) mcmc_config$adapt_delta else 0.85)
  
  iter_warmup   <- max(200L, as.integer(floor(burn_in_frac * iter_total)))
  iter_sampling <- max(200L, iter_total - iter_warmup)
  
  # Locate and compile model
  stan_model_obj <- tryCatch(stanmodels$malrebay_model, error = function(e) NULL)

  if (is.null(stan_model_obj)) {
    stan_file <- system.file("stan", "malrebay_model.stan", package = "MalReBay")
    if (file.exists(stan_file)) {
      if (verbose) message("Compiling Stan model from source (first time only)...")
      stan_model_obj <- rstan::stan_model(file = stan_file)
    } else {
      stop("Could not find pre-compiled 'malrebay_model' or .stan source file.")
    }
  }
  
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
  
  # Initialization function to prevent Stan from getting stuck at start.
  # Note: theta_recrud is not a model parameter; recrudescence prior is fixed at 0.5.
  init_fun <- function() {
    list(
      qq             = 0.1,
      qq_crossfamily = 0.001,
      d_param        = 0.5
    )
  }
  
  # Site analysis loop
  for (site in site_names) {
    if (verbose) message("\n--- Processing Site: ", site, " ---")
    
    late_site <- late_failures[late_failures$Site == site, ]
    add_site  <- additional[additional$Site == site, ]
    
    # Remove Site column for internal functions
    late_site <- late_site[, colnames(late_site) != "Site", drop = FALSE]
    add_site  <- add_site[,  colnames(add_site)  != "Site", drop = FALSE]
    
    ids  <- unique(gsub(" Day 0", "", late_site$Sample.ID[grepl("Day 0", late_site$Sample.ID)]))
    if (length(ids) == 0) next
    
    # A. Define Alleles (Site-specific)
    allele_definitions <- suppressMessages(define_alleles(rbind(late_site, add_site), marker_info))
    locinames <- names(allele_definitions)
    nloci     <- length(locinames)
    if (nloci == 0) next
    
    marker_cols <- grep("_\\d+$", colnames(late_site), value = TRUE)
    maxMOI      <- if (length(marker_cols) > 0) max(as.integer(gsub(".*_", "", marker_cols)), na.rm = TRUE) else 1L
    
    # B. Locus Comparability Matrix
    locus_summary <- data.frame(patient_id=ids, n_available_d0=0L, n_available_df=0L, n_comparable_loci=0L)
    is_locus_comparable <- matrix(FALSE, nrow=length(ids), ncol=nloci, dimnames=list(ids, locinames))
    
    for (i in seq_along(ids)) {
      pid <- ids[i]
      d0_row <- late_site[grepl(paste0("\\b", pid, " Day 0\\b"), late_site$Sample.ID), ]
      df_row <- late_site[grepl(paste0("\\b", pid, " recurrence\\b"), late_site$Sample.ID), ]
      if (nrow(d0_row) == 0 || nrow(df_row) == 0) next
      for (ln in locinames) {
        lc <- grep(paste0("^", ln, "_"), colnames(late_site), value = TRUE)
        if (any(!is.na(d0_row[, lc]))) locus_summary$n_available_d0[i] <- locus_summary$n_available_d0[i] + 1L
        if (any(!is.na(df_row[, lc]))) locus_summary$n_available_df[i] <- locus_summary$n_available_df[i] + 1L
        if (any(!is.na(d0_row[, lc])) && any(!is.na(df_row[, lc]))) {
          locus_summary$n_comparable_loci[i] <- locus_summary$n_comparable_loci[i] + 1L
          is_locus_comparable[pid, ln] <- TRUE
        }
      }
    }
    
    # C. Prepare Stan Data (Now includes additional_counts)
    if (verbose) message("  Preparing Stan data...")
    stan_input <- prepare_stan_data(
      late_failures_site  = late_site,
      additional_site     = add_site,
      allele_definitions  = allele_definitions,
      marker_info         = marker_info,
      ids                 = ids,
      locinames           = locinames,
      maxMOI              = maxMOI,
      is_locus_comparable = is_locus_comparable
    )

    if (!validate_stan_data(stan_input)) next

    # D. Run Stan Sampling (CALLED ONCE)
    if (verbose) message("  Running Stan (", n_chains, " chains, ", iter_sampling, " samples)...")

    fit <- tryCatch({
      rstan::sampling(
        object  = stan_model_obj,
        data    = stan_data_only(stan_input),
        chains  = n_chains,
        cores   = min(n_chains, parallel::detectCores(logical = FALSE)),
        iter    = iter_warmup + iter_sampling,
        warmup  = iter_warmup,
        init    = init_fun, 
        control = list(adapt_delta = adapt_delta, max_treedepth = 12),
        seed    = base_seed,
        refresh = if (verbose) 100L else 0L
      )
    }, error = function(e) {
      warning("Stan sampling failed for site '", site, "': ", e$message)
      NULL
    })
    
    if (is.null(fit)) next
    
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

#' Extract posterior draws from a fitted Stan model (length-polymorphic)
#'
#' Pulls \code{p_recrud} draws and per-chain log-posterior traces from a
#' \code{stanfit} object returned by \code{run_stan_sites}.
#'
#' @param fit       A \code{stanfit} object.
#' @param ids       Character vector of patient IDs.
#' @param locinames Character vector of locus names.
#' @param nloci     Integer number of loci.
#' @param nids      Integer number of patients.
#'
#' @return A list with \code{p_recrud_draws} (draws x patients matrix),
#'   \code{loglik_chains} (list of per-chain lp__ vectors), \code{locus_lrs},
#'   and \code{locus_dists} (both \code{NA}-filled arrays for now).
#' @noRd
extract_stan_results <- function(fit, ids, locinames, nloci, nids) {

  # 1. Pull the continuous probability of recrudescence
  p_recrud_draws <- rstan::extract(fit, pars = "p_recrud")$p_recrud

  # 2. Pull the log-likelihood (lp__) for diagnostics
  # permuted=FALSE gives array [iterations, chains, 1]; use dim() to get
  # n_chains before subsetting so a single chain doesn't drop to a vector.
  lp_array <- tryCatch(
    rstan::extract(fit, pars = "lp__", permuted = FALSE),
    error = function(e) NULL
  )

  if (!is.null(lp_array) && length(dim(lp_array)) == 3) {
    n_chains      <- dim(lp_array)[2]
    lp_matrix     <- matrix(lp_array[, , 1], ncol = n_chains)
    loglik_chains <- lapply(seq_len(n_chains), function(ch) lp_matrix[, ch])
  } else {
    loglik_chains <- list()
  }
  
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