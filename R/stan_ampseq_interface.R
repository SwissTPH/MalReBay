# =============================================================================
# R/stan_ampseq_interface.R
# Stan interface for amplicon sequencing data (ampseq)
# Parallel structure to stan_interface.R for length-polymorphic data
# =============================================================================

#' Run Stan MCMC for amplicon-sequencing marker sites
#'
#' Loops over sites in the data, prepares Stan input for the ampseq model,
#' runs HMC sampling, and collects results.
#'
#' @param late_failures  Data frame of late-failure patient rows (all sites).
#' @param additional     Data frame of additional background allele data (all sites).
#' @param marker_info    Data frame of marker metadata.
#' @param mcmc_config    Named list with fields \code{n_chains}, \code{iter},
#'   \code{burn_in_frac}, \code{random_seed}, and optionally \code{adapt_delta}.
#' @param verbose        Logical; print progress messages.
#'
#' @return A named list with elements \code{classifications},
#'   \code{all_chains_loglikelihood}, \code{ids}, \code{locus_summary},
#'   \code{locus_lrs}, \code{locus_dists}, \code{locinames}, and
#'   \code{stan_fits}, each a named list over sites.
#' @noRd
run_stan_sites_ampseq <- function(late_failures,
                                   additional,
                                   marker_info,
                                   mcmc_config,
                                   verbose = TRUE) {

  # Extract mcmc configuration parameters 
  n_chains     <- as.integer(mcmc_config$n_chains)
  iter_total   <- as.integer(mcmc_config$iter)
  burn_in_frac <- as.numeric(mcmc_config$burn_in_frac)
  base_seed    <- as.integer(mcmc_config$random_seed)
  adapt_delta  <- as.numeric(if (!is.null(mcmc_config$adapt_delta)) mcmc_config$adapt_delta else 0.85)

  iter_warmup   <- max(200L, as.integer(floor(burn_in_frac * iter_total)))
  iter_sampling <- max(200L, iter_total - iter_warmup)

  # Locate and compile model
  stan_model_obj <- tryCatch(stanmodels$malrebay_ampseq_model, error = function(e) NULL)

  if (is.null(stan_model_obj)) {
    stan_file <- system.file("stan", "malrebay_ampseq_model.stan", package = "MalReBay")
    if (file.exists(stan_file)) {
      if (verbose) message("Compiling Stan ampseq model from source (first time only)...")
      stan_model_obj <- rstan::stan_model(file = stan_file)
    } else {
      stop("Could not find pre-compiled 'malrebay_ampseq_model' or .stan source file.")
    }
  }

  # Output containers
  site_names          <- as.character(unique(late_failures$Site))
  out_classifications <- list()
  out_loglikelihood   <- list()
  out_ids             <- list()
  out_locus_summary   <- list()
  out_locus_lrs       <- list()
  out_locus_dists     <- list()
  out_locinames       <- list()
  out_stan_fits       <- list()

  # Site analysis loop
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
    stan_input <- prepare_stan_data_ampseq(
      late_failures_site  = late_site,
      additional_site     = add_site,
      marker_info         = marker_info,
      ids                 = ids,
      locinames           = locinames,
      maxMOI              = maxMOI,
      is_locus_comparable = is_locus_comparable
    )

    if (!validate_stan_data_ampseq(stan_input)) next

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

#' Extract posterior draws from a fitted Stan ampseq model
#'
#' Pulls \code{p_recrud} draws and per-chain log-posterior traces from a
#' \code{stanfit} object returned by \code{run_stan_sites_ampseq}.
#'
#' @param fit       A \code{stanfit} object.
#' @param ids       Character vector of patient IDs.
#' @param locinames Character vector of locus names.
#' @param nloci     Integer number of loci.
#' @param nids      Integer number of patients.
#'
#' @return A list with \code{p_recrud_draws}, \code{loglik_chains},
#'   \code{locus_lrs}, and \code{locus_dists}.
#' @noRd
extract_stan_results_ampseq <- function(fit, ids, locinames, nloci, nids) {

  p_recrud_draws <- rstan::extract(fit, pars = "p_recrud")$p_recrud

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
