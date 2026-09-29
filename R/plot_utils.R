#' Compute the windowed Gelman-Rubin shrink factor (as \code{coda::gelman.plot} would)
#'
#' @description Reimplements the windowing loop behind \code{coda::gelman.plot()}
#'   using only coda's exported API (its own version, \code{gelman.preplot()},
#'   is internal/unexported), so \code{plot_likelihood_diagnostics()} can build
#'   a ggplot2 version of the same diagnostic from the same numbers coda would
#'   plot: the running median and upper-CI shrink factor from
#'   \code{gelman.diag()}, recomputed over an expanding window of iterations.
#'
#' @param mcmc_list A \code{coda::mcmc.list} with more than one chain.
#' @param max_bins Maximum number of windows, matching \code{gelman.plot()}'s
#'   own \code{max.bins} default.
#' @param confidence Confidence level for the upper bound, matching
#'   \code{gelman.diag()}'s own default.
#' @return A data frame with one row per window: \code{last_iter}, \code{median},
#'   \code{upper}; or \code{NULL} if there are too few iterations to bin.
#' @noRd
compute_gelman_rubin_trace <- function(mcmc_list, max_bins = 50, confidence = 0.95) {
  n_iter <- coda::niter(mcmc_list)
  th     <- coda::thin(mcmc_list)
  nbin   <- min(floor((n_iter - 50) / th), max_bins)
  if (nbin < 1) return(NULL)

  binw      <- floor((n_iter - 50) / nbin)
  last_iter <- c(seq(from = stats::start(mcmc_list) + 50 * th, by = binw * th, length.out = nbin),
                 stats::end(mcmc_list))

  shrink <- t(vapply(last_iter, function(li) {
    psrf <- coda::gelman.diag(stats::window(mcmc_list, end = li),
                              confidence = confidence, autoburnin = FALSE,
                              multivariate = FALSE)$psrf
    as.numeric(psrf[1, ])
  }, numeric(2)))

  data.frame(last_iter = last_iter, median = shrink[, 1], upper = shrink[, 2])
}

#' Plot MCMC Likelihood Diagnostics
#'
#' @description
#' Calculates standard MCMC convergence diagnostics and builds four ggplot2
#' diagnostic panels: Gelman-Rubin shrink factor, Traceplot, Log-Posterior
#' Histogram, and Autocorrelation. The panels are returned as ggplot objects
#' (\code{$traceplot}, \code{$gelman_rubin}, \code{$log_posterior},
#' \code{$autocorrelation}) alongside the numeric diagnostics, so they can be
#' printed, saved, or recombined like any other plot in this package --
#' there's no need to save to PNG and reload just to view them (e.g. in a
#' notebook).
#'
#' @param all_chains_loglikelihood A list where each element is a numeric vector
#'   representing the log-likelihood history of one MCMC chain.
#' @param site_name A character string for labeling plots.
#' @param save_plot A logical. If `TRUE`, panels are also saved as PNG files.
#' @param output_folder A character string specifying the path to save plots.
#' @param verbose A logical. If `TRUE`, prints diagnostic summaries.
#' @param stan_fit An optional \code{CmdStanMCMC} object (from cmdstanr) for additional diagnostics.
#' @param combine_plots A logical. If \code{TRUE}, the returned list also
#'   includes \code{$combined}: all four panels arranged in a single 2x2
#'   \code{ggpubr::ggarrange()} grid.
#'
#' @return An invisible list with the numeric diagnostics (\code{gelman},
#'   \code{ess}, \code{rhat_rank}, \code{ess_bulk}, \code{ess_tail},
#'   \code{geweke}) and the ggplot2 panels (\code{traceplot},
#'   \code{gelman_rubin} -- \code{NULL} with a single chain,
#'   \code{log_posterior}, \code{autocorrelation}, and \code{combined} when
#'   \code{combine_plots = TRUE}).
#'
#' @examples
#' \dontrun{
#'   chains <- list(
#'     c(-10.2, -9.8, -10.5, -9.9, -10.1),
#'     c(-9.9,  -10.3, -10.0, -9.7, -10.4)
#'   )
#'   diag <- plot_likelihood_diagnostics(
#'     all_chains_loglikelihood = chains,
#'     site_name  = "TestSite",
#'     save_plot  = FALSE,
#'     verbose    = FALSE
#'   )
#'   diag$traceplot
#' }
#'
#' @importFrom coda as.mcmc.list mcmc varnames gelman.diag effectiveSize niter thin
#' @importFrom stats acf window start end qnorm
#'
#' @export
plot_likelihood_diagnostics <- function(all_chains_loglikelihood = NULL,
                                         site_name,
                                         stan_fit      = NULL,
                                         save_plot     = TRUE,
                                         output_folder = NULL,
                                         verbose       = TRUE,
                                         combine_plots = FALSE) {

  if (save_plot && is.null(output_folder))
    stop("`output_folder` must be provided when `save_plot = TRUE`.",
         call. = FALSE)

  if (!is.null(stan_fit) && inherits(stan_fit, "CmdStanMCMC") &&
      is.null(all_chains_loglikelihood)) {
    lp_arr <- stan_fit$draws("lp__", format = "draws_array")
    lp_mat <- matrix(as.numeric(lp_arr),
                     nrow = dim(lp_arr)[1],
                     ncol = dim(lp_arr)[2])
    all_chains_loglikelihood <- lapply(
      seq_len(ncol(lp_mat)), function(ch) lp_mat[, ch]
    )
  }

  if (is.null(all_chains_loglikelihood))
    stop("Provide either all_chains_loglikelihood or stan_fit.", call. = FALSE)
  safe_site_name <- gsub(" ", "_", site_name)
  site_dir <- NULL
  if (save_plot) {
    site_dir <- file.path(output_folder, "convergence_diagnosis", safe_site_name)
    if (dir.exists(site_dir)) unlink(site_dir, recursive = TRUE)
    dir.create(site_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Data cleaning
  clean_chains <- lapply(all_chains_loglikelihood, function(x) x[is.finite(x)])
  clean_chains <- Filter(function(x) length(x) >= 2, clean_chains)
  if (length(clean_chains) == 0) {
    warning("No valid chains for site '", site_name, "'. Skipping.", call. = FALSE)
    return(invisible(NULL))
  }
  n_chains <- length(clean_chains)

  loglikelihood_mcmc <- tryCatch({
    mlist  <- lapply(clean_chains, coda::mcmc)
    mclist <- coda::as.mcmc.list(mlist)
    coda::varnames(mclist) <- " "
    mclist
  }, error = function(e) NULL)

  if (is.null(loglikelihood_mcmc)) return(invisible(NULL))

  # Convergence diagnostics
  gelman_result <- NULL
  if (n_chains > 1) {
    gelman_result <- tryCatch(
      coda::gelman.diag(loglikelihood_mcmc, autoburnin = FALSE),
      error = function(e) NULL
    )
  }
  ess_result <- tryCatch(
    coda::effectiveSize(loglikelihood_mcmc),
    error = function(e) NULL
  )
  draws_matrix <- tryCatch(
    do.call(cbind, lapply(clean_chains, as.numeric)),
    error = function(e) NULL
  )
  rhat_rank <- tryCatch(posterior::rhat(draws_matrix),      error = function(e) NA_real_)
  ess_bulk  <- tryCatch(posterior::ess_bulk(draws_matrix),  error = function(e) NA_real_)
  ess_tail  <- tryCatch(posterior::ess_tail(draws_matrix),  error = function(e) NA_real_)
  geweke_result <- tryCatch(
    sapply(loglikelihood_mcmc, function(ch) coda::geweke.diag(ch)$z),
    error = function(e) NULL
  )

  if (verbose) {
    cat("Convergence Diagnostics for:", site_name, "\n")
    cat(strrep("-", 50), "\n")
    if (!is.null(gelman_result)) {
      cat("Classical Gelman-Rubin R-hat:\n"); print(gelman_result)
    } else {
      cat("Classical Gelman-Rubin: not available (near-constant or single chain)\n")
    }
    if (!is.na(rhat_rank))
      cat("\nRank-normalised R-hat (Vehtari 2021):", round(rhat_rank, 4),
          ifelse(rhat_rank < 1.01, "  [PASS]", "  [FAIL]"), "\n")
    if (!is.null(ess_result))
      cat("\nClassical ESS (coda):", round(as.numeric(ess_result), 1), "\n")
    if (!is.na(ess_bulk))
      cat("Bulk ESS:", round(ess_bulk, 1),
          ifelse(ess_bulk > 400, "  [PASS]", "  [FAIL]"), "\n")
    if (!is.na(ess_tail))
      cat("Tail ESS:", round(ess_tail, 1),
          ifelse(ess_tail > 400, "  [PASS]", "  [FAIL]"), "\n")
    if (!is.null(geweke_result)) {
      cat("\nGeweke Z-scores per chain |Z| < 1.96 = stationary:\n")
      for (i in seq_along(geweke_result)) {
        z <- geweke_result[i]
        cat(sprintf("  Chain %d: Z = %6.4f  %s\n", i, z,
                    ifelse(abs(z) < 1.96, "[PASS]", "[FAIL]")))
      }
    }
    cat(strrep("-", 50), "\n")
  }

  # ---- Build the four ggplot2 panels ----

  # Shared title styling: bold centred title, smaller grey subtitle underneath
  # describing what good convergence looks like on that panel.
  title_prefix <- paste("Site", site_name)
  title_theme  <- ggplot2::theme(
    plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, color = "grey30",
                                          lineheight = 1.1)
  )

  # Traceplot: log-posterior per iteration, one line per chain
  trace_df <- do.call(rbind, lapply(seq_len(n_chains), function(i) {
    data.frame(iteration = seq_along(clean_chains[[i]]),
              log_posterior = clean_chains[[i]],
              chain = factor(paste("Chain", i)))
  }))
  p_trace <- ggplot2::ggplot(
    trace_df,
    ggplot2::aes(x = .data$iteration, y = .data$log_posterior, color = .data$chain)
  ) +
    ggplot2::geom_line(alpha = 0.8, linewidth = 0.3) +
    ggplot2::labs(title = paste(title_prefix, "- Traceplot"),
                 subtitle = paste0("Good convergence: chains overlap and fluctuate around a stable level,\n",
                                   "like a 'caterpillar', with no trends or stuck chains"),
                 x = "Iterations", y = "Log-Posterior", color = NULL) +
    ggplot2::theme_classic(base_size = 14) +
    title_theme

  # Gelman-Rubin shrink factor: same windowed calculation coda::gelman.plot() uses
  p_gelman <- NULL
  if (n_chains > 1) {
    gr_trace <- tryCatch(compute_gelman_rubin_trace(loglikelihood_mcmc),
                         error = function(e) NULL)
    if (!is.null(gr_trace)) {
      gr_long <- data.frame(
        last_iter = rep(gr_trace$last_iter, 2),
        shrink    = c(gr_trace$median, gr_trace$upper),
        stat      = factor(rep(c("median", "97.5%"), each = nrow(gr_trace)),
                           levels = c("median", "97.5%"))
      )
      p_gelman <- ggplot2::ggplot(
        gr_long, ggplot2::aes(x = .data$last_iter, y = .data$shrink, color = .data$stat)
      ) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept = 1, color = "grey40") +
        ggplot2::geom_hline(yintercept = 1.1, color = "grey40", linetype = "dotted") +
        ggplot2::scale_color_manual(values = c(median = "black", `97.5%` = "indianred")) +
        ggplot2::labs(title = paste(title_prefix, "- Gelman-Rubin Diagnostic"),
                     subtitle = paste0("Good convergence: both lines drop towards 1 and stay close to it\n",
                                       "(shrink factor < 1.1) by the end of the chains"),
                     x = "last iteration in chain", y = "shrink factor", color = NULL) +
        ggplot2::theme_classic(base_size = 14) +
        title_theme
    }
  }

  # Log-posterior distribution, chains pooled
  p_hist <- ggplot2::ggplot(
    data.frame(log_posterior = unlist(clean_chains)),
    ggplot2::aes(x = .data$log_posterior)
  ) +
    ggplot2::geom_histogram(bins = 40, fill = "steelblue", color = "white") +
    ggplot2::labs(title = paste(title_prefix, "- Log-Posterior Distribution"),
                 subtitle = paste0("Good convergence: a smooth, single-peaked distribution\n",
                                   "without multiple modes or long isolated tails"),
                 x = "Log-Posterior", y = "Frequency") +
    ggplot2::theme_classic(base_size = 14) +
    title_theme

  # Autocorrelation, one facet per chain -- same 95% CI band formula stats::acf() uses
  acf_df <- do.call(rbind, lapply(seq_len(n_chains), function(i) {
    a <- stats::acf(clean_chains[[i]], lag.max = 50, plot = FALSE)
    data.frame(lag = as.numeric(a$lag), acf = as.numeric(a$acf),
              chain = paste("Chain", i))
  }))
  ci_df <- data.frame(
    chain = paste("Chain", seq_len(n_chains)),
    ci    = vapply(clean_chains, function(x) stats::qnorm(0.975) / sqrt(length(x)), numeric(1))
  )
  p_acf <- ggplot2::ggplot(acf_df, ggplot2::aes(x = .data$lag, y = .data$acf)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_segment(ggplot2::aes(xend = .data$lag, yend = 0)) +
    ggplot2::geom_hline(data = ci_df, ggplot2::aes(yintercept = .data$ci),
                        color = "blue", linetype = "dashed") +
    ggplot2::geom_hline(data = ci_df, ggplot2::aes(yintercept = -.data$ci),
                        color = "blue", linetype = "dashed") +
    ggplot2::facet_wrap(~ chain) +
    ggplot2::labs(title = paste(title_prefix, "- Autocorrelation"),
                 subtitle = paste0("Good convergence: autocorrelation drops quickly towards 0\n",
                                   "and stays mostly within the dashed 95% bounds"),
                 x = "Lag", y = "ACF") +
    ggplot2::theme_classic(base_size = 14) +
    title_theme +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "grey90", color = "black"))

  p_combined <- NULL
  if (combine_plots) {
    panels <- Filter(Negate(is.null), list(p_trace, p_gelman, p_hist, p_acf))
    p_combined <- ggpubr::ggarrange(plotlist = panels, ncol = 2, nrow = 2)
  }

  if (save_plot) {
    ggplot2::ggsave(file.path(site_dir, paste0(safe_site_name, "_traceplot.png")),
                    p_trace, width = 1000 / 120, height = 600 / 120, dpi = 120)
    if (!is.null(p_gelman))
      ggplot2::ggsave(file.path(site_dir, paste0(safe_site_name, "_gelman_rubin.png")),
                      p_gelman, width = 1000 / 120, height = 700 / 120, dpi = 120)
    ggplot2::ggsave(file.path(site_dir, paste0(safe_site_name, "_log_likelihood_distribution.png")),
                    p_hist, width = 1000 / 120, height = 600 / 120, dpi = 120)
    ggplot2::ggsave(file.path(site_dir, paste0(safe_site_name, "_autocorrelation.png")),
                    p_acf, width = 1200 / 120, height = 1000 / 120, dpi = 120)
  }

  invisible(list(
    gelman          = gelman_result,
    ess             = ess_result,
    rhat_rank       = rhat_rank,
    ess_bulk        = ess_bulk,
    ess_tail        = ess_tail,
    geweke          = geweke_result,
    traceplot       = p_trace,
    gelman_rubin    = p_gelman,
    log_posterior   = p_hist,
    autocorrelation = p_acf,
    combined        = p_combined
  ))
}

#' Plot Posterior Probability Histogram
#'
#' Creates and saves histograms of posterior probabilities of recrudescence:
#' one pooling all patients across sites, and one with a facet per site. Bars
#' are coloured by classification (dark red for recrudescence, light blue for
#' new infection; ColorBrewer RdBu, colourblind-safe) and a dotted vertical line
#' marks the classification threshold. This is an internal function called
#' automatically by \code{\link{MalReBay}}.
#'
#' @param summary_results A list returned by \code{\link{summarise_results}},
#'   containing a \code{posterior_probabilities} data frame with a
#'   \code{Probability} column, a \code{Site} column (needed for the per-site
#'   plot) and, optionally, the \code{prob_threshold} used for classification.
#' @param output_folder A string specifying the directory where the histogram
#'   PNGs will be saved. If \code{NULL}, the plots are displayed instead.
#' @param verbose Logical. If \code{TRUE}, prints a message when the files are
#'   saved. Defaults to \code{TRUE}.
#' @param prob_threshold Numeric in \verb{[0, 1]}. The posterior probability
#'   at or above which a recurrence counts as recrudescence. Defaults to
#'   \code{summary_results$prob_threshold}, or \code{0.5} if that is missing.
#'
#' @return An invisible list of ggplot objects: \code{all_sites} (all patients
#'   pooled) and \code{by_site} (one facet per site, with a free y-axis so
#'   small sites stay readable; \code{NULL} if there is no \code{Site}
#'   column). When \code{output_folder} is \code{NULL} the plots are also
#'   printed; otherwise they are saved as
#'   \code{recrudescence_probability_histogram.png} and
#'   \code{recrudescence_probability_histogram_by_site.png} in
#'   \code{output_folder}. Returns \code{invisible(NULL)} if no probabilities
#'   are available.
#'
#' @examples
#' summary_results <- list(
#'   posterior_probabilities = data.frame(
#'     Patient_ID  = c("P1", "P2", "P3", "P4"),
#'     Site        = c("A", "A", "B", "B"),
#'     Probability = c(0.12, 0.87, 0.45, 0.93)
#'   ),
#'   prob_threshold = 0.5
#' )
#' plots <- plot_probability_histogram(summary_results)
#' plots$by_site
#'
#' @export
plot_probability_histogram <- function(summary_results, output_folder = NULL, verbose = TRUE,
                                       prob_threshold = summary_results$prob_threshold) {

  posterior_probabilities <- summary_results$posterior_probabilities

  if (is.null(posterior_probabilities) || nrow(posterior_probabilities) == 0) {
    warning("No posterior probabilities to plot.")
    return(invisible(NULL))
  }

  if (is.null(prob_threshold)) prob_threshold <- 0.5

  probs <- as.numeric(as.character(posterior_probabilities$Probability))
  hist_df <- data.frame(
    Probability    = probs,
    Classification = factor(ifelse(probs >= prob_threshold, "Recrudescence", "New infection"),
                            levels = c("Recrudescence", "New infection"))
  )
  has_site <- "Site" %in% names(posterior_probabilities)
  if (has_site) hist_df$Site <- posterior_probabilities$Site

  # Stacked by class, so a bin straddling a threshold that isn't on a bin
  # edge is split into its two colours rather than mis-coloured. Bins are
  # left-closed to match the ">= threshold" rule: a probability exactly on the
  # threshold lands in the bin to the right of the line.
  p_all <- ggplot2::ggplot(hist_df, ggplot2::aes(x = .data$Probability, fill = .data$Classification)) +
    ggplot2::geom_histogram(breaks = seq(0, 1, by = 0.05), closed = "left",
                            color = "white", position = "stack") +
    ggplot2::geom_vline(xintercept = prob_threshold, linetype = "dotted",
                        color = "black", linewidth = 0.8) +
    ggplot2::scale_fill_manual(values = c(Recrudescence = "#B2182B", `New infection` = "#92C5DE"),
                               drop = FALSE) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
    # Patient counts: whole-number ticks only (small per-site panels otherwise get 2.5, 7.5, ...)
    ggplot2::scale_y_continuous(breaks = function(lims) unique(floor(pretty(lims)))) +
    ggplot2::labs(title    = "Posterior Probability Distribution - All Sites",
                  subtitle = paste("Dotted line: classification threshold =", prob_threshold),
                  x = "Probability of Recrudescence", y = "Number of Patients", fill = NULL) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, color = "grey30")
    )

  p_site <- NULL
  if (has_site) {
    n_sites   <- length(unique(hist_df$Site))
    facet_col <- ceiling(sqrt(n_sites))
    facet_row <- ceiling(n_sites / facet_col)
    p_site <- suppressMessages(p_all +
      ggplot2::facet_wrap(~ Site, ncol = facet_col, scales = "free_y") +
      ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25))) +
      ggplot2::labs(title = "Posterior Probability Distribution by Site") +
      ggplot2::theme(strip.background = ggplot2::element_rect(fill = "grey90", color = "black"))
  } else if (verbose) {
    message("INFO: No 'Site' column in posterior_probabilities; skipping per-site histogram.")
  }

  if (!is.null(output_folder)) {
    if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
    ggplot2::ggsave(file.path(output_folder, "recrudescence_probability_histogram.png"),
                    p_all, width = 8, height = 6, units = "in", dpi = 300)
    if (!is.null(p_site))
      ggplot2::ggsave(file.path(output_folder, "recrudescence_probability_histogram_by_site.png"),
                      p_site, width = 3.5 * facet_col + 2, height = 3 * facet_row + 1.5,
                      units = "in", dpi = 300)
    if (verbose) message("INFO: Probability histograms saved to: ", output_folder)
  } else {
    print(p_all)
    if (!is.null(p_site)) print(p_site)
  }

  invisible(list(all_sites = p_all, by_site = p_site))
}

#' Plot Multiplicity of Infection (MOI)
#'
#' This function calculates the MOI (number of distinct alleles per marker per
#' sample) and generates one violin plot per site for improved readability,
#' with a narrow boxplot (median/IQR) overlaid inside each violin and jittered
#' points for the individual samples.
#'
#' @param genotypedata A data frame containing `Sample.ID`, `Site`, and marker columns.
#' @param marker_pattern A regex pattern to identify marker columns.
#'   Defaults to standard allele suffix patterns.
#' @param output_folder Path to the directory where plots will be saved.
#'   If NULL, plots are not saved to disk.
#' @param filename_prefix Prefix for output PNG filenames. Each file will be
#'   named `<prefix>_<site>.png`.
#'
#' @return A named list of ggplot objects, one per site. The underlying MOI
#'   data is attached to each plot as the attribute `"moi_data"`.
#'
#' @examples
#' \dontrun{
#'   gdata <- data.frame(
#'     Sample.ID = c("P1 Day 0", "P1 recurrence", "P2 Day 0", "P2 recurrence"),
#'     Site      = "SiteA",
#'     TA1_1     = c(174, 177, 162, 162),
#'     TA1_2     = c(177,  NA, 171,  NA)
#'   )
#'   plot_moi(genotypedata = gdata, output_folder = NULL)
#' }
#' @export
plot_moi <- function(genotypedata,
                     marker_pattern  = "(_allele_\\d+|_\\d+)$",
                     output_folder   = NULL,
                     filename_prefix = "moi_per_marker") {

  # Input validation
  if (!"Site" %in% colnames(genotypedata)) {
    stop("Input 'genotypedata' must contain a column named 'Site'.")
  }

  genotypedata$Sample.ID <- as.character(genotypedata$Sample.ID)

  marker_cols <- grep(marker_pattern, colnames(genotypedata), value = TRUE)
  if (length(marker_cols) == 0) {
    warning("No marker columns found matching the pattern.")
    return(NULL)
  }

  # Compute MOI 
  moi_data <- genotypedata %>%
    tidyr::pivot_longer(
      cols          = dplyr::all_of(marker_cols),
      names_to      = "marker_replicate",
      values_to     = "allele",
      values_drop_na = TRUE
    ) %>%
    dplyr::mutate(
      marker_id = gsub(marker_pattern, "", .data$marker_replicate)
    ) %>%
    dplyr::group_by(.data$Sample.ID, .data$Site, .data$marker_id) %>%
    dplyr::summarise(MOI = dplyr::n_distinct(.data$allele), .groups = "drop")

  # Fill zeroes for samples / markers with no data
  all_markers   <- unique(gsub(marker_pattern, "", marker_cols))
  complete_grid <- tidyr::crossing(
    dplyr::distinct(genotypedata, .data$Sample.ID, .data$Site),
    marker_id = all_markers
  )

  moi_data <- dplyr::left_join(
    complete_grid, moi_data,
    by = c("Sample.ID", "Site", "marker_id")
  ) %>%
    dplyr::mutate(MOI = tidyr::replace_na(.data$MOI, 0))

  if (nrow(moi_data) == 0) {
    message("MOI data is empty. Skipping plot generation.")
    return(NULL)
  }

  # Define label data
  label_data <- moi_data %>%
    dplyr::group_by(.data$Site, .data$marker_id) %>%
    dplyr::summarise(mean_moi = mean(.data$MOI), .groups = "drop") %>%
    dplyr::left_join(
      moi_data %>%
        dplyr::group_by(.data$Site) %>%
        dplyr::summarise(label_y_pos = max(.data$MOI) + 0.5, .groups = "drop"),
      by = "Site"
    )

  marker_levels <- unique(label_data$marker_id)

  # Boxplot summary computed by hand (stat = "identity") instead of letting
  # geom_boxplot's own stat derive whiskers from 1.5*IQR: MOI is a small
  # integer count that's frequently constant within a marker (Q1 = Q3 =
  # median), and a zero IQR collapses the default whisker rule to nothing --
  # every non-median point then counts as an "outlier" and (with
  # outlier.shape = NA) simply vanishes, leaving a boxplot with no whiskers.
  # Drawing ymin/ymax as the true min/max keeps the box "complete" regardless.
  box_data <- moi_data %>%
    dplyr::group_by(.data$Site, .data$marker_id) %>%
    dplyr::summarise(
      ymin   = min(.data$MOI),
      lower  = stats::quantile(.data$MOI, 0.25, type = 7),
      middle = stats::median(.data$MOI),
      upper  = stats::quantile(.data$MOI, 0.75, type = 7),
      ymax   = max(.data$MOI),
      .groups = "drop"
    )

  # Output folder
  if (!is.null(output_folder) && !dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }

  # Build a plot per site 
  sites <- unique(moi_data$Site)

  plots <- lapply(stats::setNames(sites, sites), function(site) {

    site_moi    <- dplyr::filter(moi_data,   .data$Site == site)
    site_labels <- dplyr::filter(label_data, .data$Site == site)
    site_box    <- dplyr::filter(box_data,   .data$Site == site)

    site_moi$marker_id <- factor(site_moi$marker_id, levels = marker_levels)
    site_box$marker_id <- factor(site_box$marker_id, levels = marker_levels)

    p <- ggplot2::ggplot(
      site_moi,
      ggplot2::aes(
        x     = .data$marker_id,
        y     = .data$MOI,
        fill  = .data$marker_id,
        color = .data$marker_id
      )
    ) +
      ggplot2::geom_jitter(width = 0.15, height = 0.1, alpha = 0.3, size = 0.6) +
      ggplot2::geom_violin(alpha = 0.4, trim = FALSE) +
      ggplot2::geom_boxplot(
        data        = site_box,
        mapping     = ggplot2::aes(
          x      = .data$marker_id,
          ymin   = .data$ymin,
          lower  = .data$lower,
          middle = .data$middle,
          upper  = .data$upper,
          ymax   = .data$ymax
        ),
        stat        = "identity",
        width       = 0.12,
        color       = "grey30",
        fill        = "white",
        alpha       = 0.7,
        inherit.aes = FALSE
      ) +
      ggplot2::geom_text(
        data    = site_labels,
        mapping = ggplot2::aes(
          x     = .data$marker_id,
          y     = .data$label_y_pos,
          label = sprintf("%.2f", .data$mean_moi)
        ),
        size    = 3.5,
        vjust   = 0,
        color   = "black"
      ) +
      ggplot2::scale_fill_brewer(palette  = "Set2") +
      ggplot2::scale_color_brewer(palette = "Set2") +
      ggplot2::scale_y_continuous(breaks = scales::pretty_breaks()) +
      ggplot2::coord_cartesian(ylim = c(-0.5, NA)) +
      ggplot2::labs(
        title = paste("MOI by Marker \u2013", site),
        x     = "Marker",
        y     = "MOI (Number of Alleles)"
      ) +
      ggplot2::theme_classic(base_size = 14) +
      ggplot2::theme(
        legend.position  = "none",
        plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold"),
        axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
        strip.background = ggplot2::element_rect(fill = "grey90", color = "black"),
        strip.text.y     = ggplot2::element_text(angle = 0, face = "bold")
      )

    # Save individual plot
    if (!is.null(output_folder)) {
      safe_site  <- gsub("[^A-Za-z0-9_-]", "_", site)
      output_path <- file.path(
        output_folder,
        paste0(filename_prefix, "_", safe_site, ".png")
      )
      ggplot2::ggsave(output_path, plot = p, width = 12, height = 6)
    }

    attr(p, "moi_data") <- site_moi
    p
  })

  plots
}


#' Define Allele Bins for Plotting
#'
#' @description Internal helper to group raw allele lengths into bins based on
#'   the binning method defined in `marker_info`.
#'
#' @param genotypedata A data frame containing genotyping data.
#' @param marker_info A data frame with marker definitions.
#' @return A list of data frames, each defining allele bins for a marker.
#' @noRd
define_alleles_for_plotting <- function(genotypedata, marker_info) {
  alleles_definitions_bin  <- list()
  data_marker_columns      <- grep("(_allele_\\d+|_\\d+)$", colnames(genotypedata), value = TRUE)
  available_base_markers   <- unique(gsub("(_allele_\\d+|_\\d+)$", "", data_marker_columns))

  for (locus_name in available_base_markers) {
    locus_marker_info <- marker_info[marker_info$marker_id == locus_name, ]
    if (nrow(locus_marker_info) == 0) next

    binning_method    <- locus_marker_info$binning_method[1]
    locus_cols        <- grepl(paste0("^", locus_name, "(_allele_|_)\\d+"), colnames(genotypedata))
    raw_alleles       <- unlist(genotypedata[, locus_cols])
    unique_alleles    <- sort(unique(raw_alleles[!is.na(raw_alleles)]))
    if (length(unique_alleles) == 0) next

    bins <- if (binning_method == "microsatellite") {
      repeat_length <- locus_marker_info$repeatlength[1]
      ceiling((unique_alleles - unique_alleles[1] + 1) / repeat_length)
    } else if (binning_method == "msp_glurp") {
      gap_threshold <- locus_marker_info$repeatlength[1]
      breaks        <- c(0, which(diff(unique_alleles) > gap_threshold), length(unique_alleles))
      findInterval(seq_along(unique_alleles), breaks)
    } else {
      next
    }

    alleles_definitions_bin[[locus_name]] <- data.frame(
      min = tapply(unique_alleles, bins, min),
      max = tapply(unique_alleles, bins, max)
    )
  }
  alleles_definitions_bin
}

#' Process Data for Pie Chart
#'
#' @description Internal helper to prepare data for `plot_pie_chart`. It renames
#'   columns, creates labels, and orders the data.
#'
#' @param plot_data_df A data frame with Haplotype, Amount, and Frequency.
#' @param data_type The type of data (`"length_polymorphic"` or `"ampseq"`).
#' @return A processed data frame ready for plotting.
#' @noRd
process_pie_data <- function(plot_data_df, data_type) {
  colnames(plot_data_df) <- c("Haplotype", "Amount", "Frequency")
  plot_data_df <- plot_data_df[!is.na(plot_data_df$Haplotype), ]

  plot_data_df$Haplotype <- if (data_type == "length_polymorphic") {
    factor(plot_data_df$Haplotype, levels = sort(as.numeric(as.character(plot_data_df$Haplotype))))
  } else {
    as.factor(plot_data_df$Haplotype)
  }

  plot_data_df$Label <- paste0(round(plot_data_df$Frequency * 100), "%")
  plot_data_df       <- plot_data_df[order(plot_data_df$Amount, decreasing = TRUE), ]
  if (nrow(plot_data_df) > 4) plot_data_df[5:nrow(plot_data_df), "Label"] <- ""
  plot_data_df
}

#' Plot a Pie Chart for Haplotype/Allele Diversity
#'
#' @description Internal helper function to create a single pie chart.
#'
#' @param data_df A dataframe with Haplotype/Allele, Amount, and Frequency.
#' @param color_marker The base color for the chart palette.
#' @param marker_name The name of the genetic marker.
#' @param total_n The total number of infections for this marker.
#' @param data_type The type of data (`"length_polymorphic"` or `"ampseq"`).
#' @return A ggplot object representing the pie chart.
#'
#' @importFrom dplyr mutate lead if_else
#' @importFrom ggplot2 ggplot aes geom_col coord_polar scale_fill_manual geom_text theme_void theme element_text ggtitle
#' @importFrom forcats fct_inorder
#' @importFrom grDevices colorRampPalette
#' @noRd
plot_pie_chart <- function(data_df, color_marker, marker_name, total_n, data_type) {
  if (nrow(data_df) == 0) return(NULL)

  plot_data_df  <- process_pie_data(data_df, data_type)
  title_marker  <- paste0(marker_name, "\n(n=", total_n, ")")
  colfunc       <- grDevices::colorRampPalette(c(color_marker, "white"))

  plot_data_df <- plot_data_df %>%
    dplyr::mutate(
      csum = rev(cumsum(rev(.data$Amount))),
      pos  = .data$Amount / 2 + dplyr::lead(.data$csum, 1),
      pos  = dplyr::if_else(is.na(.data$pos), .data$Amount / 2, .data$pos)
    )

  ggplot2::ggplot(
    plot_data_df,
    ggplot2::aes(x = "", y = .data$Amount, fill = forcats::fct_inorder(.data$Haplotype))
  ) +
    ggplot2::geom_col(width = 1, color = "white") +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::scale_fill_manual(values = colfunc(nrow(plot_data_df))) +
    ggplot2::geom_text(ggplot2::aes(y = .data$pos, label = .data$Label), size = 4, color = "black") +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none",
      plot.title      = ggplot2::element_text(hjust = 0.5, size = 12)
    ) +
    ggplot2::ggtitle(title_marker)
}

# Internal helpers 

#' Pivot and tidy genotype data into long form
#' @noRd
.pivot_long_genotypes <- function(data, value_col) {
  data %>%
    tidyr::pivot_longer(
      cols        = dplyr::matches("_allele_\\d+$|_\\d+$"),
      names_to    = "marker_replicate",
      values_to   = value_col
    ) %>%
    dplyr::filter(!is.na(.data[[value_col]]) & .data[[value_col]] != "") %>%
    dplyr::mutate(marker_id = gsub("_allele_\\d+$|_\\d+$", "", .data$marker_replicate)) %>%
    dplyr::select("Sample.ID", "marker_id", dplyr::all_of(value_col)) %>%
    dplyr::distinct()
}

#' Compute per-marker allele frequencies from a long-format table
#' @noRd
.compute_frequencies <- function(long_data, allele_col) {
  long_data %>%
    dplyr::group_by(.data$marker_id, .data[[allele_col]]) %>%
    dplyr::summarise(Amount = dplyr::n(), .groups = "drop") %>%
    dplyr::left_join(
      long_data %>%
        dplyr::group_by(.data$marker_id) %>%
        dplyr::summarise(TotalInfections = dplyr::n(), .groups = "drop"),
      by = "marker_id"
    ) %>%
    dplyr::mutate(Frequency = .data$Amount / .data$TotalInfections)
}

#' Bin length-polymorphic alleles using pre-computed bin definitions
#' @noRd
.bin_lp_alleles <- function(long_data, alleles_definitions_bin) {
  markers_with_bins <- intersect(unique(long_data$marker_id), names(alleles_definitions_bin))

  purrr::map_dfr(markers_with_bins, function(marker) {
    bin_centers <- rowMeans(alleles_definitions_bin[[marker]], na.rm = TRUE)
    long_data %>%
      dplyr::filter(.data$marker_id == marker) %>%
      dplyr::mutate(
        true_alleles = bin_centers[
          sapply(.data$allele_length, function(x) which.min(abs(x - bin_centers)))
        ]
      )
  }) %>%
    dplyr::select("Sample.ID", "marker_id", "true_alleles") %>%
    dplyr::distinct()
}

#' Generate and Save Diversity Pie Charts
#'
#' @description Creates pie charts visualising allele or haplotype diversity
#'   across all samples combined (Day 0 and recurrences pooled). Pies are
#'   arranged with at most `max_cols` per row.
#'
#' @param genotypedata A dataframe containing the genotyping data.
#' @param data_type A string: `"length_polymorphic"` or `"ampseq"`.
#' @param marker_info A dataframe with marker definitions; required when
#'   `data_type = "length_polymorphic"`.
#' @param output_folder Path to the directory where the output PNG will be saved.
#' @param filename_prefix A string prefix for the output filename.
#' @param max_cols Maximum number of pie charts per row. Defaults to 4.
#' @return Invisibly returns the combined ggplot object, or `NULL` if no data.
#'
#' @examples
#' \dontrun{
#'   gdata <- data.frame(
#'     Sample.ID     = c("P1 Day 0", "P1 recurrence", "P2 Day 0", "P2 recurrence"),
#'     Site          = "SiteA",
#'     cpmp_allele_1 = c("HAPL_A", "HAPL_A", "HAPL_A", "HAPL_B"),
#'     cpmp_allele_2 = c(NA, NA, NA, NA)
#'   )
#'   plot_markers_diversity(genotypedata = gdata, data_type = "ampseq")
#' }
#'
#' @importFrom dplyr %>% filter mutate select distinct group_by summarise left_join all_of n
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggsave
#' @importFrom ggpubr ggarrange
#' @importFrom RColorBrewer brewer.pal
#' @importFrom purrr map_dfr
#' @export
plot_markers_diversity <- function(genotypedata,
                                   data_type,
                                   marker_info = NULL,
                                   output_folder = NULL,
                                   filename_prefix = "diversity",
                                   max_cols        = 4) {

  # Validation
  if (!data_type %in% c("length_polymorphic", "ampseq")) {
    stop("'data_type' must be \"length_polymorphic\" or \"ampseq\".")
  }
  if (data_type == "length_polymorphic" && is.null(marker_info)) {
    stop("'marker_info' must be provided when data_type = \"length_polymorphic\".")
  }

  sid_col <- grep("^sample.?id$", colnames(genotypedata), ignore.case = TRUE, value = TRUE)
  if (length(sid_col) == 0) stop("Input data must contain a 'Sample.ID' column.")
  genotypedata$Sample.ID <- as.character(genotypedata[[sid_col[1]]])

  if (!is.null(output_folder)) {
    if (!dir.exists(output_folder)) {
      dir.create(output_folder, recursive = TRUE)
    }
  }
  

  # Build combined frequency data
  if (data_type == "length_polymorphic") {
    alleles_definitions_bin <- define_alleles_for_plotting(genotypedata, marker_info)

    long_data      <- .pivot_long_genotypes(genotypedata, "allele_length")
    binned_data    <- .bin_lp_alleles(long_data, alleles_definitions_bin)
    frequency_data <- .compute_frequencies(binned_data, "true_alleles")
    allele_col     <- "true_alleles"

  } else {
    # Normalise Sample.ID tags for ampseq if needed
    if (!any(grepl(" Day ", genotypedata$Sample.ID))) {
      genotypedata$Sample.ID <- gsub("D0$",          " Day 0",       genotypedata$Sample.ID)
      genotypedata$Sample.ID <- gsub("D[1-9][0-9]*$", " Day Failure", genotypedata$Sample.ID)
    }

    long_data      <- .pivot_long_genotypes(genotypedata, "haplotype")
    frequency_data <- .compute_frequencies(long_data, "haplotype")
    allele_col     <- "haplotype"
  }

  if (nrow(frequency_data) == 0) {
    message("No frequency data could be computed. Skipping plot generation.")
    return(invisible(NULL))
  }

  # Build pie charts 
  list_markers <- unique(frequency_data$marker_id)
  n_markers    <- length(list_markers)

  base_colors <- RColorBrewer::brewer.pal(max(3, min(n_markers, 8)), "Set2")
  colors      <- rep_len(base_colors, n_markers)

  p_array <- Map(function(marker_name, color) {
    plot_data <- dplyr::filter(frequency_data, .data$marker_id == marker_name)
    plot_pie_chart(
      dplyr::select(plot_data, dplyr::all_of(c(allele_col, "Amount", "Frequency"))),
      color, marker_name, plot_data$TotalInfections[1], data_type
    )
  }, list_markers, colors)

  p_array <- Filter(Negate(is.null), p_array)
  if (length(p_array) == 0) {
    message("No pie charts were generated.")
    return(invisible(NULL))
  }

  # Arrange and save 
  num_cols     <- min(length(p_array), max_cols)
  num_rows     <- ceiling(length(p_array) / num_cols)
  final_figure <- ggpubr::ggarrange(plotlist = p_array, ncol = num_cols, nrow = num_rows)

  if (!is.null(output_folder)) {
    output_path <- file.path(output_folder, paste0(filename_prefix, "_", data_type, "_comparison.png"))
    ggplot2::ggsave(output_path, plot = final_figure,
                    width  = num_cols * 4,
                    height = num_rows * 4,
                    limitsize = FALSE)
  }

  invisible(final_figure)
}

#' Plot Allele Distribution
#'
#' @description Creates a distribution plot of raw allele values for each
#'   marker, pooling all samples and timepoints together, with the y-axis
#'   showing frequency (percentage of that marker's calls) rather than raw
#'   counts. For length-polymorphic markers (microsatellites, MSP1/MSP2/GLURP)
#'   this is a histogram of fragment sizes; for AmpSeq markers
#'   (\code{binning_method = "exact"}), which have no numeric length, this is
#'   instead a bar chart of frequency per haplotype, labelled by haplotype
#'   name. One panel per marker, arranged in a grid with at most
#'   \code{max_cols} per row.
#'
#' @param genotypedata A data frame containing the genotyping data, with a
#'   \code{Sample.ID} column and per-marker allele columns.
#' @param marker_info A data frame with marker definitions (\code{marker_id},
#'   \code{repeatlength}, \code{binning_method}). Used to size histogram bins
#'   to each marker's natural repeat unit (unless \code{binwidth} is given)
#'   and to decide which markers get a haplotype bar chart instead.
#' @param binwidth Fixed bin width to use for every length-polymorphic
#'   marker's histogram (ignored for AmpSeq markers). If \code{NULL}
#'   (default), microsatellite markers use their own \code{repeatlength}
#'   from \code{marker_info} (the true repeat-unit size) as long as it
#'   yields at most 60 bins across the observed allele range for that
#'   marker; MSP1/MSP2/GLURP markers, and any marker where
#'   \code{repeatlength} would produce more than 60 bins (e.g. because a
#'   few outlier alleles stretch the range), use a data-driven default
#'   instead, since their \code{repeatlength} either encodes a
#'   family-clustering gap threshold rather than a natural bin width, or
#'   would otherwise chop the range into slivers too thin to see.
#' @param output_folder Path to the directory where the output PNG will be
#'   saved. If \code{NULL} (default), the plot is not saved to disk.
#' @param filename_prefix A string prefix for the output filename.
#' @param max_cols Maximum number of panels per row. Defaults to 4.
#' @return Invisibly returns the combined ggplot object, or \code{NULL} if no
#'   allele data is available to plot.
#'
#' @examples
#' \dontrun{
#'   gdata <- data.frame(
#'     Sample.ID = c("P1 Day 0", "P1 recurrence", "P2 Day 0", "P2 recurrence"),
#'     Site      = "SiteA",
#'     TA1_1     = c(174, 177, 162, 171),
#'     TA1_2     = c(NA, NA, NA, NA)
#'   )
#'   marker_info <- data.frame(marker_id = "TA1", repeatlength = 3,
#'                              binning_method = "microsatellite")
#'   plot_allele_distribution(gdata, marker_info)
#' }
#'
#' @export
plot_allele_distribution <- function(genotypedata,
                                     marker_info,
                                     binwidth        = NULL,
                                     output_folder   = NULL,
                                     filename_prefix = "allele_distribution",
                                     max_cols        = 4) {

  sid_col <- grep("^sample.?id$", colnames(genotypedata), ignore.case = TRUE, value = TRUE)
  if (length(sid_col) == 0) stop("Input data must contain a 'Sample.ID' column.")
  genotypedata$Sample.ID <- as.character(genotypedata[[sid_col[1]]])

  long_data <- .pivot_long_genotypes(genotypedata, "allele")
  long_data <- long_data[long_data$marker_id %in% marker_info$marker_id, ]

  if (nrow(long_data) == 0) {
    message("No allele data found for the markers listed in 'marker_info'.")
    return(invisible(NULL))
  }

  if (!is.null(output_folder) && !dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }

  # `repeatlength` is the true repeat-unit size for microsatellites, so it
  # makes a sensible histogram bin width there. For msp_glurp markers it
  # instead encodes the family-clustering gap threshold (see
  # define_alleles_for_plotting()), which is typically far smaller than a
  # natural bin width -- using it directly produces dozens of near-empty
  # slivers that look like a blank plot. Only trust it for microsatellites;
  # everything else uses the data-driven fallback below.
  microsat_markers <- marker_info$marker_id[marker_info$binning_method == "microsatellite"]
  bin_lookup <- stats::setNames(
    suppressWarnings(as.numeric(marker_info$repeatlength[marker_info$marker_id %in% microsat_markers])),
    marker_info$marker_id[marker_info$marker_id %in% microsat_markers]
  )
  binning_lookup <- stats::setNames(marker_info$binning_method, marker_info$marker_id)

  list_markers <- unique(long_data$marker_id)
  n_markers    <- length(list_markers)

  base_colors <- RColorBrewer::brewer.pal(max(3, min(n_markers, 8)), "Set2")
  colors      <- rep_len(base_colors, n_markers)

  p_array <- Map(function(marker_name, color) {
    raw_vals <- long_data$allele[long_data$marker_id == marker_name]
    if (length(raw_vals) == 0) return(NULL)

    # AmpSeq markers: haplotype strings have no length to histogram, so show
    # a bar chart of counts per haplotype name instead.
    if (identical(binning_lookup[[marker_name]], "exact")) {
      haplo_counts <- as.data.frame(table(raw_vals), stringsAsFactors = FALSE)
      colnames(haplo_counts) <- c("haplotype", "count")
      haplo_counts <- haplo_counts[order(-haplo_counts$count), ]
      haplo_counts$haplotype <- factor(haplo_counts$haplotype, levels = haplo_counts$haplotype)
      haplo_counts$freq      <- haplo_counts$count / sum(haplo_counts$count)

      return(
        ggplot2::ggplot(haplo_counts, ggplot2::aes(x = .data$haplotype, y = .data$freq)) +
          ggplot2::geom_col(fill = color, color = "white") +
          ggplot2::scale_y_continuous(labels = scales::percent) +
          ggplot2::labs(
            title = paste0(marker_name, "\n(n=", length(raw_vals), ")"),
            x     = "Haplotype",
            y     = "Frequency"
          ) +
          ggplot2::theme_classic(base_size = 14) +
          ggplot2::theme(
            plot.title  = ggplot2::element_text(hjust = 0.5, size = 12),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
          )
      )
    }

    # Length-polymorphic markers: histogram of raw fragment sizes.
    vals <- suppressWarnings(as.numeric(raw_vals))
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) return(NULL)

    span <- diff(range(vals))
    bw <- binwidth
    if (is.null(bw)) {
      bw <- if (marker_name %in% names(bin_lookup)) bin_lookup[[marker_name]] else NA_real_
      # A few far-flung outlier alleles can stretch `span` well past the
      # bulk of the data, so a `repeatlength`-sized bin (correct as a repeat
      # unit) ends up chopping the range into dozens of near-empty slivers
      # that are invisible at typical plot sizes. Cap it: once repeatlength
      # would need more than max_bins bins to cover the observed range,
      # fall back to a bin width scaled to that range instead.
      max_bins <- 60
      if (is.na(bw) || bw <= 0 || (span > 0 && span / bw > max_bins)) {
        bw <- if (span > 0) span / 30 else 1
      }
    }

    ggplot2::ggplot(data.frame(allele = vals), ggplot2::aes(x = .data$allele)) +
      ggplot2::geom_histogram(
        ggplot2::aes(y = ggplot2::after_stat(.data$count / sum(.data$count))),
        binwidth = bw, fill = color, color = "white"
      ) +
      ggplot2::scale_y_continuous(labels = scales::percent) +
      ggplot2::labs(
        title = paste0(marker_name, "\n(n=", length(vals), ")"),
        x     = "Allele size (bp)",
        y     = "Frequency"
      ) +
      ggplot2::theme_classic(base_size = 14) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12))
  }, list_markers, colors)

  p_array <- Filter(Negate(is.null), p_array)
  if (length(p_array) == 0) {
    message("No plots were generated.")
    return(invisible(NULL))
  }

  num_cols     <- min(length(p_array), max_cols)
  num_rows     <- ceiling(length(p_array) / num_cols)
  final_figure <- ggpubr::ggarrange(plotlist = p_array, ncol = num_cols, nrow = num_rows)

  if (!is.null(output_folder)) {
    output_path <- file.path(output_folder, paste0(filename_prefix, ".png"))
    ggplot2::ggsave(output_path, plot = final_figure,
                    width  = num_cols * 4,
                    height = num_rows * 4,
                    limitsize = FALSE)
  }

  invisible(final_figure)
}

#' Combine MSP1/MSP2 family-variant results into one call per marker
#'
#' @description MSP1 and MSP2 are each genotyped across multiple family
#'   variants (e.g. K1/MAD20/RO33 for MSP1, 3D7/FC27/IC for MSP2). A patient
#'   counts as matching at "MSP1" (or "MSP2") if at least one of that
#'   marker's variants shows a match ("R") -- otherwise "NI", unless every
#'   variant is missing/errored (IND/ERR), in which case the combined call
#'   stays "IND" rather than being called a false "NI". The original
#'   per-variant columns are kept in the output alongside the new combined
#'   columns.
#'
#' @param match_counting_res Output of perform_match_counting().
#' @param msp1_cols Character vector of the MSP1 variant column names.
#' @param msp2_cols Character vector of the MSP2 variant column names.
#' @return match_counting_res with "msp1" and "msp2" columns added.
#' @noRd
combine_msp_variants <- function(match_counting_res, msp1_cols, msp2_cols) {
  collapse_variants <- function(row_vals) {
    if (any(row_vals == "R")) return("R")
    if (all(row_vals %in% c("IND", "ERR"))) return("IND")
    "NI"
  }
  
  match_counting_res$msp1 <- apply(match_counting_res[, msp1_cols, drop = FALSE], 1, collapse_variants)
  match_counting_res$msp2 <- apply(match_counting_res[, msp2_cols, drop = FALSE], 1, collapse_variants)
  
  match_counting_res   # original variant columns stay -- only msp1/msp2 get added
}

#' WHO classification rule -- N of M markers must match
#'
#' @description Applies the WHO rule to a given set of marker columns:
#'   WHO_loose requires at least round(loose_threshold * length(marker_cols))
#'   of them to match, WHO_strict requires all of them. Both are NA if fewer
#'   markers were even scored (R or NI) than the loose/strict threshold
#'   requires. Used both for full marker panels (marker_cols = every raw
#'   marker) and MSP1/MSP2 trio panels (marker_cols = the combined
#'   msp1/msp2/third columns, after combine_msp_variants() has already
#'   collapsed each family's variants).
#'
#' @param match_counting_res Output of perform_match_counting() (optionally
#'   with combine_msp_variants() already applied).
#' @param marker_cols Character vector of column names to apply the rule to.
#' @param loose_threshold Minimum proportion of markers required for
#'   WHO_loose. Default 0.70.
#' @return A list: table (match_counting_res with WHO_loose/WHO_strict
#'   columns added), who_loose_label (e.g. "Match counting 5/7"),
#'   who_strict_label (e.g. "Match counting 7/7").
#' @noRd
apply_who_rule <- function(match_counting_res, marker_cols, loose_threshold = 0.70) {
  n_markers    <- length(marker_cols)
  who_loose_n  <- round(loose_threshold * n_markers)
  who_strict_n <- n_markers
  
  cols <- match_counting_res[, marker_cols, drop = FALSE]
  n_matched <- rowSums(cols == "R")
  n_scored  <- rowSums(cols == "R" | cols == "NI")
  
  match_counting_res$WHO_loose  <- ifelse(n_scored < who_loose_n, NA, ifelse(n_matched >= who_loose_n, 1, 0))
  match_counting_res$WHO_strict <- ifelse(n_scored < who_strict_n, NA, ifelse(n_matched == who_strict_n, 1, 0))
  
  list(
    table = match_counting_res,
    who_loose_label  = paste0("Match counting ", who_loose_n, "/", n_markers),
    who_strict_label = paste0("Match counting ", who_strict_n, "/", n_markers)
  )
}

#' Build WHO classification columns using the appropriate rule for this panel
#'
#' @description Dispatches to the correct WHO rule based on panel
#'   composition: if MSP1/MSP2 family variants are detected (an MSP-trio
#'   panel, paired with either glurp or a single microsatellite marker),
#'   applies the classic 2/3 / 3/3 rule to the combined MSP1 + MSP2 + third
#'   marker. Otherwise (full microsatellite or full ampseq panel), applies
#'   the proportional 70%/100% rule to the whole panel.
#'
#' @param match_counting_res Output of perform_match_counting().
#' @param marker_info Marker metadata from import_data(), with binning_method.
#' @return A list: table, who_loose_label, who_strict_label -- see
#'   apply_who_trio_rule()/apply_who_proportional_rule() for details, since
#'   this just dispatches to one of the two.
#' @noRd
build_who_table <- function(match_counting_res, marker_info) {
  msp_variants <- detect_msp_variants(marker_info$marker_id)
  is_msp_trio  <- length(msp_variants$msp1) > 0 || length(msp_variants$msp2) > 0
  
  if (!is_msp_trio) {
    return(apply_who_rule(match_counting_res, unique(marker_info$marker_id)))
  }
  
  match_counting_res <- combine_msp_variants(match_counting_res, msp_variants$msp1, msp_variants$msp2)
  
  msp_glurp_leftover <- setdiff(marker_info$marker_id[marker_info$binning_method == "msp_glurp"],
                                c(msp_variants$msp1, msp_variants$msp2))
  microsat_markers   <- marker_info$marker_id[marker_info$binning_method == "microsatellite"]
  third_candidates   <- c(msp_glurp_leftover, microsat_markers)
  
  if (length(third_candidates) != 1) {
    stop("ERROR: Expected exactly one third marker (glurp or microsatellite) alongside MSP1/MSP2, found ",
         length(third_candidates), ": ", paste(third_candidates, collapse = ", "))
  }
  
  apply_who_rule(match_counting_res, c("msp1", "msp2", third_candidates))
}


#' Plot Match-Counting vs MalReBay Comparison Heatmap
#'
#' @description Visualizes the WHO match-counting rules (loose and strict)
#'   alongside MalReBay's posterior probability, one heatmap per site plus
#'   one combining all sites (patients grouped by site), from the
#'   bayesian_match_counting_comparison table after it's been passed through
#'   build_who_table().
#'
#' @param summary_results A list returned by \code{\link{summarise_results}},
#'   containing a \code{comparison} data frame (Bayesian probabilities plus
#'   traditional match-counting results).
#' @param marker_info Marker metadata from \code{\link{import_data}}, passed
#'   on to \code{build_who_table()} to derive the WHO_loose/WHO_strict columns.
#' @param output_folder A string specifying the directory where one PNG per
#'   site (\code{comparison_heatmap_<site>.png}) and one for all sites
#'   (\code{comparison_heatmap_all_sites.png}) will be saved. If \code{NULL},
#'   heatmaps are drawn on the current graphics device instead.
#' @param title_prefix Optional prefix before "Site <name>" / "All sites" in
#'   each plot's title.
#' @param verbose Logical. If \code{TRUE}, prints a message when each file is saved.
#' @return An invisible list of ggplot objects: \code{by_site} (a list named
#'   by site) and \code{all_sites}.
#'
#' @importFrom ggplot2 ggplot aes geom_tile facet_grid scale_fill_gradientn
#'   scale_y_discrete labs theme element_text element_blank element_rect margin
#' @export
plot_comparison_heatmap <- function(summary_results,
                                    marker_info,
                                    output_folder = NULL,
                                    title_prefix = "",
                                    verbose = TRUE) {
  comparison <- summary_results$comparison
  comparison$MalReBay <- comparison$MalReBay_probability
  who_result <- build_who_table(comparison, marker_info)

  comparison       <- who_result$table
  who_loose_label  <- who_result$who_loose_label
  who_strict_label <- who_result$who_strict_label

  block_levels <- c(who_loose_label, who_strict_label, "MalReBay")

  # by_site = TRUE facets rows by site (for the all-sites heatmap), so patients
  # stay grouped and labelled by site
  build_heatmap <- function(site_data, plot_title, by_site = FALSE) {
    # Site + ID as the row key: the same Sample.ID could appear at two sites
    row_key   <- paste(site_data$Site, site_data$Sample.ID, sep = "__")
    plot_data <- data.frame(
      row   = factor(rep(row_key, times = 3), levels = rev(unique(row_key))),
      Site  = factor(rep(site_data$Site, times = 3), levels = unique(site_data$Site)),
      block = factor(rep(block_levels, each = nrow(site_data)), levels = block_levels),
      value = c(site_data$WHO_loose, site_data$WHO_strict, site_data$MalReBay)
    )

    site_layers <- if (by_site) {
      # Sites stacked with no gap, a black line above every site but the first.
      # It's drawn on the top edge of the lower site's panel (y = n + 0.5, no
      # y expansion) because panels are drawn top to bottom: on the bottom edge
      # of the upper site, the next site's tiles would paint over half of it.
      # clip = "off" so the line isn't half-cut at the panel edge either.
      site_levels <- levels(plot_data$Site)
      n_per_site  <- table(factor(site_data$Site, levels = site_levels))
      separators  <- data.frame(Site = factor(site_levels[-1], levels = site_levels),
                                y    = as.numeric(n_per_site[-1]) + 0.5)
      list(
        ggplot2::facet_grid(Site ~ block, scales = "free_y", space = "free_y", switch = "y",
                            labeller = ggplot2::labeller(block = ggplot2::label_wrap_gen(width = 14))),
        ggplot2::geom_hline(data = separators, ggplot2::aes(yintercept = .data$y),
                            color = "black", linewidth = 0.35),
        ggplot2::scale_x_continuous(expand = c(0, 0)),
        ggplot2::scale_y_discrete(expand = c(0, 0)),
        ggplot2::coord_cartesian(clip = "off"),
        ggplot2::theme(panel.spacing.y = grid::unit(0, "mm"))
      )
    } else {
      # Wrap "Match counting N/M" onto two lines so it fits the narrow columns
      ggplot2::facet_grid(~block, labeller = ggplot2::label_wrap_gen(width = 14))
    }

    ggplot2::ggplot(plot_data, ggplot2::aes(x = 1, y = .data$row, fill = .data$value)) +
      # Thin white outline separates adjacent patients without the heavy look of grey borders
      ggplot2::geom_tile(color = "white", linewidth = 0.3) +
      ggplot2::scale_fill_gradientn(
        name    = "Outcome",
        # Pale blue (0, new infection) -> blue (0.5) -> red (1, recrudescence), RdBu hues
        colours = c("#EEF5FA", "#92C5DE", "#4393C3", "#F4A582", "#B2182B"),
        values  = c(0, 0.25, 0.5, 0.75, 1),
        limits  = c(0, 1), na.value = "grey80",
        breaks  = c(0, 0.25, 0.5, 0.75, 1),
        labels  = c("0 (New infection)", "0.25", "0.5", "0.75", "1 (Recrudescence)")
      ) +
      ggplot2::labs(x = NULL, y = "Samples", title = plot_title) +
      ggplot2::theme(
        axis.text.x     = ggplot2::element_blank(),
        axis.ticks.x    = ggplot2::element_blank(),
        axis.text.y     = ggplot2::element_blank(),
        axis.ticks.y    = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        panel.spacing   = grid::unit(5, "mm"),
        strip.background = ggplot2::element_rect(fill = "white", color = NA),
        strip.text      = ggplot2::element_text(size = 10),
        strip.text.y.left = ggplot2::element_text(size = 9, face = "bold", angle = 0, hjust = 1),
        strip.placement = "outside",
        plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12,
                                                 margin = ggplot2::margin(b = 8))
      ) +
      site_layers
  }

  # comparison already has one row per patient
  heatmap_data <- comparison

  if (!is.null(output_folder) && !dir.exists(output_folder))
    dir.create(output_folder, recursive = TRUE)

  sites   <- unique(heatmap_data$Site)
  by_site <- list()
  for (s in sites) {
    site_data    <- heatmap_data[heatmap_data$Site == s, ]
    p            <- build_heatmap(site_data, paste0(title_prefix, "Site ", s))
    by_site[[s]] <- p

    if (!is.null(output_folder)) {
      safe_site <- gsub("[^A-Za-z0-9_-]", "_", s)
      ggplot2::ggsave(
        file.path(output_folder, paste0("comparison_heatmap_", safe_site, ".png")),
        plot = p, width = 1800 / 300, height = 2200 / 300, units = "in", dpi = 300
      )
      if (verbose) message("INFO: Comparison heatmap saved for site: ", s)
    } else {
      print(p)
    }
  }

  # All sites together, patients grouped by site (in the order sites appear)
  all_data      <- heatmap_data[order(match(heatmap_data$Site, sites)), ]
  p_all         <- build_heatmap(all_data, paste0(title_prefix, "All sites"), by_site = TRUE)
  if (!is.null(output_folder)) {
    # Taller than the per-site plots so rows don't get squashed
    ggplot2::ggsave(
      file.path(output_folder, "comparison_heatmap_all_sites.png"),
      plot = p_all, width = 2100 / 300,
      height = max(2200 / 300, 0.12 * nrow(all_data) + 2), units = "in", dpi = 300, limitsize = FALSE
    )
    if (verbose) message("INFO: Comparison heatmap saved for all sites.")
  } else {
    print(p_all)
  }

  invisible(list(by_site = by_site, all_sites = p_all))
}
