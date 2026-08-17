#' Plot MCMC Likelihood Diagnostics
#'
#' @description
#' Calculates standard MCMC convergence diagnostics and generates four diagnostic
#' plots: Gelman-Rubin, Traceplot, Log-Posterior Histogram, and Autocorrelation.
#' When \code{save_plot = FALSE} (interactive use), all four panels are arranged
#' in a single 2x2 grid by default.
#'
#' @param all_chains_loglikelihood A list where each element is a numeric vector
#'   representing the log-likelihood history of one MCMC chain.
#' @param site_name A character string for labeling plots.
#' @param save_plot A logical. If `TRUE`, plots are saved as PNG files.
#' @param output_folder A character string specifying the path to save plots.
#' @param verbose A logical. If `TRUE`, prints diagnostic summaries.
#' @param stan_fit An optional \code{CmdStanMCMC} object (from cmdstanr) for additional diagnostics.
#' @param combine_plots A logical. If \code{TRUE}, all four panels are drawn in a
#'   single 2x2 grid on the current graphics device. Defaults to \code{TRUE} when
#'   \code{save_plot = FALSE}. Ignored when \code{save_plot = TRUE}.
#'
#' @return An invisible list containing the calculated diagnostic results.
#'
#' @examples
#' \dontrun{
#'   chains <- list(
#'     c(-10.2, -9.8, -10.5, -9.9, -10.1),
#'     c(-9.9,  -10.3, -10.0, -9.7, -10.4)
#'   )
#'   plot_likelihood_diagnostics(
#'     all_chains_loglikelihood = chains,
#'     site_name  = "TestSite",
#'     save_plot  = FALSE,
#'     verbose    = FALSE
#'   )
#' }
#'
#' @importFrom coda as.mcmc.list mcmc varnames gelman.diag effectiveSize gelman.plot
#' @importFrom grDevices png dev.off rainbow n2mfrow
#' @importFrom graphics matplot legend title hist par plot text abline lines layout plot.new
#' @importFrom stats acf
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

  loglikelihood_mcmc <- tryCatch({
    mlist  <- lapply(clean_chains, coda::mcmc)
    mclist <- coda::as.mcmc.list(mlist)
    coda::varnames(mclist) <- " "
    mclist
  }, error = function(e) NULL)

  if (is.null(loglikelihood_mcmc)) return(invisible(NULL))

  # When drawing to screen, arrange all four panels in a 2x2 grid
  if (combine_plots && !save_plot) {
    old_outer_par <- graphics::par(mfrow = c(2, 2), mar = c(4.1, 4.1, 3.5, 1.1))
    on.exit(graphics::par(old_outer_par), add = TRUE)
  }

  # Convergence diagnostics
  gelman_result <- NULL
  if (length(loglikelihood_mcmc) > 1) {
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

  # When no output_folder is provided, return diagnostics without attempting
  # to draw plots — avoids "figure margins too large" in small RStudio windows
  if (!save_plot) {
    return(invisible(list(
      gelman    = gelman_result,
      ess       = ess_result,
      rhat_rank = rhat_rank,
      ess_bulk  = ess_bulk,
      ess_tail  = ess_tail,
      geweke    = geweke_result
    )))
  }

  # Gelman-Rubin plot
  if (save_plot) {
    grDevices::png(file.path(site_dir,
                              paste0(safe_site_name, "_gelman_rubin.png")),
                   width = 1000, height = 700, res = 120)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  old_par <- graphics::par(mar = c(4.1, 4.1, 3.5, 1.1))
  on.exit(graphics::par(old_par), add = TRUE)
  if (length(loglikelihood_mcmc) > 1) {
    tryCatch(
      coda::gelman.plot(loglikelihood_mcmc, autoburnin = FALSE,
                        main = paste(site_name, "- Gelman-Rubin Diagnostic"),
                        col  = c("black", "indianred")),
      error = function(e) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "Gelman plot failed\n(near-constant chains)",
                       col = "red")
      }
    )
  } else {
    graphics::plot.new()
    graphics::text(0.5, 0.5, "Gelman plot not applicable\n(requires > 1 chain)")
  }
  if (save_plot) { grDevices::dev.off(); on.exit(NULL, add = FALSE) }
  graphics::par(old_par); on.exit(NULL, add = FALSE)

  # Traceplot
  if (save_plot) {
    grDevices::png(file.path(site_dir,
                              paste0(safe_site_name, "_traceplot.png")),
                   width = 1000, height = 600, res = 120)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  old_par <- graphics::par(mar = c(4.1, 4.1, 3.5, 1.1))
  on.exit(graphics::par(old_par), add = TRUE)
  colors <- grDevices::rainbow(length(loglikelihood_mcmc))
  graphics::matplot(
    do.call(cbind, lapply(loglikelihood_mcmc, as.numeric)),
    type = "l", lty = 1, col = colors,
    main = paste(site_name, "- Traceplot"),
    ylab = "Log-Posterior", xlab = "Iterations"
  )
  graphics::legend("topright",
                   legend = paste("Chain", seq_along(loglikelihood_mcmc)),
                   col = colors, lty = 1, bty = "n", cex = 0.8)
  if (save_plot) { grDevices::dev.off(); on.exit(NULL, add = FALSE) }
  graphics::par(old_par); on.exit(NULL, add = FALSE)

  # Log-posterior distribution
  if (save_plot) {
    grDevices::png(file.path(site_dir,
                              paste0(safe_site_name, "_log_likelihood_distribution.png")),
                   width = 1000, height = 600, res = 120)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  old_par <- graphics::par(mar = c(4.1, 4.1, 3.5, 1.1))
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::hist(unlist(clean_chains), breaks = 40,
                 main = paste(site_name, "- Log-Posterior Distribution"),
                 xlab = "Log-Posterior", col = "steelblue", border = "white")
  if (save_plot) { grDevices::dev.off(); on.exit(NULL, add = FALSE) }
  graphics::par(old_par); on.exit(NULL, add = FALSE)

  # Autocorrelation
  if (save_plot) {
    grDevices::png(file.path(site_dir,
                              paste0(safe_site_name, "_autocorrelation.png")),
                   width = 1200, height = 1000, res = 120)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  if (combine_plots && !save_plot) {
    # Single panel in the 2x2 grid: show chain 1 only
    stats::acf(clean_chains[[1]], lag.max = 50, main = "")
    graphics::title(main = "Autocorrelation - Chain 1", cex.main = 1.2)
  } else {
    plot_dims <- grDevices::n2mfrow(length(clean_chains))
    old_par   <- graphics::par(mfrow = plot_dims, mar = c(4, 4, 3, 1))
    on.exit(graphics::par(old_par), add = TRUE)
    for (i in seq_along(clean_chains)) {
      stats::acf(clean_chains[[i]], lag.max = 50, main = "")
      graphics::title(main = paste("Autocorrelation - Chain", i), cex.main = 1.2)
    }
    if (save_plot) { grDevices::dev.off(); on.exit(NULL, add = FALSE) }
    graphics::par(old_par); on.exit(NULL, add = FALSE)
  }

  invisible(list(
    gelman    = gelman_result,
    ess       = ess_result,
    rhat_rank = rhat_rank,
    ess_bulk  = ess_bulk,
    ess_tail  = ess_tail,
    geweke    = geweke_result
  ))
}

#' Plot Posterior Probability Histogram
#'
#' Creates and saves a histogram of posterior probabilities of recrudescence
#' for all patients. This is an internal function called automatically by
#' \code{\link{MalReBay}}.
#'
#' @param summary_results A list returned by \code{\link{summarise_results}},
#'   containing a \code{posterior_probabilities} data frame with a
#'   \code{Probability} column.
#' @param output_folder A string specifying the directory where the histogram
#'   PNG will be saved. Defaults to \code{"results"}.
#' @param verbose Logical. If \code{TRUE}, prints a message when the file is
#'   saved. Defaults to \code{TRUE}.
#'
#' @return When \code{output_folder} is \code{NULL}, displays the histogram
#'   on the current graphics device and returns \code{invisible(NULL)}.
#'   When \code{output_folder} is provided, saves a PNG and returns the file
#'   path invisibly. Returns \code{invisible(NULL)} if no probabilities are
#'   available.
#'
#' @examples
#' summary_results <- list(
#'   posterior_probabilities = data.frame(
#'     Patient_ID  = c("P1", "P2", "P3"),
#'     Probability = c(0.12, 0.87, 0.45)
#'   )
#' )
#' plot_probability_histogram(summary_results)
#'
#' @export
plot_probability_histogram <- function(summary_results, output_folder = NULL, verbose = TRUE) {

  posterior_probabilities <- summary_results$posterior_probabilities

  if (is.null(posterior_probabilities) || nrow(posterior_probabilities) == 0) {
    warning("No posterior probabilities to plot.")
    return(invisible(NULL))
  }

  probs <- as.numeric(as.character(posterior_probabilities$Probability))

  if (!is.null(output_folder)) {
    if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
    path <- file.path(output_folder, "recrudescence_probability_histogram.png")
    grDevices::png(path, width = 8, height = 6, units = "in", res = 300)
    graphics::hist(
      probs,
      breaks = seq(0, 1, by = 0.05),
      col    = "skyblue",
      main   = "Posterior Probability Distribution",
      xlab   = "Probability of Recrudescence",
      ylab   = "Number of Patients"
    )
    grDevices::dev.off()
    if (verbose) message("INFO: Probability histogram saved to: ", output_folder)
    return(invisible(path))
  }

  graphics::hist(
    probs,
    breaks = seq(0, 1, by = 0.05),
    col    = "skyblue",
    main   = "Posterior Probability Distribution",
    xlab   = "Probability of Recrudescence",
    ylab   = "Number of Patients"
  )
  invisible(NULL)
}

#' Plot Multiplicity of Infection (MOI)
#'
#' This function calculates the MOI (number of distinct alleles per marker per
#' sample) and generates one violin plot per site for improved readability.
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

  # Output folder
  if (!is.null(output_folder) && !dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }

  # Build a plot per site 
  sites <- unique(moi_data$Site)

  plots <- lapply(stats::setNames(sites, sites), function(site) {

    site_moi    <- dplyr::filter(moi_data,   .data$Site == site)
    site_labels <- dplyr::filter(label_data, .data$Site == site)

    site_moi$marker_id <- factor(site_moi$marker_id, levels = marker_levels)

    p <- ggplot2::ggplot(
      site_moi,
      ggplot2::aes(
        x     = .data$marker_id,
        y     = .data$MOI,
        fill  = .data$marker_id,
        color = .data$marker_id
      )
    ) +
      ggplot2::geom_jitter(width = 0.15, height = 0.1, alpha = 0.3) +
      ggplot2::geom_violin(alpha = 0.4, trim = FALSE) +
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
#'   columns added), who_loose_label (e.g. "WHO 5/7"), who_strict_label
#'   (e.g. "WHO 7/7").
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
    who_loose_label  = paste0("WHO ", who_loose_n, "/", n_markers),
    who_strict_label = paste0("WHO ", who_strict_n, "/", n_markers)
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
#'   alongside MalReBay's posterior probability, one heatmap per site, from
#'   the bayesian_match_counting_comparison table after it's been passed
#'   through build_who_table().
#'
#' @param comparison The comparison table (e.g. summary_results$comparison,
#'   after build_who_table() has added WHO_loose/WHO_strict, and with a
#'   "MalReBay" column -- rename from "Probability" if needed).
#' @param who_loose_label Column title for the WHO_loose block, from
#'   build_who_table()'s who_loose_label.
#' @param who_strict_label Column title for the WHO_strict block, from
#'   build_who_table()'s who_strict_label.
#' @param title_prefix Optional prefix before the site name in each plot's title.
#' @return Invisibly NULL; draws one heatmap per site as a side effect.
#' @export
plot_comparison_heatmap <- function(summary_results, 
                                    marker_info, 
                                    output_folder = NULL,    
                                    title_prefix = "",
                                    verbose = TRUE) {
  comparison <- summary_results$comparison
  comparison$MalReBay <- comparison$Probability
  who_result <- build_who_table(comparison, marker_info)
  
  comparison       <- who_result$table
  who_loose_label  <- who_result$who_loose_label
  who_strict_label <- who_result$who_strict_label
  
  col_fun <- circlize::colorRamp2(c(0, 0.25, 0.5, 0.75, 1),
                                  c("#67A9CF", "#ade8f4", "#F7F7F7", "#F4A582", "#D6604D"))
  
  block_defaults <- list(
    col = col_fun, na_col = "grey80",
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_column_names = FALSE,
    column_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    rect_gp = grid::gpar(col = "grey70", lwd = 0.4)
  )
  
  build_site_heatmap <- function(site_data) {
    mat_who_loose  <- as.matrix(site_data[, "WHO_loose", drop = FALSE])
    mat_who_strict <- as.matrix(site_data[, "WHO_strict", drop = FALSE])
    mat_mr         <- as.matrix(site_data[, "MalReBay", drop = FALSE])
    rownames(mat_who_loose) <- rownames(mat_who_strict) <- rownames(mat_mr) <- site_data$Sample.ID
    
    ht_who_loose <- do.call(ComplexHeatmap::Heatmap, c(list(
      matrix = mat_who_loose, name = "Probability", column_title = who_loose_label,
      show_row_names = FALSE, row_names_gp = grid::gpar(fontsize = 7), row_names_side = "left",
      heatmap_legend_param = list(
        title = "Outcome",
        at = c(0, 0.25, 0.5, 0.75, 1),
        labels = c("0 (NI)", "0.25", "0.5", "0.75", "1 (R)")
      )
    ), block_defaults))
    
    ht_who_strict <- do.call(ComplexHeatmap::Heatmap, c(list(
      matrix = mat_who_strict, name = "WHO_strict_leg", column_title = who_strict_label,
      show_row_names = FALSE, show_heatmap_legend = FALSE
    ), block_defaults))
    
    ht_mr <- do.call(ComplexHeatmap::Heatmap, c(list(
      matrix = mat_mr, name = "MalReBay_leg", column_title = "MalReBay",
      show_row_names = FALSE, show_heatmap_legend = FALSE
    ), block_defaults))
    
    ht_who_loose + ht_who_strict + ht_mr
  }
  
  heatmap_data <- comparison[grepl(" Day 0$", comparison$Sample.ID), ]
  heatmap_data$Sample.ID <- trimws(gsub(" Day 0$", "", heatmap_data$Sample.ID))
  
  for (s in unique(heatmap_data$Site)) {
    site_data <- heatmap_data[heatmap_data$Site == s, ]
    ht        <- build_site_heatmap(site_data)
    
    if (!is.null(output_folder)) {
      if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
      safe_site <- gsub("[^A-Za-z0-9_-]", "_", s)
      grDevices::png(
        file.path(output_folder, paste0("comparison_heatmap_", safe_site, ".png")),
        width = 1800, height = 2200, res = 300
      )
      ComplexHeatmap::draw(
        ht, column_title = paste0(title_prefix, s),
        column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
        column_title_side = "bottom", row_title = "Samples",
        row_title_side = "left", gap = grid::unit(5, "mm")
      )
      grDevices::dev.off()
      if (verbose) message("INFO: Comparison heatmap saved for site: ", s)
    } else {
      ComplexHeatmap::draw(
        ht, column_title = paste0(title_prefix, s),
        column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
        column_title_side = "bottom", row_title = "Samples",
        row_title_side = "left", gap = grid::unit(5, "mm")
      )
    }
  }
  
  invisible(NULL)
}
