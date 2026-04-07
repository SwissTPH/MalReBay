#' Plot MCMC Likelihood Diagnostics
#'
#' @description
#' Calculates standard MCMC convergence diagnostics and generates diagnostic plots
#' in a multi-stage process: 1) Gelman-Rubin, 2) Traceplot & Histogram, 3) All ACF plots.
#'
#' @param all_chains_loglikelihood A list where each element is a numeric vector
#'   representing the log-likelihood history of one MCMC chain.
#' @param site_name A character string for labeling plots.
#' @param save_plot A logical. If `TRUE`, plots are saved as PNG files.
#' @param output_folder A character string specifying the path to save plots.
#' @param verbose A logical. If `TRUE`, prints diagnostic summaries.
#'
#' @return An invisible list containing the calculated diagnostic results.
#'
#' @importFrom coda as.mcmc.list mcmc varnames gelman.diag effectiveSize gelman.plot
#' @importFrom grDevices png dev.off rainbow n2mfrow
#' @importFrom graphics matplot legend title hist par plot text abline lines layout plot.new
#' @importFrom stats acf
#'
#' @export
plot_likelihood_diagnostics <- function(all_chains_loglikelihood,
                                        site_name,
                                        save_plot     = TRUE,
                                        output_folder = NULL,
                                        verbose       = TRUE) {
  
  if (save_plot && is.null(output_folder)) {
    stop("`output_folder` must be provided when `save_plot = TRUE`.", call. = FALSE)
  }
  
  safe_site_name <- gsub(" ", "_", site_name)
  
  site_dir <- NULL
  if (save_plot) {
    site_dir <- file.path(output_folder, "convergence_diagnosis", safe_site_name)
    if (dir.exists(site_dir)) unlink(site_dir, recursive = TRUE)
    dir.create(site_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # --- Data cleaning --------------------------------------------------------
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
  
  # Gelman-Rubin
  gelman_result <- NULL
  if (length(loglikelihood_mcmc) > 1) {
    gelman_result <- tryCatch(
      coda::gelman.diag(loglikelihood_mcmc, autoburnin = FALSE),
      error = function(e) NULL
    )
  }
  
  # ESS 
  ess_result <- tryCatch(
    coda::effectiveSize(loglikelihood_mcmc),
    error = function(e) NULL
  )
  
  # Rank-normalised R-hat + bulk and tail ESS 
  draws_matrix <- tryCatch(
    do.call(cbind, lapply(clean_chains, as.numeric)),
    error = function(e) NULL
  )
  
  rhat_rank <- tryCatch(
    posterior::rhat(draws_matrix),
    error = function(e) NA_real_
  )
  ess_bulk <- tryCatch(
    posterior::ess_bulk(draws_matrix),
    error = function(e) NA_real_
  )
  ess_tail <- tryCatch(
    posterior::ess_tail(draws_matrix),
    error = function(e) NA_real_
  )
  
  # Geweke Z-scores per chain
  geweke_result <- tryCatch(
    sapply(loglikelihood_mcmc, function(ch) coda::geweke.diag(ch)$z),
    error = function(e) NULL
  )
  
  if (verbose) {
    cat("Convergence Diagnostics for:", site_name, "\n")
    cat(strrep("-", 50), "\n")
    
    # Classical Gelman-Rubin
    if (!is.null(gelman_result)) {
      cat("Classical Gelman-Rubin R-hat (Gelman & Rubin 1992):\n")
      print(gelman_result)
    } else {
      cat("Classical Gelman-Rubin: not available (near-constant chains or single chain)\n")
    }
    
    # Rank-normalised R-hat
    if (!is.na(rhat_rank)) {
      cat("\nRank-normalised R-hat (Vehtari et al. 2021):", round(rhat_rank, 4),
          ifelse(rhat_rank < 1.01, "  [PASS < 1.01]", "  [FAIL >= 1.01]"), "\n")
    }
    
    # ESS -- classical, bulk, and tail
    if (!is.null(ess_result)) {
      cat("\nClassical ESS (coda):", round(as.numeric(ess_result), 1), "\n")
    }
    if (!is.na(ess_bulk)) {
      cat("Bulk ESS (Vehtari et al. 2021):", round(ess_bulk, 1),
          ifelse(ess_bulk > 400, "  [PASS > 400]", "  [FAIL <= 400]"), "\n")
    }
    if (!is.na(ess_tail)) {
      cat("Tail ESS (Vehtari et al. 2021):", round(ess_tail, 1),
          ifelse(ess_tail > 400, "  [PASS > 400]", "  [FAIL <= 400]"), "\n")
    }
    
    # Geweke
    if (!is.null(geweke_result)) {
      cat("\nGeweke Z-scores per chain (Geweke 1992) -- |Z| < 1.96 indicates stationarity:\n")
      for (i in seq_along(geweke_result)) {
        z <- geweke_result[i]
        cat(sprintf("  Chain %d: Z = %6.4f  %s\n", i, z,
                    ifelse(abs(z) < 1.96, "[PASS]", "[FAIL]")))
      }
    }
    cat(strrep("-", 50), "\n")
  }
  
  # Gelman-Rubin diagnostic plot
  if (save_plot) {
    grDevices::png(
      file.path(site_dir, paste0(safe_site_name, "_gelman_rubin.png")),
      width = 1000, height = 700, res = 120
    )
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  old_par <- par(mar = c(4.1, 4.1, 3.5, 1.1))
  on.exit(par(old_par), add = TRUE)
  if (length(loglikelihood_mcmc) > 1) {
    tryCatch(
      coda::gelman.plot(
        loglikelihood_mcmc,
        autoburnin = FALSE,
        main       = paste(site_name, "- Gelman-Rubin Diagnostic"),
        col        = c("black", "indianred")
      ),
      error = function(e) {
        plot.new()
        text(0.5, 0.5, "Gelman plot failed\n(near-constant chains)", col = "red")
      }
    )
  } else {
    plot.new()
    text(0.5, 0.5, "Gelman plot not applicable\n(requires > 1 chain)")
  }
  if (save_plot) { grDevices::dev.off(); on.exit(NULL, add = FALSE) }
  par(old_par); on.exit(NULL, add = FALSE)
  
  # Traceplot
  if (save_plot) {
    grDevices::png(
      file.path(site_dir, paste0(safe_site_name, "_traceplot.png")),
      width = 1000, height = 600, res = 120
    )
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  old_par <- par(mar = c(4.1, 4.1, 3.5, 1.1))
  on.exit(par(old_par), add = TRUE)
  colors <- grDevices::rainbow(length(loglikelihood_mcmc))
  matplot(
    do.call(cbind, lapply(loglikelihood_mcmc, as.numeric)),
    type = "l", lty = 1, col = colors,
    main = paste(site_name, "- Traceplot"),
    ylab = "Log-Likelihood",
    xlab = "Iterations"
  )
  legend("topright", legend = paste("Chain", seq_along(loglikelihood_mcmc)),
         col = colors, lty = 1, bty = "n", cex = 0.8)
  if (save_plot) { grDevices::dev.off(); on.exit(NULL, add = FALSE) }
  par(old_par); on.exit(NULL, add = FALSE)
  
  # Log-likelihood distribution
  if (save_plot) {
    grDevices::png(
      file.path(site_dir, paste0(safe_site_name, "_log_likelihood_distribution.png")),
      width = 1000, height = 600, res = 120
    )
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  old_par <- par(mar = c(4.1, 4.1, 3.5, 1.1))
  on.exit(par(old_par), add = TRUE)
  hist(unlist(clean_chains), breaks = 40,
       main = paste(site_name, "- Log-Likelihood Distribution"),
       xlab = "Log-Likelihood", col = "steelblue", border = "white")
  if (save_plot) { grDevices::dev.off(); on.exit(NULL, add = FALSE) }
  par(old_par); on.exit(NULL, add = FALSE)
  
  # Autocorrelation plot 
  if (save_plot) {
    grDevices::png(
      file.path(site_dir, paste0(safe_site_name, "_autocorrelation.png")),
      width = 1200, height = 1000, res = 120
    )
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  plot_dims <- grDevices::n2mfrow(length(clean_chains))
  old_par   <- par(mfrow = plot_dims, mar = c(4, 4, 3, 1))
  on.exit(par(old_par), add = TRUE)
  for (i in seq_along(clean_chains)) {
    stats::acf(clean_chains[[i]], lag.max = 50, main = "")
    graphics::title(main = paste("Autocorrelation - Chain", i), cex.main = 1.2)
  }
  if (save_plot) { grDevices::dev.off(); on.exit(NULL, add = FALSE) }
  par(old_par); on.exit(NULL, add = FALSE)
  
  
  invisible(list(
    gelman       = gelman_result,     
    ess          = ess_result,         
    rhat_rank    = rhat_rank,          
    ess_bulk     = ess_bulk,           
    ess_tail     = ess_tail,           
    geweke       = geweke_result       
  ))
}

#' Plot Posterior Probability Histogram
#'
#' Creates and saves a histogram of posterior probabilities of recrudescence
#' for all patients. Also called automatically by \code{\link{MalReBay}}.
#'
#' @param summary_results A list returned by \code{\link{MalReBay}},
#'   containing a \code{posterior_probabilities} data frame with a
#'   \code{Probability} column.
#' @param output_folder A string specifying the directory where the histogram
#'   PNG will be saved. Defaults to \code{"results"}.
#' @param verbose Logical. If \code{TRUE}, prints a message when the file is
#'   saved. Defaults to \code{TRUE}.
#'
#' @return Invisibly returns the file path of the saved PNG. Returns
#'   \code{invisible(NULL)} if no probabilities are available to plot.
#'
#' @export
plot_probability_histogram <- function(summary_results, output_folder = "results", verbose = TRUE) {
  
  posterior_probabilities <- summary_results$posterior_probabilities
  
  if (is.null(posterior_probabilities) || nrow(posterior_probabilities) == 0) {
    warning("No posterior probabilities to plot.")
    return(invisible(NULL))
  }
  
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  
  path <- file.path(output_folder, "recrudescence_probability_histogram.png")
  grDevices::png(path, width = 8, height = 6, units = "in", res = 300)
  graphics::hist(
    as.numeric(as.character(posterior_probabilities$Probability)),
    breaks = seq(0, 1, by = 0.05),
    col    = "skyblue",
    main   = "Posterior Probability Distribution",
    xlab   = "Probability of Recrudescence",
    ylab   = "Number of Patients"
  )
  grDevices::dev.off()
  
  if (verbose) message("INFO: Probability histogram saved to: ", output_folder)
  invisible(path)
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
#' @importFrom ggplot2 ggplot aes geom_col coord_polar scale_fill_manual geom_text
#'   theme_void theme element_text ggtitle
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
#' @importFrom dplyr %>% filter mutate select distinct group_by summarise
#'   left_join all_of n
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggsave
#' @importFrom ggpubr ggarrange
#' @importFrom RColorBrewer brewer.pal
#' @importFrom purrr map_dfr
#' @export
plot_markers_diversity <- function(genotypedata,
                                   data_type,
                                   marker_info     = NULL,
                                   output_folder,
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
  
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  
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
  
  output_path <- file.path(output_folder, paste0(filename_prefix, "_", data_type, "_comparison.png"))
  ggplot2::ggsave(output_path, plot = final_figure,
                  width  = num_cols * 4,
                  height = num_rows * 4,
                  limitsize = FALSE)
  
  invisible(final_figure)
}