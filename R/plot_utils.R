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

  # --- Convergence diagnostics ----------------------------------------------

  # Classical Gelman-Rubin (Gelman & Rubin 1992, Brooks & Gelman 1998)
  gelman_result <- NULL
  if (length(loglikelihood_mcmc) > 1) {
    gelman_result <- tryCatch(
      coda::gelman.diag(loglikelihood_mcmc, autoburnin = FALSE),
      error = function(e) NULL
    )
  }

  # Classical ESS (coda)
  ess_result <- tryCatch(
    coda::effectiveSize(loglikelihood_mcmc),
    error = function(e) NULL
  )

  # Rank-normalised R-hat + bulk and tail ESS (Vehtari et al. 2021)
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

  # Geweke Z-scores per chain (Geweke 1992)
  # Z within (-1.96, 1.96) indicates stationarity at the 95% level
  geweke_result <- tryCatch(
    sapply(loglikelihood_mcmc, function(ch) coda::geweke.diag(ch)$z),
    error = function(e) NULL
  )

  # --- Print diagnostics ----------------------------------------------------
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

    # ESS — classical, bulk, and tail
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
      cat("\nGeweke Z-scores per chain (Geweke 1992) — |Z| < 1.96 indicates stationarity:\n")
      for (i in seq_along(geweke_result)) {
        z <- geweke_result[i]
        cat(sprintf("  Chain %d: Z = %6.4f  %s\n", i, z,
                    ifelse(abs(z) < 1.96, "[PASS]", "[FAIL]")))
      }
    }
    cat(strrep("-", 50), "\n")
  }

  # --- Gelman-Rubin diagnostic plot ----------------------------------------
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

  # --- Traceplot -----------------------------------------------------------
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

  # --- Log-likelihood distribution -----------------------------------------
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

  # --- Autocorrelation plot ------------------------------------------------
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

  # --- Return all diagnostic values ----------------------------------------
  invisible(list(
    gelman       = gelman_result,     # classical Gelman-Rubin object
    ess          = ess_result,         # classical ESS from coda
    rhat_rank    = rhat_rank,          # rank-normalised R-hat (Vehtari et al. 2021)
    ess_bulk     = ess_bulk,           # bulk ESS (Vehtari et al. 2021)
    ess_tail     = ess_tail,           # tail ESS (Vehtari et al. 2021)
    geweke       = geweke_result       # Geweke Z-scores per chain
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
#' @return Invisibly returns the file path of the saved PNG. Returns
#'   \code{invisible(NULL)} if no probabilities are available to plot.
#'
#' @keywords internal
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
#' This function calculates the MOI (number of distinct alleles per marker per sample) 
#' and generates a faceted violin plot across different sites.
#'
#' @param genotypedata A data frame containing `Sample.ID`, `Site`, and marker columns.
#' @param marker_pattern A regex pattern to identify marker columns. 
#'   Defaults to standard allele suffix patterns.
#' @param output_folder Path to the directory where the plot will be saved. 
#'   If NULL, the plot is not saved to disk.
#' @param filename The name for the output PDF file.
#' 
#' @return A ggplot object. The underlying MOI data is attached as an 
#'   attribute "moi_data".
#' @export
plot_moi <- function(genotypedata,
                     marker_pattern = "(_allele_\\d+|_\\d+)$",
                     output_folder = NULL,
                     filename = "moi_per_marker_by_site.png") {
  
  # Validating if the site column is available
  if (!"Site" %in% colnames(genotypedata)) {
    stop("Input 'genotypedata' must contain a column named 'Site'.")
  }
  
  # Calculating the MOI for each day of observation per sample
  genotypedata$Sample.ID <- as.character(genotypedata$Sample.ID)
  
  # Identify marker columns
  marker_cols <- grep(marker_pattern, colnames(genotypedata), value = TRUE)
  if (length(marker_cols) == 0) {
    warning("No marker columns found matching the pattern.")
    return(NULL)
  }
  
  moi_data <- genotypedata %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(marker_cols),
      names_to = "marker_replicate",
      values_to = "allele",
      values_drop_na = TRUE
    ) %>%
    dplyr::mutate(marker_id = gsub(marker_pattern, "", .data$marker_replicate)) %>%
    dplyr::group_by(.data$Sample.ID, .data$Site, .data$marker_id) %>%
    dplyr::summarise(MOI = dplyr::n_distinct(.data$allele), .groups = "drop")
  
  # Fill in zeroes for samples/markers with no data
  all_sample_sites <- dplyr::distinct(genotypedata, .data$Sample.ID, .data$Site)
  all_markers <- unique(gsub(marker_pattern, "", marker_cols))
  complete_grid <- tidyr::crossing(all_sample_sites, marker_id = all_markers)
  
  moi_data <- dplyr::left_join(complete_grid, moi_data, by = c("Sample.ID", "Site", "marker_id")) %>%
    dplyr::mutate(MOI = tidyr::replace_na(.data$MOI, 0))

  if (nrow(moi_data) == 0) {
    message("MOI data is empty. Skipping plot generation.")
    return(NULL)
  }

  # Prepare plotting labels
  label_data <- moi_data %>%
    dplyr::group_by(.data$Site, .data$marker_id) %>%
    dplyr::summarise(mean_moi = mean(.data$MOI), .groups = "drop")
  
  site_max_moi <- moi_data %>%
    dplyr::group_by(.data$Site) %>%
    dplyr::summarise(label_y_pos = max(.data$MOI) + 0.5, .groups = "drop")
  
  label_data <- dplyr::left_join(label_data, site_max_moi, by = "Site")
  moi_data$marker_id <- factor(moi_data$marker_id, levels = unique(label_data$marker_id))
  
  # Generating the MOI plots
  moi_plot <- ggplot2::ggplot(moi_data, ggplot2::aes(x = .data$marker_id, y = .data$MOI, fill = .data$marker_id, color = .data$marker_id)) + 
    ggplot2::geom_jitter(width = 0.15, height = 0.1, alpha = 0.3) + 
    ggplot2::geom_violin(alpha = 0.4, trim = FALSE) + 
    ggplot2::geom_text(
      data = label_data,
      ggplot2::aes(x = .data$marker_id, y = .data$label_y_pos, label = sprintf("%.2f", .data$mean_moi)),
      size = 3.5, vjust = 0, color = "black"
    ) +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::scale_color_brewer(palette = "Set2") +
    ggplot2::facet_grid(.data$Site ~ ., scales = "free_y") +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks()) +
    ggplot2::coord_cartesian(ylim = c(-0.5, NA)) +
    ggplot2::labs(
      title = "Multiplicity of Infection (MOI) by Marker and Site",
      x = "Marker", y = "MOI (Number of Alleles)"
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      legend.position = "none", 
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.background = ggplot2::element_rect(fill = "grey90", color = "black"),
      strip.text.y = ggplot2::element_text(angle = 0, face = "bold")
    )
  
  # Save and return
  if (!is.null(output_folder)) {
    if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
    output_path <- file.path(output_folder, filename)
    n_sites <- length(unique(moi_data$Site))
    ggplot2::ggsave(output_path, plot = moi_plot, width = 12, height = 3 + (1.5 * n_sites))
    message(paste("Successfully saved MOI plot to:", output_path))
  }
  
  # Attach data as attribute so it's not "lost"
  attr(moi_plot, "moi_data") <- moi_data
  
  return(moi_plot)
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
  alleles_definitions_bin <- list()
  data_marker_columns <- grep("(_allele_\\d+|_\\d+)$", colnames(genotypedata), value = TRUE)
  available_base_markers <- unique(gsub("(_allele_\\d+|_\\d+)$", "", data_marker_columns))
  
  for (locus_name in available_base_markers) {
    locus_marker_info <- marker_info[marker_info$marker_id == locus_name, ]
    if (nrow(locus_marker_info) == 0) next
    
    binning_method <- locus_marker_info$binning_method[1]
    
    if (binning_method == "microsatellite") {
      locus_cols_pattern <- paste0("^", locus_name, "(_allele_|_)\\d+")
      locus_cols_logical <- grepl(locus_cols_pattern, colnames(genotypedata))
      raw_alleles_vector <- unlist(genotypedata[, locus_cols_logical])
      unique_alleles <- sort(unique(raw_alleles_vector[!is.na(raw_alleles_vector)]))
      
      if (length(unique_alleles) > 0) {
        repeat_length <- locus_marker_info$repeatlength[1]
        bins <- ceiling((unique_alleles - unique_alleles[1] + 1) / repeat_length)
        alleles_definitions_bin[[locus_name]] <- data.frame(min = tapply(unique_alleles, bins, min),
                                                            max = tapply(unique_alleles, bins, max))
      }
    } else if (binning_method == "msp_glurp") {
      locus_cols_pattern <- paste0("^", locus_name, "(_allele_|_)\\d+")
      locus_cols_logical <- grepl(locus_cols_pattern, colnames(genotypedata))
      raw_alleles_vector <- unlist(genotypedata[, locus_cols_logical])
      unique_alleles <- sort(unique(raw_alleles_vector[!is.na(raw_alleles_vector)]))
      
      if (length(unique_alleles) > 0) {
        gap_threshold <- locus_marker_info$repeatlength[1]
        breaks <- c(0, which(diff(unique_alleles) > gap_threshold), length(unique_alleles))
        bins <- findInterval(seq_along(unique_alleles), breaks)
        alleles_definitions_bin[[locus_name]] <- data.frame(min = tapply(unique_alleles, bins, min),
                                                            max = tapply(unique_alleles, bins, max))
      }
    }
  }
  return(alleles_definitions_bin)
}


#' Process Data for Pie Chart
#'
#' @description Internal helper to prepare data for `plot_pie_chart`. It renames
#'   columns, creates labels, and orders the data.
#'
#' @param plot_data_df A data frame with Haplotype, Amount, and Frequency.
#' @param data_type The type of data ("lp" or "ampseq").
#' @return A processed data frame ready for plotting.
#' @noRd

# Merged data processing function for plotting
process_pie_data <- function(plot_data_df, data_type) {
  colnames(plot_data_df) <- c("Haplotype", "Amount", "Frequency")
  plot_data_df <- plot_data_df[!is.na(plot_data_df$Haplotype), ]
  
  if (data_type == "lp") {
    plot_data_df$Haplotype <- factor(plot_data_df$Haplotype, levels = sort(as.numeric(as.character(plot_data_df$Haplotype))))
  } else {
    plot_data_df$Haplotype <- as.factor(plot_data_df$Haplotype)
  }
  
  plot_data_df$Label <- paste0(round(plot_data_df$Frequency * 100), "%")
  plot_data_df <- plot_data_df[order(plot_data_df$Amount, decreasing = TRUE), ]
  
  if (nrow(plot_data_df) > 4) {
    plot_data_df[5:nrow(plot_data_df), "Label"] <- ""
  }
  return(plot_data_df)
}



#' Plot a Pie Chart for Haplotype/Allele Diversity
#'
#' @description Internal helper function to create a single pie chart.
#'
#' @param data_df A dataframe with Haplotype/Allele, Amount, and Frequency.
#' @param color_marker The base color for the chart palette.
#' @param marker_name The name of the genetic marker.
#' @param total_n The total number of infections for this marker.
#' @param data_type The type of data ("length_polymorphic" or "ampseq").
#' @return A ggplot object representing the pie chart.
#'
#' @importFrom dplyr mutate lead if_else
#' @importFrom ggplot2 ggplot aes geom_col coord_polar scale_fill_manual geom_text theme_void theme element_text ggtitle
#' @importFrom forcats fct_inorder
#' @importFrom grDevices colorRampPalette
#' @noRd
plot_pie_chart <- function(data_df, color_marker, marker_name, total_n, data_type) {
  if (nrow(data_df) == 0) return(NULL)
  
  # Process data for plotting
  plot_data_df <- process_pie_data(data_df, data_type)
  title_marker <- paste0(marker_name, "\n(n=", total_n, ")")
  
  plot_data_df <- plot_data_df %>%
    dplyr::mutate(
      csum = rev(cumsum(rev(.data$Amount))),
      pos = .data$Amount/2 + dplyr::lead(.data$csum, 1),
      pos = dplyr::if_else(is.na(.data$pos), .data$Amount/2, .data$pos)
    )
  
  colfunc <- grDevices::colorRampPalette(c(color_marker, "white"))
  
  pie_chart <- ggplot2::ggplot(plot_data_df, ggplot2::aes(x = "", y = .data$Amount, fill = forcats::fct_inorder(.data$Haplotype))) +
    ggplot2::geom_col(width = 1, color = "white") +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::scale_fill_manual(values = colfunc(nrow(plot_data_df))) +
    ggplot2::geom_text(ggplot2::aes(y = .data$pos, label = .data$Label), size = 4, color = "black") +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12)
    ) +
    ggplot2::ggtitle(title_marker)
  
  return(pie_chart)
}


#' Generate and Save Diversity Pie Charts
#'
#' @description This function creates pie charts to visualize allele or haplotype
#' diversity at Day 0 and at the day of recurrence. It saves the combined plot
#' as an image.
#'
#' @param genotypedata A dataframe containing the genotyping data.
#' @param data_type A string, either "length_polymorphic" or "ampseq".
#' @param marker_info A dataframe with information about the markers, required
#'   if `data_type` is "length_polymorphic".
#' @param output_folder The path to the directory where the output PDF will be saved.
#' @param filename_prefix A string prefix for the output PDF filename.
#' @return Invisibly returns the final combined ggplot object.
#'
#' @importFrom dplyr %>% filter mutate select distinct group_by summarize summarise matches left_join bind_rows n
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggsave
#' @importFrom ggpubr ggarrange annotate_figure text_grob
#' @importFrom RColorBrewer brewer.pal
#' @export

plot_markers_diversity <- function(genotypedata,
                                     data_type,
                                     marker_info = NULL,
                                     output_folder,
                                     filename_prefix = "diversity") {
  
  if (!data_type %in% c("length_polymorphic", "ampseq")) {
    stop("Invalid data_type specified. Please use 'length_polymorphic' or 'ampseq'.")
  }
  
  if (data_type == "length_polymorphic" && is.null(marker_info)) {
    stop("marker_info must be provided for data_type = 'length_polymorphic'.")
  }
  
  if (!"Sample.ID" %in% colnames(genotypedata)) {
    # Check for case-insensitive match
    sid_col <- grep("^sample.?id$", colnames(genotypedata), ignore.case = TRUE, value = TRUE)
    if(length(sid_col) > 0) {
      genotypedata$Sample.ID <- genotypedata[[sid_col[1]]]
    } else {
      stop("Input data must contain a 'Sample.ID' column.")
    }
  }
  
  if (data_type == "ampseq") {
    if (!any(grepl(" Day ", genotypedata$Sample.ID))) {
      genotypedata$Sample.ID <- gsub("D0$", " Day 0", genotypedata$Sample.ID)
      genotypedata$Sample.ID <- gsub("D[1-9][0-9]*$", " Day Failure", genotypedata$Sample.ID)
    }
  }
  
  day0_data <- genotypedata[grepl(" Day 0", genotypedata$Sample.ID), ]
  recurrence_data <- genotypedata[!grepl(" Day 0", genotypedata$Sample.ID), ]
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  figure_day0 <- NULL
  figure_recurrence <- NULL
  
  # Variables to track how many rows of pie charts are in each section
  rows_d0 <- 0
  rows_rec <- 0
  
  # day 0
  if (nrow(day0_data) > 0) {
    if (data_type == "length_polymorphic") {
      alleles_definitions_bin <- define_alleles_for_plotting(genotypedata = genotypedata, marker_info = marker_info)
      long_data_d0 <- day0_data %>%
        tidyr::pivot_longer(cols = matches("_allele_\\d+$|_\\d+$"), names_to = "marker_replicate", values_to = "allele_length") %>%
        dplyr::filter(!is.na(.data$allele_length)) %>%
        dplyr::mutate(marker_id = gsub("_allele_\\d+$|_\\d+$", "", .data$marker_replicate)) %>%
        dplyr::select("Sample.ID", "marker_id", "allele_length") %>%
        dplyr::distinct()
      
      if(nrow(long_data_d0) > 0) {
        binned_alleles_list_d0 <- list()
        for (marker in unique(long_data_d0$marker_id)) {
          if (!marker %in% names(alleles_definitions_bin)) { next }
          marker_data <- long_data_d0 %>% filter(.data$marker_id == marker)
          current_definitions <- alleles_definitions_bin[[marker]]
          bin_centers <- rowMeans(current_definitions, na.rm = TRUE)
          marker_data$true_alleles <- sapply(marker_data$allele_length, function(x) bin_centers[which.min(abs(x - bin_centers))])
          binned_alleles_list_d0[[marker]] <- marker_data
        }
        final_alleles_d0 <- bind_rows(binned_alleles_list_d0) %>% select("Sample.ID", "marker_id", "true_alleles") %>% distinct()
        frequency_data_d0 <- final_alleles_d0 %>% group_by(.data$marker_id, .data$true_alleles) %>% summarize(Amount = n(), .groups = 'drop')
        total_infections_d0 <- frequency_data_d0 %>% group_by(.data$marker_id) %>% summarise(TotalInfections = sum(.data$Amount), .groups = 'drop')
        frequency_data_d0 <- left_join(frequency_data_d0, total_infections_d0, by = "marker_id") %>% mutate(Frequency = .data$Amount / .data$TotalInfections)
        
        list_markers_d0 <- unique(frequency_data_d0$marker_id)
        p_array_d0 <- vector('list', length(list_markers_d0))
        colors <- RColorBrewer::brewer.pal(max(3, length(list_markers_d0)), "Set2")
        if (length(list_markers_d0) > length(colors)) colors <- rep(colors, length.out = length(list_markers_d0))
        
        for (i in seq_along(list_markers_d0)) {
          marker_name <- list_markers_d0[i]
          plot_data <- frequency_data_d0 %>% filter(.data$marker_id == marker_name)
          p_array_d0[[i]] <- plot_pie_chart(plot_data %>% select("true_alleles", "Amount", "Frequency"), colors[i], marker_name, plot_data$TotalInfections[1], "length_polymorphic")
        }
        
        p_array_d0 <- p_array_d0[!sapply(p_array_d0, is.null)]
        if(length(p_array_d0) > 0) {
          num_cols_d0 <- min(length(p_array_d0), 7)
          rows_d0 <- ceiling(length(p_array_d0) / num_cols_d0) 
          figure_day0 <- ggpubr::ggarrange(plotlist = p_array_d0, ncol = num_cols_d0, nrow = rows_d0)
          figure_day0 <- ggpubr::annotate_figure(figure_day0, left = ggpubr::text_grob("Day 0", rot = 90, size = 18))
        }
      }
    } else if (data_type == "ampseq") {
      long_data_d0 <- day0_data %>%
        tidyr::pivot_longer(cols = matches("_allele_\\d+$|_\\d+$"), names_to = "marker_replicate", values_to = "haplotype") %>%
        dplyr::filter(!is.na(.data$haplotype) & .data$haplotype != "") %>%
        dplyr::mutate(marker_id = gsub("_allele_\\d+$|_\\d+$", "", .data$marker_replicate)) %>%
        dplyr::select("Sample.ID", "marker_id", "haplotype") %>%
        dplyr::distinct()
      
      if(nrow(long_data_d0) > 0) {
        frequency_data_d0 <- long_data_d0 %>% group_by(.data$marker_id, .data$haplotype) %>% summarize(Amount = n(), .groups = 'drop')
        total_infections_d0 <- frequency_data_d0 %>% group_by(.data$marker_id) %>% summarise(TotalInfections = sum(.data$Amount), .groups = 'drop')
        frequency_data_d0 <- left_join(frequency_data_d0, total_infections_d0, by = "marker_id") %>% mutate(Frequency = .data$Amount / .data$TotalInfections)
        
        list_markers_d0 <- unique(frequency_data_d0$marker_id)
        p_array_d0 <- vector('list', length(list_markers_d0))
        colors <- RColorBrewer::brewer.pal(max(3, length(list_markers_d0)), "Set2")
        if (length(list_markers_d0) > length(colors)) colors <- rep(colors, length.out = length(list_markers_d0))
        
        for (i in seq_along(list_markers_d0)) {
          marker_name <- list_markers_d0[i]
          plot_data <- frequency_data_d0 %>% filter(.data$marker_id == marker_name)
          p_array_d0[[i]] <- plot_pie_chart(plot_data %>% select("haplotype", "Amount", "Frequency"), colors[i], marker_name, plot_data$TotalInfections[1], "ampseq")
        }
        
        p_array_d0 <- p_array_d0[!sapply(p_array_d0, is.null)]
        if(length(p_array_d0) > 0) {
          num_cols_d0 <- min(length(p_array_d0), 7)
          rows_d0 <- ceiling(length(p_array_d0) / num_cols_d0) 
          figure_day0 <- ggpubr::ggarrange(plotlist = p_array_d0, ncol = num_cols_d0, nrow = rows_d0)
          figure_day0 <- ggpubr::annotate_figure(figure_day0, left = ggpubr::text_grob("Day 0", rot = 90, size = 18))
        }
      }
    }
  } else { message("No Day 0 data found. Skipping.") }
  
  # Day of recurrence
  if (nrow(recurrence_data) > 0) {
    if (data_type == "length_polymorphic") {
      if (!exists("alleles_definitions_bin")) {
        alleles_definitions_bin <- define_alleles_for_plotting(genotypedata = genotypedata, marker_info = marker_info)
      }
      long_data_rec <- recurrence_data %>%
        tidyr::pivot_longer(cols = matches("_allele_\\d+$|_\\d+$"), names_to = "marker_replicate", values_to = "allele_length") %>%
        dplyr::filter(!is.na(.data$allele_length)) %>%
        dplyr::mutate(marker_id = gsub("_allele_\\d+$|_\\d+$", "", .data$marker_replicate)) %>%
        dplyr::select("Sample.ID", "marker_id", "allele_length") %>%
        dplyr::distinct()
      
      if(nrow(long_data_rec) > 0) {
        binned_alleles_list_rec <- list()
        for (marker in unique(long_data_rec$marker_id)) {
          if (!marker %in% names(alleles_definitions_bin)) { next }
          marker_data <- long_data_rec %>% filter(.data$marker_id == marker)
          current_definitions <- alleles_definitions_bin[[marker]]
          bin_centers <- rowMeans(current_definitions, na.rm = TRUE)
          marker_data$true_alleles <- sapply(marker_data$allele_length, function(x) bin_centers[which.min(abs(x - bin_centers))])
          binned_alleles_list_rec[[marker]] <- marker_data
        }
        final_alleles_rec <- bind_rows(binned_alleles_list_rec) %>% select("Sample.ID", "marker_id", "true_alleles") %>% distinct()
        
        frequency_data_rec <- final_alleles_rec %>% group_by(.data$marker_id, .data$true_alleles) %>% summarize(Amount = n(), .groups = 'drop')
        total_infections_rec <- frequency_data_rec %>% group_by(.data$marker_id) %>% summarise(TotalInfections = sum(.data$Amount), .groups = 'drop')
        frequency_data_rec <- left_join(frequency_data_rec, total_infections_rec, by = "marker_id") %>% mutate(Frequency = .data$Amount / .data$TotalInfections)
        
        list_markers_rec <- unique(frequency_data_rec$marker_id)
        p_array_rec <- vector('list', length(list_markers_rec))
        colors <- RColorBrewer::brewer.pal(max(3, length(list_markers_rec)), "Set2")
        if (length(list_markers_rec) > length(colors)) colors <- rep(colors, length.out = length(list_markers_rec))
        
        for (i in seq_along(list_markers_rec)) {
          marker_name <- list_markers_rec[i]
          plot_data <- frequency_data_rec %>% filter(.data$marker_id == marker_name)
          p_array_rec[[i]] <- plot_pie_chart(plot_data %>% select("true_alleles", "Amount", "Frequency"), colors[i], marker_name, plot_data$TotalInfections[1], "length_polymorphic")
        }
        
        p_array_rec <- p_array_rec[!sapply(p_array_rec, is.null)]
        if(length(p_array_rec) > 0) {
          num_cols_rec <- min(length(p_array_rec), 7)
          rows_rec <- ceiling(length(p_array_rec) / num_cols_rec)
          figure_recurrence <- ggpubr::ggarrange(plotlist = p_array_rec, ncol = num_cols_rec, nrow = rows_rec)
          figure_recurrence <- ggpubr::annotate_figure(figure_recurrence, left = ggpubr::text_grob("Day of Recurrence", rot = 90, size = 18))
        }
      }
    } else if (data_type == "ampseq") {
      long_data_rec <- recurrence_data %>%
        tidyr::pivot_longer(cols = matches("_allele_\\d+$|_\\d+$"), names_to = "marker_replicate", values_to = "haplotype") %>%
        dplyr::filter(!is.na(.data$haplotype) & .data$haplotype != "") %>%
        dplyr::mutate(marker_id = gsub("_allele_\\d+$|_\\d+$", "", .data$marker_replicate)) %>%
        dplyr::select("Sample.ID", "marker_id", "haplotype") %>%
        dplyr::distinct()
      
      if(nrow(long_data_rec) > 0) {
        frequency_data_rec <- long_data_rec %>% group_by(.data$marker_id, .data$haplotype) %>% summarize(Amount = n(), .groups = 'drop')
        total_infections_rec <- frequency_data_rec %>% group_by(.data$marker_id) %>% summarise(TotalInfections = sum(.data$Amount), .groups = 'drop')
        frequency_data_rec <- left_join(frequency_data_rec, total_infections_rec, by = "marker_id") %>% mutate(Frequency = .data$Amount / .data$TotalInfections)
        
        list_markers_rec <- unique(frequency_data_rec$marker_id)
        p_array_rec <- vector('list', length(list_markers_rec))
        colors <- RColorBrewer::brewer.pal(max(3, length(list_markers_rec)), "Set2")
        if (length(list_markers_rec) > length(colors)) colors <- rep(colors, length.out = length(list_markers_rec))
        
        for (i in seq_along(list_markers_rec)) {
          marker_name <- list_markers_rec[i]
          plot_data <- frequency_data_rec %>% filter(.data$marker_id == marker_name)
          p_array_rec[[i]] <- plot_pie_chart(plot_data %>% select("haplotype", "Amount", "Frequency"), colors[i], marker_name, plot_data$TotalInfections[1], "ampseq")
        }
        
        p_array_rec <- p_array_rec[!sapply(p_array_rec, is.null)]
        if(length(p_array_rec) > 0) {
          num_cols_rec <- min(length(p_array_rec), 7)
          rows_rec <- ceiling(length(p_array_rec) / num_cols_rec)
          figure_recurrence <- ggpubr::ggarrange(plotlist = p_array_rec, ncol = num_cols_rec, nrow = rows_rec)
          figure_recurrence <- ggpubr::annotate_figure(figure_recurrence, left = ggpubr::text_grob("Day of Recurrence", rot = 90, size = 18))
        }
      }
    }
  } else { message("No recurrence data found. Skipping.") }
  
  final_plot_list <- list(figure_day0, figure_recurrence)
  final_plot_list <- final_plot_list[!sapply(final_plot_list, is.null)]
  plot_heights <- c()
  if (!is.null(figure_day0)) plot_heights <- c(plot_heights, rows_d0)
  if (!is.null(figure_recurrence)) plot_heights <- c(plot_heights, rows_rec)
  
  if (length(final_plot_list) > 0) {
    final_figure <- ggpubr::ggarrange(plotlist = final_plot_list, 
                                      ncol = 1, 
                                      nrow = length(final_plot_list),
                                      heights = plot_heights)
    
    # output_path <- file.path(output_folder, paste0(filename_prefix, "_", data_type, "_comparison.pdf"))
    output_path <- file.path(output_folder, paste0(filename_prefix, "_", data_type, "_comparison.png"))
    calc_height <- max(4, sum(plot_heights) * 3.5)
    
    ggsave(output_path, plot = final_figure, width = 24, height = calc_height, limitsize = FALSE)
    message(paste("Successfully saved combined diversity plots to:", output_path))
    invisible(final_figure)
  } else {
    message("No plots were generated.")
    invisible(NULL)
  }
}