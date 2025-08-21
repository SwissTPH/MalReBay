#' Generate and save an allele frequency distribution plot
#'
#' This function creates a faceted bar plot showing the frequency of each allele
#' length for Day 0 vs. Day of Failure samples. It can generate plots for a
#' single site or an overall plot for all sites combined.
#'
#' @param raw_data_df A data frame containing raw genotyping data. For a site-specific
#'   plot, this should be filtered for one site. For an overall plot, it should
#'   contain data from all sites.
#' @param site_name The name of the site for the plot title and filename.
#'   If `NULL` (the default), an "Overall" plot is generated.
#' @param output_folder The directory where the plot PDF will be saved.
#' @param plots_per_row The number of marker plots to display per row.
#' @param plot_width The width of the output PDF file.
#' @param plot_height The height of the output PDF file.
#'
#' @importFrom dplyr mutate across all_of if_else group_by summarise n ungroup
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_bar facet_wrap scale_fill_manual labs theme_bw theme element_rect element_text ggsave
#' @importFrom stats frequency
#' @importFrom rlang .data
#'
generate_allele_frequency_plot <- function(raw_data_df,
                                           site_name = NULL,
                                           output_folder,
                                           plots_per_row = 4,
                                           plot_width = 12,
                                           plot_height = 10) {
  
  # Set up plot title and filename based on whether it's overall or site-specific
  if (is.null(site_name)) {
    message("Generating OVERALL raw allele frequency plot for all sites.")
    plot_title <- "Overall Raw Allele Frequency Distribution (All Sites)"
    plot_filename_part <- "Overall"
  } else {
    message(paste("Generating raw allele frequency plot for site:", site_name))
    plot_title <- paste("Raw Allele Frequency Distribution for Site:", site_name)
    plot_filename_part <- gsub("[^a-zA-Z0-9_]", "_", site_name) # Sanitize site name for filename
  }
  
  # Data Preparation
  marker_cols <- grep("_", colnames(raw_data_df), value = TRUE)
  if (length(marker_cols) == 0) {
    warning("No marker columns found for plotting. Skipping plot generation.")
    return(invisible(NULL))
  }
  
  plot_data_prepared <- raw_data_df %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(marker_cols), as.numeric)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(marker_cols),
      names_to = "allele_name",
      values_to = "raw_length",
      values_drop_na = TRUE
    ) %>%
    dplyr::mutate(
      marker_id = sub("_.*", "", .data$allele_name),
      Day_Type = dplyr::if_else(grepl("Day 0", .data$Sample.ID), "Day 0", "Day Failure"), 
      Day_Type = factor(.data$Day_Type, levels = c("Day 0", "Day Failure"))
    )
  
  # Calculate Frequencies
  frequency_data <- plot_data_prepared %>%
    dplyr::group_by(.data$marker_id, .data$Day_Type, .data$raw_length) %>%
    dplyr::summarise(count = n(), .groups = 'drop') %>%
    dplyr::group_by(.data$marker_id, .data$Day_Type) %>% 
    dplyr::mutate(frequency = .data$count / sum(.data$count)) %>%
    dplyr::ungroup()
  
  # Create the Plot
  allele_plot <- ggplot2::ggplot(data = frequency_data, ggplot2::aes(x = .data$raw_length, y = .data$frequency, fill = .data$Day_Type)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge", color = "black", width = 1) +
    ggplot2::facet_wrap(~ marker_id, ncol = plots_per_row, scales = "free_x") + 
    ggplot2::scale_fill_manual(name = "Sample day", values = c("Day 0" = "black", "Day of Failure" = "white")) +
    ggplot2::labs(
      title = plot_title,
      subtitle = "Comparison of Day 0 and Day of Failure Samples",
      x = "Allele Length (bp)",
      y = "Frequency"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      strip.background = ggplot2::element_rect(fill = "grey90", color = "black"),
      strip.text = ggplot2::element_text(face = "bold", size = 10),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  
  # 5. Save the Plot
  plot_filename <- file.path(output_folder, paste0("allele_freq_plot_", plot_filename_part, ".pdf"))
  ggplot2::ggsave(
    filename = plot_filename,
    plot = allele_plot,
    width = plot_width,
    height = plot_height,
    dpi = 300
  )
  message("Plot saved to: ", plot_filename)
  
  return(invisible(allele_plot))
}