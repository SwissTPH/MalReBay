#' Calculate Marker-Specific MOI
#'
#' @param genotypedata A data frame with `Sample.ID`, `Site`, and marker columns.
#' @param marker_pattern A regex pattern to identify marker columns.
#' @return A data frame summarizing MOI for each sample and marker.
#' @export

calculate_marker_moi <- function(genotypedata, marker_pattern = "_\\d+$") {
  
  if (!"Site" %in% colnames(genotypedata)) {
    stop("Input 'genotypedata' must contain a column named 'Site'.")
  }
  
  genotypedata$Sample.ID <- as.character(genotypedata$Sample.ID)
  
  moi_data <- genotypedata %>%
    tidyr::pivot_longer(
      cols = dplyr::matches(marker_pattern),
      names_to = "marker_replicate",
      values_to = "allele",
      values_drop_na = TRUE
    ) %>%
    dplyr::mutate(marker_id = gsub(marker_pattern, "", .data$marker_replicate)) %>%
    dplyr::group_by(.data$Sample.ID, .data$Site, .data$marker_id) %>%
    dplyr::summarise(MOI = dplyr::n_distinct(.data$allele), .groups = "drop")
  
  all_sample_sites <- dplyr::distinct(genotypedata, .data$Sample.ID, .data$Site)
  all_markers <- unique(gsub(marker_pattern, "", colnames(genotypedata)[grep(marker_pattern, colnames(genotypedata))]))
  
  complete_grid <- tidyr::crossing(all_sample_sites, marker_id = all_markers)
  
  moi_data_complete <- dplyr::left_join(complete_grid, moi_data, by = c("Sample.ID", "Site", "marker_id")) %>%
    dplyr::mutate(MOI = tidyr::replace_na(.data$MOI, 0))
  
  return(moi_data_complete)
}


#' Generate and Save an MOI Plot
#'
#' @param moi_data A data frame of MOI values, typically from `calculate_marker_moi`.
#' @param output_folder Path to the directory for saving the plot.
#' @param filename The name for the output PDF file.
#' @return The ggplot object of the MOI plot.
#' @export
generate_marker_moi_plot <- function(moi_data,
                                     output_folder,
                                     filename = "moi_per_marker_by_site.pdf") {
  
  if (nrow(moi_data) == 0) {
    message("MOI data is empty. Skipping plot generation.")
    return(NULL)
  }
  
  label_data <- moi_data %>%
    dplyr::group_by(.data$Site, .data$marker_id) %>%
    dplyr::summarise(mean_moi = mean(.data$MOI), .groups = "drop")
  
  site_max_moi <- moi_data %>%
    dplyr::group_by(.data$Site) %>%
    dplyr::summarise(label_y_pos = max(.data$MOI) + 0.5, .groups = "drop")
  
  label_data <- dplyr::left_join(label_data, site_max_moi, by = "Site")
  
  moi_data$marker_id <- factor(moi_data$marker_id, levels = unique(label_data$marker_id))
  
  moi_plot <- ggplot2::ggplot(moi_data, ggplot2::aes(x = .data$marker_id, y = .data$MOI, fill = .data$marker_id, color = .data$marker_id)) + 
    ggplot2::geom_jitter(width = 0.15, height = 0.1, alpha = 0.3) + 
    ggplot2::geom_violin(alpha = 0.4, trim = FALSE, quantiles = 0.5) + 
    ggplot2::geom_text(
      data = label_data,
      ggplot2::aes(x = .data$marker_id, y = .data$label_y_pos, label = sprintf("%.2f", .data$mean_moi)),
      size = 3.5,
      vjust = 0,
      color = "black"
    ) +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::scale_color_brewer(palette = "Set2") +
    ggplot2::facet_grid(
      .data$Site ~ ., 
      scales = "free_y"
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks(n = max(moi_data$MOI, 1))
    ) +
    ggplot2::coord_cartesian(ylim = c(-0.5, NA)) +
    ggplot2::labs(
      title = "Multiplicity of Infection (MOI) by Marker and Site",
      x = "Marker",
      y = "MOI (Number of Alleles)"
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      legend.position = "none", 
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.background = ggplot2::element_rect(fill = "grey90", color = "black"),
      strip.text.y = ggplot2::element_text(angle = 0, face = "bold")
    )
  
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  output_path <- file.path(output_folder, filename)
  ggplot2::ggsave(output_path, plot = moi_plot, width = 12, height = 3 + 1.5 * length(unique(moi_data$Site)))
  message(paste("Successfully saved faceted MOI plot to:", output_path))
  
  return(moi_plot)
}