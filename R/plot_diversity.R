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
  data_marker_columns <- grep("_\\d+$", colnames(genotypedata), value = TRUE)
  available_base_markers <- unique(gsub("_\\d+$", "", data_marker_columns))
  
  for (locus_name in available_base_markers) {
    locus_marker_info <- marker_info[marker_info$marker_id == locus_name, ]
    if (nrow(locus_marker_info) == 0) next
    
    binning_method <- locus_marker_info$binning_method[1]
    
    if (binning_method == "microsatellite") {
      locus_cols_pattern <- paste0("^", locus_name, "_")
      locus_cols_logical <- grepl(locus_cols_pattern, colnames(genotypedata))
      raw_alleles_vector <- unlist(genotypedata[, locus_cols_logical])
      unique_alleles <- sort(unique(raw_alleles_vector[!is.na(raw_alleles_vector)]))
      
      if (length(unique_alleles) > 0) {
        repeat_length <- locus_marker_info$repeatlength[1]
        bins <- ceiling((unique_alleles - unique_alleles[1] + 1) / repeat_length)
        alleles_definitions_bin[[locus_name]] <- data.frame(min = tapply(unique_alleles, bins, min),
                                                        max = tapply(unique_alleles, bins, max))
      }
    } else if (binning_method == "cluster") {
      locus_cols_pattern <- paste0("^", locus_name, "_")
      locus_cols_logical <- grepl(locus_cols_pattern, colnames(genotypedata))
      raw_alleles_vector <- unlist(genotypedata[, locus_cols_logical])
      unique_alleles <- sort(unique(raw_alleles_vector[!is.na(raw_alleles_vector)]))
      
      if (length(unique_alleles) > 0) {
        gap_threshold <- locus_marker_info$cluster_gap_threshold[1]
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
    # Smaller labels for better proportion with pie charts
    ggplot2::geom_text(ggplot2::aes(y = .data$pos, label = .data$Label), size = 3, color = "black") +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold")
    ) +
    ggplot2::ggtitle(title_marker)
  
  return(pie_chart)
}


#' Generate and Save Diversity Pie Charts
#'
#' @description This function creates pie charts to visualize allele or haplotype
#' diversity at Day 0 and at the day of recurrence. It saves the combined plot
#' as a PDF file.
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

generate_diversity_plots <- function(genotypedata,
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
  
  # day 0
  if (nrow(day0_data) > 0) {
    if (data_type == "length_polymorphic") {
      alleles_definitions_bin <- define_alleles_for_plotting(genotypedata = genotypedata, marker_info = marker_info)
      long_data_d0 <- day0_data %>%
        pivot_longer(cols = matches("_\\d+$"), names_to = "marker_replicate", values_to = "allele_length") %>%
        filter(!is.na(.data$allele_length)) %>%
        mutate(marker_id = gsub("_\\d+$", "", .data$marker_replicate)) %>%
        select(.data$Sample.ID, .data$marker_id, .data$allele_length) %>%
        distinct()
      
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
        final_alleles_d0 <- bind_rows(binned_alleles_list_d0) %>% select(.data$Sample.ID, .data$marker_id, .data$true_alleles) %>% distinct()
        
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
          p_array_d0[[i]] <- plot_pie_chart(plot_data %>% select(.data$true_alleles, .data$Amount, .data$Frequency), colors[i], marker_name, plot_data$TotalInfections[1], "length_polymorphic")
        }
        
        p_array_d0 <- p_array_d0[!sapply(p_array_d0, is.null)]
        if(length(p_array_d0) > 0) {
          num_cols_d0 <- min(length(p_array_d0), 7)
          figure_day0 <- ggpubr::ggarrange(plotlist = p_array_d0, ncol = num_cols_d0, nrow = ceiling(length(p_array_d0) / num_cols_d0))
          figure_day0 <- ggpubr::annotate_figure(figure_day0, top = ggpubr::text_grob("Allele Diversity at Day 0", color = "black", face = "bold", size = 20))
        }
      }
    } else if (data_type == "ampseq") {
      long_data_d0 <- day0_data %>%
        pivot_longer(cols = matches("_allele_\\d+$"), names_to = "marker_replicate", values_to = "haplotype") %>%
        filter(!is.na(.data$haplotype) & .data$haplotype != "") %>%
        mutate(marker_id = gsub("_allele_\\d+$", "", .data$marker_replicate)) %>%
        select(.data$Sample.ID, .data$marker_id, .data$haplotype) %>%
        distinct()
      
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
          p_array_d0[[i]] <- plot_pie_chart(plot_data %>% select(.data$haplotype, .data$Amount, .data$Frequency), colors[i], marker_name, plot_data$TotalInfections[1], "ampseq")
        }
        
        p_array_d0 <- p_array_d0[!sapply(p_array_d0, is.null)]
        if(length(p_array_d0) > 0) {
          num_cols_d0 <- min(length(p_array_d0), 7)
          figure_day0 <- ggpubr::ggarrange(plotlist = p_array_d0, ncol = num_cols_d0, nrow = ceiling(length(p_array_d0) / num_cols_d0))
          figure_day0 <- ggpubr::annotate_figure(figure_day0, top = ggpubr::text_grob("Haplotype Diversity at Day 0", color = "black", face = "bold", size = 20))
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
        pivot_longer(cols = matches("_\\d+$"), names_to = "marker_replicate", values_to = "allele_length") %>%
        filter(!is.na(.data$allele_length)) %>%
        mutate(marker_id = gsub("_\\d+$", "", .data$marker_replicate)) %>%
        select(.data$Sample.ID, .data$marker_id, .data$allele_length) %>%
        distinct()
      
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
        final_alleles_rec <- bind_rows(binned_alleles_list_rec) %>% select(.data$Sample.ID, .data$marker_id, .data$true_alleles) %>% distinct()
        
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
          p_array_rec[[i]] <- plot_pie_chart(plot_data %>% select(.data$true_alleles, .data$Amount, .data$Frequency), colors[i], marker_name, plot_data$TotalInfections[1], "length_polymorphic")
        }
        
        p_array_rec <- p_array_rec[!sapply(p_array_rec, is.null)]
        if(length(p_array_rec) > 0) {
          num_cols_rec <- min(length(p_array_rec), 7)
          figure_recurrence <- ggpubr::ggarrange(plotlist = p_array_rec, ncol = num_cols_rec, nrow = ceiling(length(p_array_rec) / num_cols_rec))
          figure_recurrence <- ggpubr::annotate_figure(figure_recurrence, top = ggpubr::text_grob("Allele Diversity at Recurrence", color = "black", face = "bold", size = 20))
        }
      }
    } else if (data_type == "ampseq") {
      long_data_rec <- recurrence_data %>%
        pivot_longer(cols = matches("_allele_\\d+$"), names_to = "marker_replicate", values_to = "haplotype") %>%
        filter(!is.na(.data$haplotype) & .data$haplotype != "") %>%
        mutate(marker_id = gsub("_allele_\\d+$", "", .data$marker_replicate)) %>%
        select(.data$Sample.ID, .data$marker_id, .data$haplotype) %>%
        distinct()
      
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
          p_array_rec[[i]] <- plot_pie_chart(plot_data %>% select(.data$haplotype, .data$Amount, .data$Frequency), colors[i], marker_name, plot_data$TotalInfections[1], "ampseq")
        }
        
        p_array_rec <- p_array_rec[!sapply(p_array_rec, is.null)]
        if(length(p_array_rec) > 0) {
          num_cols_rec <- min(length(p_array_rec), 7)
          figure_recurrence <- ggpubr::ggarrange(plotlist = p_array_rec, ncol = num_cols_rec, nrow = ceiling(length(p_array_rec) / num_cols_rec))
          figure_recurrence <- ggpubr::annotate_figure(figure_recurrence, top = ggpubr::text_grob("Haplotype Diversity at Recurrence", color = "black", face = "bold", size = 20))
        }
      }
    }
  } else { message("No recurrence data found. Skipping.") }
  
  final_plot_list <- list(figure_day0, figure_recurrence)
  final_plot_list <- final_plot_list[!sapply(final_plot_list, is.null)]
  
  if (length(final_plot_list) > 0) {
    final_figure <- ggpubr::ggarrange(plotlist = final_plot_list, ncol = 1, nrow = length(final_plot_list))
    output_path <- file.path(output_folder, paste0(filename_prefix, "_", data_type, "_comparison.pdf"))
    ggsave(output_path, plot = final_figure, width = 24, height = 5 * length(final_plot_list), limitsize = FALSE)
    message(paste("Successfully saved combined diversity plots to:", output_path))
    invisible(final_figure)
  } else {
    message("No plots were generated.")
    invisible(NULL)
  }
}