#' Import and process genotyping data
#'
#' @param filepath Path to the input Excel file.
#' @return A list containing processed data frames and metadata.
#' @noRd
#' 
#' 
# This script imports genotyping data from an Excel file, detects the data type (length-polymorphic or amplicon sequencing)

detect_data_type <- function(input_df) {
  meta_cols_guess <- min(3, ncol(input_df) - 1)
  allele_cols <- (meta_cols_guess + 1):ncol(input_df)
  if (length(allele_cols) == 0) {
    stop("Error: No allele columns found to determine data type.")
  }
  all_allele_values <- unlist(input_df[, allele_cols])
  first_value <- all_allele_values[!is.na(all_allele_values) & all_allele_values != ""][1]
  if (is.na(first_value)) {
    stop("Error: All allele columns are empty. Cannot determine data type.")
  }
  if (!is.na(suppressWarnings(as.numeric(first_value)))) {
    return("length_polymorphic")
  } else {
    return("ampseq")
  }
}


# Import and process genotyping data from an Excel file, detecting the data type and preparing marker information.

import_data <- function(filepath) {
  sheet_names <- try(readxl::excel_sheets(filepath), silent = TRUE)
  if (inherits(sheet_names, "try-error")) {
    stop("ERROR: Cannot read file. Ensure it's a valid Excel file: ", filepath)
  }
  temp_df <- as.data.frame(readxl::read_excel(filepath, sheet = sheet_names[1]))
  data_type <- detect_data_type(temp_df)
  message("INFO: Detected '", data_type, "' data format.")
  
  if (data_type == "length_polymorphic") {
    marker_info <- data.frame(
      marker_id = c("313", "383", "TA1", "POLYA", "PolyA", "PFPK2", "2490", "TA109",
                    "TA42", "TA81", "TA87", "TA40", "ARAII", "PFG377", "TA60",
                    "K1", "3D7", "MAD20", "RO33", "FC27", "IC", "glurp"),
      markertype = c(rep("microsatellite", 15), "msp1", "msp1", "msp1", "msp1", "msp2", "msp2", "glurp"),
      repeatlength = c(2, 2, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, NA, NA, NA, NA, NA, NA, NA), 
      binning_method = c(rep("microsatellite", 15), rep("cluster", 7)),
      cluster_gap_threshold = c(rep(NA, 15), 10, 10, 10, 10, 10, 10, 50)
    )
  } else { # data_type == "ampseq"
    allele_cols <- grep("_allele_\\d+$", colnames(temp_df), value = TRUE)
    base_markers <- unique(gsub("_allele_\\d+$", "", allele_cols))
    
    marker_info <- data.frame(
      marker_id = base_markers,
      markertype = "amplicon",
      repeatlength = NA, 
      binning_method = "exact",
      cluster_gap_threshold = NA,
      stringsAsFactors = FALSE
    )
  }
  
  late_failures_df <- NULL
  additional_df <- NULL
  missing_markers <- c("N/A", "-", "NA", "na", "", " ", "Failed", "failed", "0")
  
  if (data_type == "ampseq") {
    message("INFO: Following amplicon sequencing data processing pipeline.")
    colnames(temp_df)[1:2] <- c("Sample.ID", "Site")
    temp_df[temp_df %in% missing_markers] <- NA
    
    # Clean Sample IDs and check for pairs
    temp_df$Sample.ID <- gsub("D0$", " Day 0", temp_df$Sample.ID)
    temp_df$Sample.ID <- gsub("D[0-9]+$", " Day Failure", temp_df$Sample.ID)
    day0_ids <- unique(gsub(" Day 0", "", temp_df$Sample.ID[grepl(" Day 0", temp_df$Sample.ID)]))
    day_failure_ids <- unique(gsub(" Day Failure", "", temp_df$Sample.ID[grepl(" Day Failure", temp_df$Sample.ID)]))
    
    if (!all(day0_ids %in% day_failure_ids)) {
      missing_pairs <- day0_ids[!day0_ids %in% day_failure_ids]
      stop("ERROR: Ampseq data has 'Day 0' samples missing their 'Day Failure' pair. Missing pairs for: ", 
           paste(missing_pairs, collapse=", "))
    }
    
    late_failures_df <- temp_df
    additional_df <- temp_df[0, ]
    
  } else { # data_type == "length_polymorphic"
    message("INFO: Following length-polymorphic data processing pipeline.")
    
    colnames_cleaned <- tolower(colnames(temp_df)[1:2])
    is_cleaned_format <- all(c("sample.id", "site") %in% colnames_cleaned)
    
    if (is_cleaned_format) {
      late_failures_df <- temp_df
    } else {
      late_failures_df <- data.frame(
        Sample.ID = paste0(temp_df[[1]], "D", temp_df[[2]]),
        Site = temp_df[[3]],
        temp_df[, -c(1, 2, 3), drop = FALSE],
        check.names = FALSE
      )
    }
    
    if (length(sheet_names) > 1) {
      additional_df_raw <- as.data.frame(readxl::read_excel(filepath, sheet = sheet_names[2]))
      if (nrow(additional_df_raw) > 0) {
        add_cleaned <- all(c("sample.id", "site") %in% tolower(colnames(additional_df_raw)[1:2]))
        if(add_cleaned) {
          additional_df <- additional_df_raw
        } else {
          additional_df <- data.frame(
            Sample.ID = paste0(additional_df_raw[[1]], "D", additional_df_raw[[2]]),
            Site = additional_df_raw[[3]],
            additional_df_raw[, -c(1, 2, 3), drop = FALSE],
            check.names = FALSE
          )
        }
      } else {
        additional_df <- late_failures_df[0, ] 
      }
    } else {
      additional_df <- late_failures_df[0, ] 
    }
    
    late_failures_df[late_failures_df %in% missing_markers] <- NA
    additional_df[additional_df %in% missing_markers] <- NA
    
    if (ncol(late_failures_df) > 2) {
      late_failures_df[, 3:ncol(late_failures_df)] <- lapply(late_failures_df[, 3:ncol(late_failures_df)], function(x) as.numeric(as.character(x)))
    }
    if (ncol(additional_df) > 2) {
      additional_df[, 3:ncol(additional_df)] <- lapply(additional_df[, 3:ncol(additional_df)], function(x) as.numeric(as.character(x)))
    }
    
    late_failures_df$Sample.ID <- gsub("D0$", " Day 0", late_failures_df$Sample.ID)
    late_failures_df$Sample.ID <- gsub("D[0-9]+$", " Day Failure", late_failures_df$Sample.ID)

    # Removing uninformative pairs
    all_ids <- unique(gsub(" Day 0| Day Failure", "", late_failures_df$Sample.ID))
    allele_cols <- 3:ncol(late_failures_df) 
    ids_to_remove <- c()
    
    for (id in all_ids) {
      day_failure_row_idx <- which(late_failures_df$Sample.ID == paste(id, "Day Failure"))
      if (length(day_failure_row_idx) > 0) {
        day_failure_alleles <- late_failures_df[day_failure_row_idx, allele_cols]
        if (all(is.na(day_failure_alleles))) {
          ids_to_remove <- c(ids_to_remove, id)
        }
      }
    }
    
    if (length(ids_to_remove) > 0) {
      message("INFO: Removing ", length(ids_to_remove), " patient(s) with no allele data on Day of Failure: ", 
              paste(ids_to_remove, collapse=", "))
      late_failures_df <- late_failures_df[!grepl(paste(ids_to_remove, collapse="|"), late_failures_df$Sample.ID), ]
    }

    day0_ids <- unique(gsub(" Day 0", "", late_failures_df$Sample.ID[grepl(" Day 0", late_failures_df$Sample.ID)]))
    day_failure_ids <- unique(gsub(" Day Failure", "", late_failures_df$Sample.ID[grepl(" Day Failure", late_failures_df$Sample.ID)]))
    
    if (!all(day0_ids %in% day_failure_ids)) {
      missing_pairs <- day0_ids[!day0_ids %in% day_failure_ids]
      stop("ERROR: Length data has 'Day 0' samples missing their 'Day Failure' pair. Missing pairs for: ", 
           paste(missing_pairs, collapse=", "))
    }
  }
  return(
    list(
      late_failures = late_failures_df,
      additional = additional_df,
      marker_info = marker_info,
      data_type = data_type
    )
  )
}