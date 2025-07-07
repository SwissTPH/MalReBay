# This file defines a function to import and clean the genotyping data.
# It takes a file path as input and returns a list with two data frames:
# 'late_failures' and 'additional'.

import_and_clean_data <- function(filepath) {
  sheet_names <- try(readxl::excel_sheets(filepath), silent = TRUE)
  if (inherits(sheet_names, "try-error") || length(sheet_names) == 0) {
    stop("ERROR: Cannot read sheets from the file. Please ensure it is a valid Excel file: ", filepath)
  }
  
  # Read Late Treatment ailures Data as first sheet
  late_failures_sheet <- sheet_names[1]
  late_failures_df <- as.data.frame(readxl::read_excel(filepath, sheet = late_failures_sheet, skip = 3))
  
  # Robustly Rename ID Column (First Column)
  if (ncol(late_failures_df) < 1) {
    stop("The first sheet '", late_failures_sheet, "' appears to be empty or has no columns.")
  }
  original_id_col_name <- names(late_failures_df)[1]
  names(late_failures_df)[1] <- "Sample.ID"
  message("INFO: Renamed first column from '", original_id_col_name, "' to 'Sample.ID'")
  
  # Read Additional Data from the second Sheet, if it exists
  if (length(sheet_names) > 1) {
    additional_sheet <- sheet_names[2]
    additional_df <- as.data.frame(readxl::read_excel(filepath, sheet = additional_sheet, skip = 3))
    
    if (ncol(additional_df) > 0) {
      original_add_id_col <- names(additional_df)[1]
      names(additional_df)[1] <- "Sample.ID"
      message("INFO: Renamed first column of additional data from '", original_add_id_col, "' to 'Sample.ID'")
    }
    
  } else {
    message("INFO: No second sheet found for additional data. Creating an empty placeholder.")
    additional_df <- late_failures_df[0, , drop = FALSE]
  }
  missing_markers <- c(0, "0", "N/A", "-", "NA", "na", "", " ", "Failed", "failed")
  
  # Clean late failures data
  late_failures_df[late_failures_df %in% missing_markers] <- NA
  late_failures_df$Sample.ID <- sub("D0$", " Day 0", late_failures_df$Sample.ID)
  late_failures_df$Sample.ID <- sub("D[0-9]+$", " Day Failure", late_failures_df$Sample.ID)
  
  # Validate pairing
  day0_ids <- unique(unlist(strsplit(late_failures_df$Sample.ID[grepl("Day 0", late_failures_df$Sample.ID)], " Day 0")))
  dayF_ids <- unique(unlist(strsplit(late_failures_df$Sample.ID[grepl("Day Failure", late_failures_df$Sample.ID)], " Day Failure")))
  
  if (any(!paste(day0_ids, "Day Failure") %in% late_failures_df$Sample.ID)) {
    stop("ERROR: Some 'Day 0' samples are missing their 'Day Failure' pair. Please check your data.")
  }
  if (any(!paste(dayF_ids, "Day 0") %in% late_failures_df$Sample.ID)) {
    stop("ERROR: Some 'Day Failure' samples are missing their 'Day 0' pair. Please check your data.")
  }
  
  # Clean additional data (if it has rows)
  if (nrow(additional_df) > 0) {
    additional_df[additional_df %in% missing_markers] <- NA
    additional_df$Sample.ID <- sub("_D0", " Day 0", additional_df$Sample.ID)
    additional_df$Sample.ID <- sub("_D[0-9]*", " Day Failure", additional_df$Sample.ID)
  }
  
  # Convert marker columns to numeric; assuming first two columns are non-numeric (e.g., Sample.ID, Site)
  if (ncol(late_failures_df) > 2) {
    cols_to_convert_late <- colnames(late_failures_df)[-c(1, 2)]
    late_failures_df[, cols_to_convert_late] <- suppressWarnings(lapply(late_failures_df[, cols_to_convert_late], as.numeric))
  }
  if (nrow(additional_df) > 0 && ncol(additional_df) > 2) {
    cols_to_convert_add <- colnames(additional_df)[-c(1, 2)]
    additional_df[, cols_to_convert_add] <- suppressWarnings(lapply(additional_df[, cols_to_convert_add], as.numeric))
  }
  
  return(
    list(
      late_failures = late_failures_df,
      additional = additional_df
    )
  )
}