#' Import Genotyping Data from Excel
#'
#' @description Reads data from a specified Excel file, automatically detecting
#'   the data type and separating sheets into a structured list.
#'
#' @param filepath The full path to the input Excel file.
#' @param verbose Logical. If TRUE, prints progress and data-cleaning messages.
#' @param marker_filepath Path to Excel file containing marker metadata (optional if marker_info sheet is present)
#' @return A list containing the imported data.
#' @export
import_data <- function(filepath, marker_filepath = NULL, verbose = TRUE) {
  # Load Workbook and Sheets
  sheet_names <- try(readxl::excel_sheets(filepath), silent = TRUE)
  if (inherits(sheet_names, "try-error")) stop("ERROR: Cannot read file: ", filepath)
  
  # Load Primary Data (Sheet 1)
  temp_df <- as.data.frame(readxl::read_excel(filepath, sheet = sheet_names[1]))
  if (ncol(temp_df) < 3) stop("ERROR: Input data must have at least 3 columns (Metadata + Alleles).")

  # Detect data type based on first allele value
  # Look at columns 4 onwards to find the first actual allele value
  sample_values <- unlist(temp_df[, 4:ncol(temp_df)])
  first_val <- sample_values[!is.na(sample_values) & sample_values != ""][1]
  
  if (is.na(first_val)) stop("ERROR: All allele columns are empty. cannot detect data type.")
  
  # If the first value is numeric, it's length_polymorphic; otherwise, it's ampseq
  data_type <- if (!is.na(suppressWarnings(as.numeric(first_val)))) "length_polymorphic" else "ampseq"
  if (verbose) message("INFO: Detected '", data_type, "' data format.")

  # Load Marker Metadata
  if (!is.null(marker_filepath) && file.exists(marker_filepath)) {
    marker_info <- as.data.frame(readxl::read_excel(marker_filepath))
  } else if ("marker_info" %in% sheet_names) {
    marker_info <- as.data.frame(readxl::read_excel(filepath, sheet = "marker_info"))
  } else {
    stop("ERROR: Marker information not found. Provide marker_filepath or add a 'marker_info' sheet to your Excel file.")
  }

  # Ensure marker_id is always a string
  marker_info$marker_id <- as.character(marker_info$marker_id)
  marker_info$repeatlength <- as.numeric(as.character(marker_info$repeatlength))

  # Standardize metadata for the columns (1, 2, and 3)
  # Column 1 = Sample.ID, Column 2 = Site, Columns 3+ = Alleles
  cols_lower <- tolower(colnames(temp_df))
  
  if ("day" %in% cols_lower) {
    # If "Day" column exists, merge ID and Day, then use Site as 2nd col
    id_idx <- which(cols_lower %in% c("sample.id", "sample_id", "sampleid", "id"))[1]
    day_idx <- which(cols_lower == "day")[1]
    site_idx <- which(cols_lower == "site")[1]
    
    # Create the standardized dataframe
    late_failures_df <- data.frame(
      Sample.ID = paste0(temp_df[[id_idx]], " Day ", ifelse(temp_df[[day_idx]] == 0, "0", "Failure")),
      Site = temp_df[[site_idx]],
      temp_df[, -c(id_idx, day_idx, site_idx), drop = FALSE],
      check.names = FALSE
    )
  } else {
    # No Day column: Assume Col 1 is ID, Col 2 is Site (standard format)
    colnames(temp_df)[1:2] <- c("Sample.ID", "Site")
    
    # Clean ID names to ensure "Day 0" or "Day Failure" suffix
    temp_df$Sample.ID <- gsub("_?D0$| Day 0$", " Day 0", temp_df$Sample.ID)
    is_day0 <- grepl(" Day 0$", temp_df$Sample.ID)
    temp_df$Sample.ID[!is_day0] <- gsub("_?D[1-9][0-9]*$| Day Failure$", " Day Failure", temp_df$Sample.ID[!is_day0])
    
    # If it wasn't already suffixed, force suffix based on common logic
    unlabelled <- !grepl(" Day (0|Failure)$", temp_df$Sample.ID)
    temp_df$Sample.ID[unlabelled] <- paste(temp_df$Sample.ID[unlabelled], "Day 0")
    
    late_failures_df <- temp_df
  }

  # Process the additional sheet (Sheet 2)
  additional_df <- late_failures_df[0, ] # Default empty
  if (data_type == "length_polymorphic" && length(sheet_names) > 1) {
    raw_add <- as.data.frame(readxl::read_excel(filepath, sheet = sheet_names[2]))
    if (nrow(raw_add) > 0) {
      # Apply same 3-column metadata logic to sheet 2
      if (all(c("sample.id", "site") %in% tolower(colnames(raw_add)[1:2]))) {
        additional_df <- raw_add
      } else {
        additional_df <- data.frame(
          Sample.ID = paste0(raw_add[[1]], " Day ", ifelse(raw_add[[2]] == 0, "0", "Failure")),
          Site = raw_add[[3]], raw_add[, -c(1, 2, 3), drop = FALSE], check.names = FALSE
        )
      }
    }
  }

  # Data cleaning
  missing_alleles <- c("N/A", "-", "NA", "na", "", " ", "Failed", "failed", "0")
  
  clean_data <- function(df) {
    if (nrow(df) == 0) return(df)
    df[df %in% missing_alleles] <- NA
    # Convert allele columns (3 onwards) to numeric if length polymorphic
    if (data_type == "length_polymorphic" && ncol(df) > 2) {
      df[, 3:ncol(df)] <- lapply(df[, 3:ncol(df)], function(x) as.numeric(as.character(x)))
    }
    return(df)
  }

  late_failures_df <- clean_data(late_failures_df)
  additional_df <- clean_data(additional_df)

  # 6. Remove "Failure" rows that are entirely empty (No allele data)
  allele_cols <- 3:ncol(late_failures_df)
  if (length(allele_cols) > 0) {
    is_failure <- grepl("Day Failure", late_failures_df$Sample.ID)
    empty_alleles <- rowSums(!is.na(late_failures_df[, allele_cols, drop = FALSE])) == 0
    to_remove_ids <- unique(gsub(" Day Failure", "", late_failures_df$Sample.ID[is_failure & empty_alleles]))
    
    if (length(to_remove_ids) > 0) {
      if (verbose) message("INFO: Removing ", length(to_remove_ids), " patient(s) with no data on Day Failure.")
      pattern <- paste(to_remove_ids, collapse = "|")
      late_failures_df <- late_failures_df[!grepl(pattern, late_failures_df$Sample.ID), ]
    }
  }

  # Check if all Day 0 samples have corresponding Day Failure samples
  day0_ids <- gsub(" Day 0", "", late_failures_df$Sample.ID[grepl(" Day 0", late_failures_df$Sample.ID)])
  fail_ids <- gsub(" Day Failure", "", late_failures_df$Sample.ID[grepl(" Day Failure", late_failures_df$Sample.ID)])
  
  missing_pairs <- setdiff(day0_ids, fail_ids)
  if (length(missing_pairs) > 0 && verbose) {
    warning("WARNING: 'Day 0' samples missing their 'Day Failure' pair: ", paste(missing_pairs, collapse=", "))
  }

  return(list(
    late_failures = late_failures_df,
    additional = additional_df,
    marker_info = marker_info,
    data_type = data_type
  ))
}