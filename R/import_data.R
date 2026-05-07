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
import_data <- function(filepath        = system.file("extdata", "Angola_2021_TES_7NMS.xlsx", package = "MalReBay"),
                        marker_filepath = system.file("extdata", "makers_details.xlsx",        package = "MalReBay"),
                        verbose         = TRUE) {
  # Load Workbook and Sheets
  sheet_names <- try(readxl::excel_sheets(filepath), silent = TRUE)
  if (inherits(sheet_names, "try-error")) stop("ERROR: Cannot read file: ", filepath)
  
  # Load Primary Data (Sheet 1)
  temp_df <- as.data.frame(readxl::read_excel(filepath, sheet = sheet_names[1]))
  
  # Need ID, Site/Day, and at least 2 allele columns to be valid for processing
  if (ncol(temp_df) < 4) {
    stop("ERROR: Input data must have at least 4 columns (e.g., ID, Site, and at least 2 allele columns).")
  }

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
  marker_info$repeatlength <- suppressWarnings(as.numeric(as.character(marker_info$repeatlength)))

  # Standardize metadata for the columns (1, 2, and 3)
  # Column 1 = Sample.ID, Column 2 = Site, Columns 3+ = Alleles
  cols_lower <- tolower(colnames(temp_df))
  
  if ("day" %in% cols_lower) {
    # If "Day" column exists, merge ID and Day, then use Site as second column
    id_idx   <- which(cols_lower %in% c("sample.id", "sample_id", "sampleid",
                                        "id", "patientid", "patient_id"))[1]
    day_idx  <- which(cols_lower == "day")[1]
    site_idx <- which(cols_lower == "site")[1]
    
    # Validate all required columns were found before subsetting
    if (is.na(id_idx))   stop("ERROR: Cannot find ID column.")
    if (is.na(site_idx)) stop("ERROR: Cannot find Site column.")
    
    # Create the standardized dataframe
    late_failures_df <- data.frame(
      Sample.ID = paste0(temp_df[[id_idx]], ifelse(temp_df[[day_idx]] == 0, " Day 0", " recurrence")),
      Site      = temp_df[[site_idx]],
      temp_df[, -c(id_idx, day_idx, site_idx), drop = FALSE],
      check.names = FALSE
    )
  } else {
    # No Day column: Assume Col 1 is ID, Col 2 is Site (standard format)
    colnames(temp_df)[1:2] <- c("Sample.ID", "Site")
    
    # Clean ID names to ensure " Day 0" or " recurrence" suffix
    temp_df$Sample.ID <- gsub("_?D0$| D0$| Day 0$", " Day 0", temp_df$Sample.ID)
    is_day0 <- grepl(" Day 0$", temp_df$Sample.ID)
    temp_df$Sample.ID[!is_day0] <- gsub("(_D[0-9A-Za-z]+|D[0-9]+)$| recurrence$| Day Failure$", " recurrence", temp_df$Sample.ID[!is_day0])
    
    # Unlabelled check
    unlabelled <- !grepl(" Day 0$| recurrence$", temp_df$Sample.ID)
    temp_df$Sample.ID[unlabelled] <- paste(temp_df$Sample.ID[unlabelled], "Day 0")
    
    late_failures_df <- temp_df
  }
  
  # Cleaning missing values
  missing_alleles <- c("N/A", "-", "NA", "na", "", " ", "Failed", "failed", "0")
  
  clean_data <- function(df) {
    if (nrow(df) == 0) return(df)
    df[df %in% missing_alleles] <- NA
    # Convert allele columns (3 onwards) to numeric if length polymorphic
    if (data_type == "length_polymorphic" && ncol(df) > 2) {
      df[, 3:ncol(df)] <- lapply(df[, 3:ncol(df)], function(x) suppressWarnings(as.numeric(as.character(x))))
    }
    return(df)
  }
  
  late_failures_df <- clean_data(late_failures_df)

  # Process the additional sheet (Sheet 2)
  additional_df <- late_failures_df[0, ] # Default empty
  if (data_type == "length_polymorphic" && length(sheet_names) > 1) {
    raw_add <- as.data.frame(readxl::read_excel(filepath, sheet = sheet_names[2]))
    if (nrow(raw_add) > 0) {
      # Apply same 3-column metadata logic to sheet 2
      if (all(c("sample.id", "site") %in% tolower(colnames(raw_add)[1:2]))) {
        additional_df <- raw_add
      } else {
        additional_clean <- data.frame(
          Sample.ID = paste0(raw_add[[1]], " Day ", ifelse(raw_add[[2]] == 0, "0", "recurrence")),
          Site = raw_add[[3]], raw_add[, -c(1, 2, 3), drop = FALSE], check.names = FALSE
        )
        additional_df <- clean_data(additional_clean)
      }
    }
  }

  # Remove "Failure" rows that are entirely empty (No allele data)
  allele_idx <- 3:ncol(late_failures_df)
  
  if (length(allele_idx) > 0) {
    is_failure <- grepl("recurrence", late_failures_df$Sample.ID)
    empty_alleles <- rowSums(!is.na(late_failures_df[, allele_idx, drop = FALSE])) == 0
    to_remove_ids <- unique(gsub(" recurrence", "", late_failures_df$Sample.ID[is_failure & empty_alleles]))
    
    if (length(to_remove_ids) > 0) {
      if (verbose) message("INFO: Removing ", length(to_remove_ids), " patient(s) with no data on recurrence.")
      pattern <- paste(to_remove_ids, collapse = "|")
      late_failures_df <- late_failures_df[!grepl(pattern, late_failures_df$Sample.ID), ]
    }
  }

  # Check if all Day 0 samples have corresponding recurrence samples
  day0_ids <- gsub(" Day 0", "", late_failures_df$Sample.ID[grepl(" Day 0", late_failures_df$Sample.ID)])
  fail_ids <- gsub(" recurrence", "", late_failures_df$Sample.ID[grepl(" recurrence", late_failures_df$Sample.ID)])
  
  missing_pairs <- setdiff(day0_ids, fail_ids)
  if (length(missing_pairs) > 0 && verbose) {
    warning("WARNING: 'Day 0' samples missing their 'recurrence' pair: ", paste(missing_pairs, collapse=", "))
  }
  
  # Filter to markers that are actually present in the data
  marker_suffix_regex <- "(_allele_|_)\\d+$"
  allele_colnames <- colnames(late_failures_df)[3:ncol(late_failures_df)]
  base_names_in_data <- gsub(marker_suffix_regex, "", allele_colnames)
  
  markers_to_use <- intersect(marker_info$marker_id, unique(base_names_in_data))
  if (length(markers_to_use) == 0) stop("No matching markers found between metadata and data.")
  
  # Subset dataframes to only include columns belonging to valid markers
  valid_cols <- allele_colnames[base_names_in_data %in% markers_to_use]
  
  if (nrow(additional_df) > 0) {
    missing_cols <- setdiff(valid_cols, colnames(additional_df))
    if (length(missing_cols) > 0) {
      stop("ERROR: Additional data sheet is missing marker columns found in the main sheet: ", 
           paste(missing_cols, collapse = ", "))
    }
    
    extra_cols <- setdiff(colnames(additional_df)[3:ncol(additional_df)], valid_cols)
    if (length(extra_cols) > 0 && verbose) {
      message("INFO: Dropping ", length(extra_cols), " extra column(s) from additional data: ", 
              paste(extra_cols, collapse = ", "))
    }
    # Sync additional_df columns
    additional_df <- additional_df[, c(colnames(additional_df)[1:2], valid_cols)]
  }
  
  # Subset main data and metadata to final selection
  late_failures_df <- late_failures_df[, c(colnames(late_failures_df)[1:2], valid_cols)]
  marker_info <- marker_info[marker_info$marker_id %in% markers_to_use, ]
  
  if (verbose) {
    message("INFO: Using ", length(markers_to_use), " markers: ", paste(markers_to_use, collapse = ", "))
  }

  return(list(
    late_failures = late_failures_df,
    additional = additional_df,
    marker_info = marker_info,
    data_type = data_type
  ))
}