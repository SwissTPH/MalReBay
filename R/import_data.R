#' Import Genotyping Data from Excel
#'
#' @description Reads data from a specified Excel file, automatically detecting
#'   the data type and separating sheets into a structured list.
#'
#' @param filepath The full path to the input Excel file.
#' @param verbose Logical. If TRUE, prints progress and data-cleaning messages.
#' @param marker_filepath Path to Excel file containing marker metadata 
#' @param additional_filepath Optional path to a separate additional/background
#'   data file (csv or xlsx). Ignored if the main file already has 2 sheets.
#' (optional if marker_info sheet is present)
#' @return A list containing the imported data.
#'
#' @examples
#' \dontrun{
#'   data_file <- system.file("extdata", 
#'                            "Angola_2021_TES_7NMS.xlsx", 
#'                            package = "MalReBay")
#'   marker_file <- system.file("extdata", 
#'                              "makers_details.xlsx", 
#'                              package = "MalReBay")
#'   imported    <- import_data(filepath = data_file, 
#'                              marker_filepath = marker_file)
#' }
#' @export
import_data <- function(
    filepath = system.file("extdata", 
                           "Angola_2021_TES_7NMS.xlsx", 
                           package = "MalReBay"),
    marker_filepath = system.file("extdata", 
                                  "makers_details.xlsx", 
                                  package = "MalReBay"),
    additional_filepath = NULL,
    verbose = TRUE) {
  
  # Load Workbook and Sheets
  ext <- tools::file_ext(filepath)
  
  if (ext == "csv") {
    temp_df <- as.data.frame(readr::read_csv(filepath, show_col_types = FALSE))
    sheet_names <- NULL
  } else {
    sheet_names <- try(readxl::excel_sheets(filepath), silent = TRUE)
    if (inherits(sheet_names, "try-error")) stop("ERROR: Cannot read file: ", filepath)
    temp_df <- as.data.frame(readxl::read_excel(filepath, sheet = sheet_names[1]))
  }
  
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
  
  # Load Marker Metadata
  if (!is.null(marker_filepath) && file.exists(marker_filepath)) {
    marker_info <- as.data.frame(readxl::read_excel(marker_filepath))
  } else if (!is.null(sheet_names) && "marker_info" %in% sheet_names) {
    marker_info <- as.data.frame(readxl::read_excel(filepath, sheet = "marker_info"))
  } else {
    stop("ERROR: Marker information not found. Provide marker_filepath or add a 'marker_info' sheet to your Excel file.")
  }
  
  marker_info$marker_id <- as.character(marker_info$marker_id)
  marker_info$repeatlength <- suppressWarnings(as.numeric(as.character(marker_info$repeatlength)))
  

  # Process additional/background data. Two possible sources:
  #  1. An embedded second sheet in the main xlsx file (existing behaviour).
  #  2. A separately supplied additional_filepath (csv or xlsx).
  # If both are present, the embedded sheet wins -- merging two independently
  # prepared background sources automatically risks double-counting overlapping
  # samples, which would bias freq[] estimation, so we pick one deterministically
  # and tell the user which, rather than guessing at a merge.
  additional_df <- late_failures_df[0, ] # Default empty
  
  if (!is.null(sheet_names) && length(sheet_names) > 1) {
    raw_add <- as.data.frame(readxl::read_excel(filepath, sheet = sheet_names[2]))
    if (nrow(raw_add) > 0) {
      # Apply same 3-column metadata logic to sheet 2
      if (all(c("sample.id", "site") %in% tolower(colnames(raw_add)[1:2]))) {
        
        colnames(raw_add)[1:2] <- c("Sample.ID", "Site")
        raw_add$Sample.ID <- gsub("_?D0$| D0$| Day 0$", " Day 0", raw_add$Sample.ID)
        is_day0 <- grepl(" Day 0$", raw_add$Sample.ID)
        raw_add$Sample.ID[!is_day0] <- gsub("(_D[0-9A-Za-z]+|D[0-9]+)$| recurrence$| Day Failure$", " recurrence", raw_add$Sample.ID[!is_day0])
        unlabelled <- !grepl(" Day 0$| recurrence$", raw_add$Sample.ID)
        raw_add$Sample.ID[unlabelled] <- paste(raw_add$Sample.ID[unlabelled], "Day 0")
        
        additional_df <- clean_data(raw_add)
      } else {
        additional_clean <- data.frame(
          Sample.ID = paste0(raw_add[[1]], " Day ", ifelse(raw_add[[2]] == 0, "0", "recurrence")),
          Site = raw_add[[3]], raw_add[, -c(1, 2, 3), drop = FALSE], check.names = FALSE
        )
        additional_df <- clean_data(additional_clean)
      }
    }
    
    if (verbose) {
      message("INFO: Detected 2 sheets in ", basename(filepath), ": ",
              "Sheet 1 (late failures) has ", nrow(late_failures_df), " sample(s), ",
              "Sheet 2 (additional) has ", nrow(additional_df), " sample(s).")
      if (!is.null(additional_filepath)) {
        message("INFO: Ignoring additional_filepath since Sheet 2 already ",
                "supplies additional data -- combine your data into Sheet 2 ",
                "if you need both sources included.")
      }
    }
    
  } else if (!is.null(additional_filepath)) {
    if (!file.exists(additional_filepath)) {
      stop("ERROR: additional_filepath does not exist: ", additional_filepath)
    }
    if (normalizePath(additional_filepath) == normalizePath(filepath)) {
      stop("ERROR: additional_filepath is the same file as filepath -- ",
           "these must be different files.")
    }
    
    add_ext <- tools::file_ext(additional_filepath)
    raw_add <- if (add_ext == "csv") {
      as.data.frame(readr::read_csv(additional_filepath, show_col_types = FALSE))
    } else {
      add_sheet_names <- try(readxl::excel_sheets(additional_filepath), silent = TRUE)
      if (inherits(add_sheet_names, "try-error")) stop("ERROR: Cannot read file: ", additional_filepath)
      if (length(add_sheet_names) > 1 && verbose) {
        message("INFO: additional_filepath has ", length(add_sheet_names),
                " sheet(s) -- only using the first sheet ('", add_sheet_names[1], "').")
      }
      as.data.frame(readxl::read_excel(additional_filepath, sheet = add_sheet_names[1]))
    }
    
    if (nrow(raw_add) > 0) {
      if (all(c("sample.id", "site") %in% tolower(colnames(raw_add)[1:2]))) {
        
        colnames(raw_add)[1:2] <- c("Sample.ID", "Site")
        raw_add$Sample.ID <- gsub("_?D0$| D0$| Day 0$", " Day 0", raw_add$Sample.ID)
        is_day0 <- grepl(" Day 0$", raw_add$Sample.ID)
        raw_add$Sample.ID[!is_day0] <- gsub("(_D[0-9A-Za-z]+|D[0-9]+)$| recurrence$| Day Failure$", " recurrence", raw_add$Sample.ID[!is_day0])
        unlabelled <- !grepl(" Day 0$| recurrence$", raw_add$Sample.ID)
        raw_add$Sample.ID[unlabelled] <- paste(raw_add$Sample.ID[unlabelled], "Day 0")
        
        additional_df <- clean_data(raw_add)
      } else {
        additional_clean <- data.frame(
          Sample.ID = paste0(raw_add[[1]], " Day ", ifelse(raw_add[[2]] == 0, "0", "recurrence")),
          Site = raw_add[[3]], raw_add[, -c(1, 2, 3), drop = FALSE], check.names = FALSE
        )
        additional_df <- clean_data(additional_clean)
      }
      
      if (!all(grepl(" Day 0$", additional_df$Sample.ID))) {
        stop("ERROR: additional_filepath must contain only Day 0 samples -- ",
             "found recurrence/non-Day-0 records. Check that you're pointing ",
             "at background/additional data, not a late-failures file.")
      }
    }
    
    if (verbose) {
      message("INFO: Loaded additional data from ", basename(additional_filepath),
              ": ", nrow(additional_df), " sample(s).")
    }
  }
  
  # Site is an identifier, not a quantity -- some source files store it as a
  # bare number (e.g. a numeric site code) rather than a name. Downstream
  # code (e.g. posterior_probabilities$Site, built from site names used as
  # list keys, which R always coerces to character) always treats Site as
  # character, so keep it consistent here too or joins on Site fail with a
  # type mismatch.
  late_failures_df$Site <- as.character(late_failures_df$Site)
  additional_df$Site    <- as.character(additional_df$Site)

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
  detected_markers <- unique(base_names_in_data)

  markers_to_use <- intersect(marker_info$marker_id, detected_markers)
  if (length(markers_to_use) == 0) stop("No matching markers found between metadata and data.")

  unmatched_markers <- setdiff(detected_markers, markers_to_use)

  if (verbose) {
    message("")
    message("INFO: Detected ", length(detected_markers), " marker(s) in the data file.")
    if (length(unmatched_markers) > 0) {
      message("INFO: ", length(unmatched_markers),
              " marker(s) will NOT be processed (no matching entry in the marker metadata):")
      message("      ", paste(unmatched_markers, collapse = ", "))
    }
    message("INFO: Proceeding with analysis using ", length(markers_to_use), " marker(s):")
    message("      ", paste(markers_to_use, collapse = ", "))
    message("")
  }

  # Subset dataframes to only include columns belonging to valid markers
  valid_cols <- allele_colnames[base_names_in_data %in% markers_to_use]
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
  
  # Subset main data and metadata to final selection
  late_failures_df <- late_failures_df[, c(colnames(late_failures_df)[1:2], valid_cols)]
  marker_info <- marker_info[marker_info$marker_id %in% markers_to_use, ]
  
  return(list(
    late_failures = late_failures_df,
    additional = additional_df,
    marker_info = marker_info,
    data_type = data_type
  ))
}