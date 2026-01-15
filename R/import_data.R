#' Detect Allele Data Type
#'
#' @description
#' Examines allele data to determine if it is length-polymorphic (numeric
#' fragments) or amplicon sequencing (character alleles). It intelligently
#' skips metadata columns to inspect only the allele values.
#'
#' @param input_df A data frame containing genotyping data, as read from the
#'   input file.
#'
#' @return A character string: either `"length_polymorphic"` or `"ampseq"`.
#'
#' @keywords internal
#' @noRd
#'
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

#' Import Genotyping Data from Excel
#'
#' @description Reads data from a specified Excel file, automatically detecting
#'   the data type and separating sheets into a structured list.
#'
#' @param filepath The full path to the input Excel file.
#' @param verbose Logical. If TRUE, prints progress and data-cleaning messages.
#' @return A list containing the imported data.
#' @export
import_data <- function(filepath, verbose = TRUE) {
  sheet_names <- try(readxl::excel_sheets(filepath), silent = TRUE)
  if (inherits(sheet_names, "try-error")) {
    stop("ERROR: Cannot read file. Ensure it's a valid Excel file: ", filepath)
  }
  
  temp_df <- as.data.frame(readxl::read_excel(filepath, sheet = sheet_names[1]))
  data_type <- detect_data_type(temp_df) 
  
  if (verbose) message("INFO: Detected '", data_type, "' data format.")
  
  if (data_type == "length_polymorphic") {
    marker_info <- data.frame(
      marker_id = c("313", "383", "TA1", "POLYA", "PolyA", "PFPK2", "2490", "TA109",
                    "TA42", "TA81", "TA87", "TA40", "ARAII", "PFG377", "TA60",
                    "m1", "m2", "m3", "m4", "m5", 
                    "K1", "MAD20", "RO33", "R033", "FC27", "IC", "3D7", "glurp"),
      markertype = c(rep("microsatellite", 20), "msp1", "msp1", "msp1", "msp1", "msp2", "msp2", "msp2", "glurp"),
      repeatlength = c(rep(2, 2), 3, 3, 3, 3, 3, 3, rep(2, 7), rep(3, 5), rep(NA, 8)),
      binning_method = c(rep("microsatellite", 15), rep("microsatellite", 5), rep("cluster", 8)),
      cluster_gap_threshold = c(rep(NA, 20), rep(2, 7), 2)  
    )
  } else if (data_type == "ampseq") {
    marker_info <- data.frame(
      marker_id = c("cpmp", "cpp", "amaD3", "ama1-D3", "csp", "msp7", "SERA2", "TRAP3", "cpp2", "csp2"),
      markertype = "amplicon",
      repeatlength = NA,
      binning_method = "exact",
      cluster_gap_threshold = NA
    )
  }
  
  late_failures_df <- NULL
  additional_df <- NULL
  missing_markers <- c("N/A", "-", "NA", "na", "", " ", "Failed", "failed", "0")
  
  if (data_type == "ampseq") {
    if (verbose) message("INFO: Following amplicon sequencing data processing pipeline.")
    cols_lower <- tolower(colnames(temp_df))
    
    if ("day" %in% cols_lower) {
      
      id_idx <- which(cols_lower %in% c("sample.id", "sample_id", "sampleid", "id"))[1]
      day_idx <- which(cols_lower == "day")[1]
      site_idx <- which(cols_lower == "site")[1]
      
      if (is.na(id_idx) || is.na(day_idx) || is.na(site_idx)) {
        stop("ERROR: 'Day' column detected but could not locate 'Sample.ID' or 'Site' columns.")
      }
      
      colnames(temp_df)[c(id_idx, day_idx, site_idx)] <- c("Sample.ID", "Day", "Site")
      temp_df$Sample.ID <- ifelse(temp_df$Day == 0, 
                                  paste0(temp_df$Sample.ID, " Day 0"), 
                                  paste0(temp_df$Sample.ID, " Day Failure"))
      
      allele_cols <- setdiff(names(temp_df), c("Sample.ID", "Day", "Site"))
      temp_df <- temp_df[, c("Sample.ID", "Site", allele_cols)]
      
    } else {
      colnames(temp_df)[1:2] <- c("Sample.ID", "Site")
      temp_df$Sample.ID <- gsub("_?D0$", " Day 0", temp_df$Sample.ID)    
      is_day0 <- grepl(" Day 0$", temp_df$Sample.ID)
      temp_df$Sample.ID[!is_day0] <- gsub("_?D([^0].*|0.+)$", " Day Failure", temp_df$Sample.ID[!is_day0])
    }
    
    temp_df[temp_df %in% missing_markers] <- NA
    
    day0_ids <- unique(gsub(" Day 0", "", temp_df$Sample.ID[grepl(" Day 0", temp_df$Sample.ID)]))
    day_failure_ids <- unique(gsub(" Day Failure", "", temp_df$Sample.ID[grepl(" Day Failure", temp_df$Sample.ID)]))
    
    if (verbose) {
      message("DEBUG: Found ", length(day0_ids), " Day 0 IDs and ", length(day_failure_ids), " Failure IDs.")
    }
    
    if (!all(day0_ids %in% day_failure_ids)) {
      missing_pairs <- day0_ids[!day0_ids %in% day_failure_ids]
      warning("WARNING: Ampseq data has 'Day 0' samples missing their 'Day Failure' pair: ",
              paste(missing_pairs, collapse=", "))
    }
    late_failures_df <- temp_df
    additional_df <- temp_df[0, ]
    
  } else { 
    
    if (verbose) message("INFO: Following length-polymorphic data processing pipeline.")
    
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
      } else { additional_df <- late_failures_df[0, ] }
    } else { additional_df <- late_failures_df[0, ] }
    
    late_failures_df[late_failures_df %in% missing_markers] <- NA
    additional_df[additional_df %in% missing_markers] <- NA
    
    if (ncol(late_failures_df) > 2) {
      late_failures_df[, 3:ncol(late_failures_df)] <- suppressWarnings(
        lapply(late_failures_df[, 3:ncol(late_failures_df)], function(x) as.numeric(as.character(x)))
      )
    }
    if (ncol(additional_df) > 2) {
      additional_df[, 3:ncol(additional_df)] <- suppressWarnings(
        lapply(additional_df[, 3:ncol(additional_df)], function(x) as.numeric(as.character(x)))
      )
    }
    
    late_failures_df$Sample.ID <- gsub("D0$", " Day 0", late_failures_df$Sample.ID)
    late_failures_df$Sample.ID <- gsub("D[0-9]+$", " Day Failure", late_failures_df$Sample.ID)
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
      if (verbose) {
        message("INFO: Removing ", length(ids_to_remove), " patient(s) with no allele data on Day of Failure.")
      }
      late_failures_df <- late_failures_df[!grepl(paste(ids_to_remove, collapse="|"), late_failures_df$Sample.ID), ]
    }
    
    day0_ids <- unique(gsub(" Day 0", "", late_failures_df$Sample.ID[grepl(" Day 0", late_failures_df$Sample.ID)]))
    day_failure_ids <- unique(gsub(" Day Failure", "", late_failures_df$Sample.ID[grepl(" Day Failure", late_failures_df$Sample.ID)]))
    
    if (!all(day0_ids %in% day_failure_ids)) {
      missing_pairs <- day0_ids[!day0_ids %in% day_failure_ids]
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