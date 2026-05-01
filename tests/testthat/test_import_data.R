library(testthat)
library(MalReBay)
#   angola_imported_data, create_mock_microsat_data(),
#   create_mock_markers(), create_mock_xlsx(), mock_allele_cols()

# ============================================================
# Happy path — real Angola dataset (pre-imported from .rds)
# These tests inspect the known-good imported object directly,
# avoiding system.file() / path resolution issues during testing.
# ============================================================

test_that("import_data: returns correct list structure", {
  expect_type(angola_imported_data, "list")
  expect_named(angola_imported_data,
               c("late_failures", "additional", "marker_info", "data_type"))
})

test_that("import_data: detects length_polymorphic data type correctly", {
  expect_equal(angola_imported_data$data_type, "length_polymorphic")
})

test_that("import_data: late_failures has required columns", {
  expect_true("Sample.ID" %in% colnames(angola_imported_data$late_failures))
  expect_true("Site"      %in% colnames(angola_imported_data$late_failures))
  expect_gte(ncol(angola_imported_data$late_failures), 3)
})

test_that("import_data: Sample.ID contains only 'Day 0' or 'recurrence' suffixes", {
  ids   <- angola_imported_data$late_failures$Sample.ID
  valid <- grepl(" Day 0$| recurrence$", ids)
  expect_true(all(valid),
              info = paste("Invalid IDs:", paste(ids[!valid], collapse = ", ")))
})

test_that("import_data: every recurrence row has a matching Day 0 row", {
  ids       <- angola_imported_data$late_failures$Sample.ID
  day0_ids  <- gsub(" Day 0$",      "", ids[grepl(" Day 0$",      ids)])
  recur_ids <- gsub(" recurrence$", "", ids[grepl(" recurrence$", ids)])
  missing   <- setdiff(recur_ids, day0_ids)
  expect_equal(length(missing), 0,
               info = paste("Recurrence rows without Day 0:", paste(missing, collapse = ", ")))
})

test_that("import_data: marker_info contains required columns", {
  expect_true(all(c("marker_id", "repeatlength", "binning_method") %in%
                    colnames(angola_imported_data$marker_info)))
})

test_that("import_data: marker_info only contains markers present in data", {
  marker_regex <- "(_allele_|_)\\d+$"
  cols_in_data <- gsub(marker_regex, "",
                       colnames(angola_imported_data$late_failures)[-(1:2)])
  expect_true(all(angola_imported_data$marker_info$marker_id %in%
                    unique(cols_in_data)))
})

test_that("import_data: allele columns are numeric for length_polymorphic data", {
  allele_df <- angola_imported_data$late_failures[, -(1:2), drop = FALSE]
  expect_true(all(sapply(allele_df, is.numeric)))
})

test_that("import_data: additional and late_failures have identical columns", {
  expect_equal(colnames(angola_imported_data$late_failures),
               colnames(angola_imported_data$additional))
})

test_that("import_data: no recurrence rows with entirely missing alleles", {
  recur_rows  <- angola_imported_data$late_failures[
    grepl(" recurrence$", angola_imported_data$late_failures$Sample.ID), ]
  allele_cols <- recur_rows[, -(1:2), drop = FALSE]
  expect_equal(sum(rowSums(!is.na(allele_cols)) == 0), 0)
})

# ============================================================
# Known-input checks — mock microsatellite data
# ============================================================

test_that("import_data: mock microsatellite data loads correctly", {
  files <- create_mock_xlsx(create_mock_microsat_data())
  on.exit(unlink(c(files$data, files$marker)))
  
  result <- import_data(files$data, files$marker, verbose = FALSE)
  
  expect_equal(result$data_type, "length_polymorphic")
  expect_equal(nrow(result$late_failures), 6)
  expect_equal(
    result$late_failures$Sample.ID,
    c("BD21-002 Day 0", "BD21-002 recurrence",
      "BD21-040 Day 0", "BD21-040 recurrence",
      "BD21-041 Day 0", "BD21-041 recurrence")
  )
})

test_that("import_data: mock data marker_info matches expected", {
  files <- create_mock_xlsx(create_mock_microsat_data())
  on.exit(unlink(c(files$data, files$marker)))
  
  result <- import_data(files$data, files$marker, verbose = FALSE)
  
  expect_equal(result$marker_info$marker_id,      c("TA1", "POLYA", "PFPK2", "TA109"))
  expect_equal(result$marker_info$binning_method, rep("microsatellite", 4))
  expect_equal(result$marker_info$repeatlength,   c(3, 3, 3, 3))
})

# ============================================================
# ID standardisation
# ============================================================

test_that("import_data: _D0/_DX suffix format standardised correctly", {
  files <- create_mock_xlsx(cbind(
    data.frame(Sample.ID = c("S1_D0", "S1_DX", "S2_D0", "S2_DX"),
               Site      = "A",
               stringsAsFactors = FALSE),
    mock_allele_cols(4)
  ))
  on.exit(unlink(c(files$data, files$marker)))
  
  ids <- import_data(files$data, files$marker, verbose = FALSE)$late_failures$Sample.ID
  
  expect_true(all(grepl(" Day 0$| recurrence$", ids)))
  expect_true(any(grepl(" Day 0$",      ids)))
  expect_true(any(grepl(" recurrence$", ids)))
})

test_that("import_data: numeric Day column format standardised correctly", {
  files <- create_mock_xlsx(cbind(
    data.frame(PatientID = c("U001", "U001", "U002", "U002"),
               Day       = c(0, 28, 0, 21),
               Site      = "A",
               stringsAsFactors = FALSE),
    mock_allele_cols(4)
  ))
  on.exit(unlink(c(files$data, files$marker)))
  
  ids <- import_data(files$data, files$marker, verbose = FALSE)$late_failures$Sample.ID
  
  expect_true(all(grepl(" Day 0$| recurrence$", ids)))
})

test_that("import_data: numeric day suffix D42 standardised to recurrence", {
  files <- create_mock_xlsx(cbind(
    data.frame(Sample.ID = c("BD21-002D0", "BD21-002D42"),
               Site      = "Benguela",
               stringsAsFactors = FALSE),
    mock_allele_cols(2)
  ))
  on.exit(unlink(c(files$data, files$marker)))
  
  ids <- import_data(files$data, files$marker, verbose = FALSE)$late_failures$Sample.ID
  
  expect_true(any(grepl(" Day 0$",      ids)))
  expect_true(any(grepl(" recurrence$", ids)))
})

# ============================================================
# Bad input — file path errors
# ============================================================

test_that("import_data: non-existent file throws error", {
  expect_error(
    import_data("nonexistent_file.xlsx", "any_marker.xlsx"),
    regexp = "Cannot read file"
  )
})

test_that("import_data: missing marker file throws error", {
  files <- create_mock_xlsx(create_mock_microsat_data())
  on.exit(unlink(c(files$data, files$marker)))
  
  expect_error(
    import_data(files$data, "nonexistent_markers.xlsx"),
    regexp = "Marker information not found"
  )
})

# ============================================================
# Bad input — malformed data
# ============================================================

test_that("import_data: fewer than 4 columns throws error", {
  files <- create_mock_xlsx(data.frame(ID = "S1", Site = "A", marker1 = 180))
  on.exit(unlink(c(files$data, files$marker)))
  
  expect_error(import_data(files$data, files$marker),
               regexp = "at least 4 columns")
})

test_that("import_data: all allele columns empty throws error", {
  files <- create_mock_xlsx(cbind(
    data.frame(Sample.ID = c("S1_D0", "S1_DX"),
               Site      = "A",
               stringsAsFactors = FALSE),
    data.frame(
      TA1_1   = c(NA_real_, NA_real_), TA1_2   = c(NA_real_, NA_real_),
      POLYA_1 = c(NA_real_, NA_real_), POLYA_2 = c(NA_real_, NA_real_),
      PFPK2_1 = c(NA_real_, NA_real_), PFPK2_2 = c(NA_real_, NA_real_),
      TA109_1 = c(NA_real_, NA_real_), TA109_2 = c(NA_real_, NA_real_)
    )
  ))
  on.exit(unlink(c(files$data, files$marker)))
  
  expect_error(import_data(files$data, files$marker),
               regexp = "All allele columns are empty")
})

# ============================================================
# Sheet 2 — additional background data (length_polymorphic)
# ============================================================

test_that("import_data: loads additional data from Sheet 2", {
  main_data <- create_mock_microsat_data()

  # Additional background samples — same column layout as main data
  additional_data <- data.frame(
    Sample.ID = c("BG-001 Day 0", "BG-001 recurrence"),
    Site      = "Benguela",
    TA1_1     = c(174, 174), TA1_2   = c(NA,  NA),
    POLYA_1   = c(162, 162), POLYA_2 = c(NA,  NA),
    PFPK2_1   = c(168, 168), PFPK2_2 = c(NA,  NA),
    TA109_1   = c(184, 184), TA109_2 = c(NA,  NA),
    check.names      = FALSE,
    stringsAsFactors = FALSE
  )

  tmp_data   <- tempfile(fileext = ".xlsx")
  tmp_marker <- tempfile(fileext = ".xlsx")
  writexl::write_xlsx(list(Sheet1 = main_data, Sheet2 = additional_data), tmp_data)
  writexl::write_xlsx(create_mock_markers(), tmp_marker)
  on.exit(unlink(c(tmp_data, tmp_marker)))

  result <- suppressMessages(suppressWarnings(
    import_data(tmp_data, tmp_marker, verbose = FALSE)
  ))

  expect_gt(nrow(result$additional), 0L)
  expect_equal(colnames(result$additional), colnames(result$late_failures))
})

test_that("import_data: no matching markers throws error", {
  files <- create_mock_xlsx(
    data = data.frame(
      Sample.ID = c("S1_D0", "S1_DX"),
      Site      = "A",
      unknown_1 = c(180, 180),
      unknown_2 = c(181, 181),
      unknown_3 = c(182, 182),
      unknown_4 = c(183, 183),
      stringsAsFactors = FALSE
    ),
    markers = data.frame(
      marker_id      = "POLYA",  # deliberately mismatches "unknown" columns
      markertype     = "microsatellite",
      binning_method = "microsatellite",
      repeatlength   = 3,
      stringsAsFactors = FALSE
    )
  )
  on.exit(unlink(c(files$data, files$marker)))
  
  expect_error(import_data(files$data, files$marker),
               regexp = "No matching markers found")
})