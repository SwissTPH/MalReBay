library(testthat)
library(MalReBay)

zaire_imported <- readRDS(system.file("extdata", "imported_data.rds", package = "MalReBay"))

late      <- zaire_imported$late_failures
markers   <- zaire_imported$marker_info
locinames <- markers$marker_id
ids       <- unique(gsub(" Day 0$| recurrence$", "", late$Sample.ID))

test_that("import_data returns correct list structure", {
  expect_type(imported, "list")
  expect_named(imported, c("late_failures", "additional", "marker_info", "data_type"))
})

test_that("import_data detects length_polymorphic Zaire data and converts alleles to numeric", {
  expect_equal(imported$data_type, "length_polymorphic")
  allele_cols <- late[, setdiff(colnames(late), c("Sample.ID", "Site"))]
  expect_true(all(sapply(allele_cols, is.numeric)))
})

test_that("import_data late_failures has required columns", {
  expect_true(all(c("Sample.ID", "Site") %in% colnames(late)))
  expect_gte(ncol(late), 3)
})

test_that("import_data Sample.ID contains only 'Day 0' or 'recurrence' suffixes", {
  expect_true(all(grepl(" Day 0$| recurrence$", late$Sample.ID)))
})

test_that("import_data recurrence rows are properly paired and non-empty", {
  day0_ids  <- gsub(" Day 0$",      "", late$Sample.ID[grepl(" Day 0$",      late$Sample.ID)])
  recur_ids <- gsub(" recurrence$", "", late$Sample.ID[grepl(" recurrence$", late$Sample.ID)])
  expect_equal(setdiff(recur_ids, day0_ids), character(0))
  
  recur_rows  <- late[grepl(" recurrence$", late$Sample.ID), ]
  allele_cols <- recur_rows[, setdiff(colnames(recur_rows), c("Sample.ID", "Site"))]
  expect_equal(sum(rowSums(!is.na(allele_cols)) == 0), 0)
})

test_that("import_data marker_info contains required columns", {
  expect_true(all(c("marker_id", "repeatlength", "binning_method") %in% colnames(markers)))
})

test_that("import_data marker_info only contains markers present in the data", {
  allele_cols  <- setdiff(colnames(late), c("Sample.ID", "Site"))
  cols_in_data <- unique(gsub("(_allele_|_)\\d+$", "", allele_cols))
  expect_true(all(markers$marker_id %in% cols_in_data))
})

test_that("import_data additional and late_failures have identical columns", {
  expect_equal(colnames(late), colnames(add))
})

test_that("import_data missing marker file throws error", {
  expect_error(import_data(data_file, "nonexistent_markers.xlsx"), regexp = "Marker information not found")
})

test_that("import_data no matching markers throws error", {
  bad_markers <- markers
  bad_markers$marker_id <- paste0("UNKNOWN_", bad_markers$marker_id)
  
  tmp_marker <- tempfile(fileext = ".xlsx")
  writexl::write_xlsx(bad_markers, tmp_marker)
  on.exit(unlink(tmp_marker))
  
  expect_error(import_data(data_file, tmp_marker), regexp = "No matching markers found")
})

test_that("import_data drops patients with no recurrence data and messages which ones", {
  broken <- late
  first_recur <- which(grepl(" recurrence$", broken$Sample.ID))[1]
  broken[first_recur, setdiff(colnames(broken), c("Sample.ID", "Site"))] <- NA
  
  tmp_data <- tempfile(fileext = ".xlsx")
  writexl::write_xlsx(list(Sheet1 = broken, Sheet2 = add), tmp_data)
  on.exit(unlink(tmp_data))
  
  expect_message(import_data(tmp_data, marker_file, verbose = TRUE), regexp = "Removing 1 patient")
})

test_that("data_quality_check confirms Zaire meets the minimum sample requirement", {
  quality       <- suppressWarnings(data_quality_check(imported, min_paired_samples = 10))
  zaire_quality <- quality[quality$Site == "Zaire", ]
  expect_true(zaire_quality$n_paired >= 10)
  expect_true(zaire_quality$viable)
})

test_that("detect_msp_variants finds no MSP families in Zaire's marker panel", {
  variants <- detect_msp_variants(markers$marker_id)
  expect_length(variants$msp1, 0)
  expect_length(variants$msp2, 0)
})
