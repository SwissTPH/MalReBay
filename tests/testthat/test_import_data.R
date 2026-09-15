library(testthat)
library(MalReBay)

imported <- readRDS(system.file("extdata", "imported_data.rds", package = "MalReBay"))

late      <- imported$late_failures
add       <- imported$additional
markers   <- imported$marker_info
locinames <- markers$marker_id
ids       <- unique(gsub(" Day 0$| recurrence$", "", late$Sample.ID))

test_that("import_data returns correct list structure", {
  expect_type(imported, "list")
  expect_named(imported, c("late_failures", "additional", "marker_info", "data_type"))
})

test_that("import_data detects length_polymorphic and converts alleles to numeric", {
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

test_that("detect patients with no recurrence data in `late`", {
  missing_observation <- late
  
  first_recur <- which(grepl("\\brecurrence\\s*$", missing_observation$Sample.ID, ignore.case = TRUE))[1]
  expect_true(!is.na(first_recur), info = "no 'recurrence' Sample.ID found in `late` dataset")
  cols_to_na <- setdiff(colnames(missing_observation), c("Sample.ID", "Site"))
  missing_observation[first_recur, cols_to_na] <- NA
  
  base_ids <- unique(gsub(" Day 0$| recurrence$", "", missing_observation$Sample.ID))
  is_recur_missing <- vapply(base_ids, function(id) {
    recur_idx <- which(missing_observation$Sample.ID == paste0(id, " recurrence"))
    if (length(recur_idx) == 0) return(TRUE)
    all_na <- all(is.na(missing_observation[recur_idx, cols_to_na, drop = FALSE]))
    all_na
  }, logical(1))
  
  expect_equal(sum(is_recur_missing), 1)
  
  mutated_base <- gsub(" (Day 0|recurrence)$", "", missing_observation$Sample.ID[first_recur])
  expect_true(mutated_base %in% base_ids[is_recur_missing])
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
