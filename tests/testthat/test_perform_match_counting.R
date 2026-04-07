library(testthat)
library(MalReBay)
# (create_mock_microsat_data, create_mock_markers)
#
# Verified match outcomes — all markers repeat=3:
#
#   BD21-002  R  / NI / NI / NI  -> Number_Matches=1, Number_Loci_Compared=4
#             TA1: 174 vs 177, diff=3 == repeat=3 -> R
#
#   BD21-040  NI / R  / NI / NI  -> Number_Matches=1, Number_Loci_Compared=4
#             TA1: 162 vs 171, diff=9 > repeat=3  -> NI
#
#   BD21-041  R  / R  / R  / R   -> Number_Matches=4, Number_Loci_Compared=4

# ============================================================
# Output structure
# ============================================================

test_that("perform_match_counting: returns correct columns", {
  result <- perform_match_counting(create_mock_microsat_data(), create_mock_markers())
  expect_true(all(
    c("Sample.ID", "Number_Matches", "Number_Loci_Compared",
      "TA1", "POLYA", "PFPK2", "TA109") %in% colnames(result)
  ))
})

test_that("perform_match_counting: one row per patient not per sample", {
  result <- perform_match_counting(create_mock_microsat_data(), create_mock_markers())
  expect_equal(nrow(result), 3)
  expect_true(all(c("BD21-002", "BD21-040", "BD21-041") %in% result$Sample.ID))
})

# ============================================================
# Per-locus outcomes — BD21-002
# ============================================================

test_that("perform_match_counting: BD21-002 matches only at TA1", {
  result  <- perform_match_counting(create_mock_microsat_data(), create_mock_markers())
  bd21_02 <- result[result$Sample.ID == "BD21-002", ]
  
  expect_equal(bd21_02$TA1,   "R")   # 174 vs 177, diff=3 == repeat=3
  expect_equal(bd21_02$POLYA, "NI")  # 153 vs 105, diff=48
  expect_equal(bd21_02$PFPK2, "NI")  # 159 vs 183, diff=24
  expect_equal(bd21_02$TA109, "NI")  # 172 vs 184, diff=12
  expect_equal(bd21_02$Number_Matches,       1)
  expect_equal(bd21_02$Number_Loci_Compared, 4)
})

# ============================================================
# Per-locus outcomes — BD21-040
# ============================================================

test_that("perform_match_counting: BD21-040 matches only at POLYA", {
  result  <- perform_match_counting(create_mock_microsat_data(), create_mock_markers())
  bd21_40 <- result[result$Sample.ID == "BD21-040", ]
  
  expect_equal(bd21_40$TA1,   "NI")  # 162 vs 171, diff=9 > repeat=3
  expect_equal(bd21_40$POLYA, "R")   # 159 vs 159, exact match
  expect_equal(bd21_40$PFPK2, "NI")  # 162/171 vs 183, no overlap within repeat=3
  expect_equal(bd21_40$TA109, "NI")  # 160 vs 178, diff=18
  expect_equal(bd21_40$Number_Matches,       1)
  expect_equal(bd21_40$Number_Loci_Compared, 4)
})

# ============================================================
# Per-locus outcomes — BD21-041
# ============================================================

test_that("perform_match_counting: BD21-041 matches at all 4 loci", {
  result  <- perform_match_counting(create_mock_microsat_data(), create_mock_markers())
  bd21_41 <- result[result$Sample.ID == "BD21-041", ]
  
  expect_equal(bd21_41$TA1,   "R")   # 177 vs 177
  expect_equal(bd21_41$POLYA, "R")   # 162/165 vs 162/165
  expect_equal(bd21_41$PFPK2, "R")   # 168/171 vs 168/171
  expect_equal(bd21_41$TA109, "R")   # 184 vs 184
  expect_equal(bd21_41$Number_Matches,       4)
  expect_equal(bd21_41$Number_Loci_Compared, 4)
})

# ============================================================
# Edge cases
# ============================================================

test_that("perform_match_counting: IND returned when recurrence alleles all NA", {
  mock <- data.frame(
    Sample.ID = c("S1 Day 0", "S1 recurrence"),
    Site      = "Benguela",
    TA1_1     = c(174,      NA), TA1_2   = c(NA, NA),
    POLYA_1   = c(153,      NA), POLYA_2 = c(NA, NA),
    PFPK2_1   = c(159,      NA), PFPK2_2 = c(NA, NA),
    TA109_1   = c(172,      NA), TA109_2 = c(NA, NA),
    check.names      = FALSE,
    stringsAsFactors = FALSE
  )
  result <- perform_match_counting(mock, create_mock_markers())
  
  expect_equal(result$TA1,   "IND")
  expect_equal(result$POLYA, "IND")
  expect_equal(result$PFPK2, "IND")
  expect_equal(result$TA109, "IND")
  expect_equal(result$Number_Loci_Compared, 0)
})

test_that("perform_match_counting: match at boundary of repeat length tolerance", {
  # diff = repeat length = 3 -> should be R (boundary is inclusive)
  mock <- data.frame(
    Sample.ID = c("S1 Day 0", "S1 recurrence"),
    Site      = "Benguela",
    TA1_1     = c(174, 177),  # diff = 3 == repeat = 3
    TA1_2     = c(NA,  NA),
    POLYA_1   = c(NA,  NA),   POLYA_2 = c(NA, NA),
    PFPK2_1   = c(NA,  NA),   PFPK2_2 = c(NA, NA),
    TA109_1   = c(NA,  NA),   TA109_2 = c(NA, NA),
    check.names      = FALSE,
    stringsAsFactors = FALSE
  )
  result <- perform_match_counting(mock, create_mock_markers())
  expect_equal(result$TA1, "R")
})