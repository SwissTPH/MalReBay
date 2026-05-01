library(testthat)
library(MalReBay)

# ============================================================
# define_alleles()
# ============================================================

test_that("define_alleles: microsatellite bins cover observed allele values", {
  defs <- suppressMessages(
    define_alleles(create_mock_microsat_data(), create_mock_markers())
  )
  expect_true(any(defs$POLYA[, "lower"] <= 153 & defs$POLYA[, "upper"] >= 153))
  expect_true(any(defs$POLYA[, "lower"] <= 105 & defs$POLYA[, "upper"] >= 105))
  expect_equal(colnames(defs$POLYA), c("lower", "upper"))
  expect_equal(defs$POLYA[, "lower"], sort(defs$POLYA[, "lower"]))
})

test_that("define_alleles: maxk filters to top k alleles", {
  defs_full <- suppressMessages(define_alleles(create_mock_microsat_data(), create_mock_markers()))
  defs_k2   <- suppressMessages(define_alleles(create_mock_microsat_data(), create_mock_markers(), maxk = 2))
  
  expect_lte(nrow(defs_k2$POLYA), 2)
  expect_gte(nrow(defs_full$POLYA), nrow(defs_k2$POLYA))
})

# ============================================================
# recodeallele()
# ============================================================

test_that("recodeallele: observed values map to valid and distinct bin indices", {
  bins <- suppressMessages(define_alleles(create_mock_microsat_data(), create_mock_markers()))$POLYA
  # With repeatlength = 3, expected bins from mock POLYA values (105, 153, 159, 162, 165):
  #   bin 1: ~103.5–106.5  (covers 105)
  #   bin 2: ~151.5–154.5  (covers 153)
  #   bin 3: ~157.5–163.5  (covers 159, 162)
  #   bin 4: ~163.5–166.5  (covers 165)
  
  idx_105 <- recodeallele(bins, 105)
  idx_153 <- recodeallele(bins, 153)
  
  expect_true(!is.na(idx_105))
  expect_true(!is.na(idx_153))
  expect_true(idx_105 != idx_153)
})

test_that("recodeallele: returns NA_integer_ for NA, empty bins, or beyond max distance", {
  bins       <- suppressMessages(define_alleles(create_mock_microsat_data(), create_mock_markers()))$POLYA
  empty_bins <- matrix(NA, ncol = 2, nrow = 0, dimnames = list(NULL, c("lower", "upper")))
  
  expect_equal(recodeallele(bins, NA),                            NA_integer_)
  expect_equal(recodeallele(empty_bins, 153),                     NA_integer_)
  expect_equal(recodeallele(bins, 300, max_distance_allowed = 5), NA_integer_)
})

# ============================================================
# calculate_frequencies()
# ============================================================
test_that("calculate_frequencies: microsatellite frequencies sum to 1", {
  defs  <- suppressMessages(define_alleles(create_mock_microsat_data(), create_mock_markers()))
  freqs <- calculate_frequencies(create_mock_microsat_data(), defs, create_mock_markers())

  expect_equal(sum(freqs$freq_matrix["POLYA", ]), 1, tolerance = 1e-6)
  expect_equal(unname(freqs$n_alleles["POLYA"]), nrow(defs$POLYA))
})

test_that("calculate_frequencies: exact method frequencies sum to 1", {
  # ampseq data: character haplotype strings, binning_method = "exact"
  # create_mock_ampseq_markers() uses binning_method="ampseq" which is the
  # import-level label; calculate_frequencies expects "exact" for this branch.
  exact_markers <- data.frame(
    marker_id      = c("cpmp", "cpp"),
    markertype     = "amplicon",
    binning_method = "exact",
    repeatlength   = c(NA_real_, NA_real_),
    stringsAsFactors = FALSE
  )
  freqs <- calculate_frequencies(
    genotypedata        = create_mock_ampseq_data(),
    alleles_definitions = list(),         # not used by exact branch
    marker_info         = exact_markers
  )

  # cpmp: HAPL_A appears 5×, HAPL_B 1× -> 2 distinct alleles, freqs sum to 1
  expect_equal(sum(freqs$freq_matrix["cpmp", ]), 1, tolerance = 1e-6)
  expect_equal(unname(freqs$n_alleles["cpmp"]), 2L)
})

test_that("calculate_frequencies: exact method variability is expected heterozygosity", {
  # cpmp: HAPL_A (5/6) and HAPL_B (1/6)
  #   EH = 1 - ((5/6)^2 + (1/6)^2) = 10/36 ≈ 0.2778
  exact_markers <- data.frame(
    marker_id      = c("cpmp", "cpp"),
    markertype     = "amplicon",
    binning_method = "exact",
    repeatlength   = c(NA_real_, NA_real_),
    stringsAsFactors = FALSE
  )
  freqs <- calculate_frequencies(
    genotypedata        = create_mock_ampseq_data(),
    alleles_definitions = list(),
    marker_info         = exact_markers
  )
  expected_eh_cpmp <- 1 - ((5/6)^2 + (1/6)^2)
  expect_equal(unname(freqs$variability["cpmp"]), expected_eh_cpmp, tolerance = 1e-6)
})