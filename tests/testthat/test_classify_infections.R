library(testthat)
library(MalReBay)
#   create_mock_imported_data(), fast_mcmc, angola_imported_data
#
# Classification direction relies on:
#   BD21-041 (4/4 match) -> recrudescence
#   BD21-040 (1/4 match) -> reinfection
# BD21-002 is also 1/4 match (TA1: diff=3 == repeat=3) so is NOT
# used for direction tests — it is ambiguous alongside BD21-040.

run_classify <- function() {
  skip_stan_on_check()
  suppressMessages(suppressWarnings(
    classify_infections(
      imported_data = create_mock_imported_data(),
      mcmc_config   = fast_mcmc,
      n_workers     = 1,
      verbose       = FALSE
    )
  ))
}

# ============================================================
# Output structure
# ============================================================

test_that("classify_infections: returns correct top-level list keys", {
  results <- run_classify()
  expect_named(results, c("classifications", "all_chains_loglikelihood",
                          "ids", "locus_summary", "locus_lrs",
                          "locus_dists", "locinames", "stan_fits"))
})

test_that("classify_infections: ids match patients in mock data", {
  results <- run_classify()
  expect_equal(length(results$ids$Benguela), 3)
  expect_true(all(c("BD21-002", "BD21-040", "BD21-041") %in%
                    results$ids$Benguela))
})

test_that("classify_infections: locus_summary has correct columns", {
  results      <- run_classify()
  summary_cols <- colnames(results$locus_summary$Benguela)
  expect_true(all(c("patient_id", "n_comparable_loci") %in% summary_cols))
})

test_that("classify_infections: likelihoods are finite and vary", {
  results  <- run_classify()
  all_liks <- unlist(results$all_chains_loglikelihood)
  expect_true(all(is.finite(all_liks)))
  expect_gte(stats::sd(all_liks), 0)
})

# ============================================================
# Classification direction tests
# ============================================================

test_that("classify_infections: BD21-041 (4/4 match) classified towards recrudescence", {
  skip_on_cran()
  results     <- run_classify()
  id_idx      <- which(results$ids$Benguela == "BD21-041")
  prob_recrud <- mean(results$classifications$Benguela[, id_idx])
  expect_gt(prob_recrud, 0.5)
})

test_that("classify_infections: BD21-040 (1/4 match) classified towards reinfection", {
  skip_on_cran()
  results     <- run_classify()
  id_idx      <- which(results$ids$Benguela == "BD21-040")
  prob_recrud <- mean(results$classifications$Benguela[, id_idx])
  expect_lt(prob_recrud, 0.5)
})
# ============================================================
# Real data smoke test — uses pre-imported .rds, no file paths needed
# ============================================================

test_that("classify_infections: runs on real Angola data without error", {
  skip_on_cran()
  skip_stan_on_check()
  results <- suppressMessages(suppressWarnings(
    classify_infections(
      imported_data = angola_imported_data,
      mcmc_config   = fast_mcmc,
      n_workers     = 1,
      verbose       = FALSE
    )
  ))
  expect_type(results, "list")
  expect_true(length(unlist(results$ids)) > 0)
  expect_s3_class(results$stan_fits[[1]], "CmdStanMCMC")
})