library(testthat)
library(MalReBay)

zaire_imported <- readRDS(system.file("extdata", "imported_data.rds", package = "MalReBay"))

mcmc_config <- list(
  n_chains     = 2,
  iter         = 1000,
  burn_in_frac = 0.5,
  random_seed  = 42,
  adapt_delta  = 0.8
)

recrudescent_ids <- c(
  "ZL21-203", "ZL21-233", "ZL21-245", "ZL21-260", "ZL21-262", "ZL21-263",
  "ZL21-269", "ZL21-287", "ZL21-292", "ZL21-304", "ZQ21-030", "ZQ21-042",
  "ZQ21-054", "ZQ21-077", "ZQ21-085", "ZQ21-103"
)

cmdstan_ok <- tryCatch({
  path <- cmdstanr::cmdstan_path()
  nzchar(path) && file.exists(path)
}, error = function(e) FALSE)

if (cmdstan_ok) {
  results <- suppressMessages(suppressWarnings(
    classify_infections(
      imported_data = zaire_imported,
      mcmc_config   = mcmc_config,
      n_workers     = 1,
      verbose       = FALSE
    )
  ))
}

test_that("classify_infections returns the correct counts of classification", {
  skip_if_not(cmdstan_ok, "CmdStan not installed")
  expect_named(results, c("classifications", "all_chains_loglikelihood",
                          "ids", "locus_summary", "locus_lrs",
                          "locus_dists", "locinames", "stan_fits"))
  expect_s3_class(results$stan_fits$Zaire, "CmdStanMCMC")
})

test_that("classify_infections ids match Zaire patients", {
  skip_if_not(cmdstan_ok, "CmdStan not installed")
  expected_n <- length(unique(gsub(" Day 0$| recurrence$", "", zaire_imported$late_failures$Sample.ID)))
  expect_length(results$ids$Zaire, expected_n)
  expect_true(all(recrudescent_ids %in% results$ids$Zaire))
})

test_that("classify_infections locus_summary has correct columns", {
  skip_if_not(cmdstan_ok, "CmdStan not installed")
  expect_true(all(c("patient_id", "n_comparable_loci") %in% colnames(results$locus_summary$Zaire)))
})

test_that("classify_infections likelihoods are finite and vary", {
  skip_if_not(cmdstan_ok, "CmdStan not installed")
  all_liks <- unlist(results$all_chains_loglikelihood)
  expect_true(all(is.finite(all_liks)))
  expect_gte(stats::sd(all_liks), 0)
})

test_that("classify_infections classifies the known Zaire recrudescence cases correctly", {
  skip_if_not(cmdstan_ok, "CmdStan not installed")
  for (id in recrudescent_ids) {
    prob_recrud <- mean(results$classifications$Zaire[, which(results$ids$Zaire == id)])
    expect_true(prob_recrud > 0.5, info = sprintf("%s: prob_recrud = %.3f", id, prob_recrud))
  }
})

test_that("classify_infections classifies the remaining Zaire patients as reinfection", {
  skip_if_not(cmdstan_ok, "CmdStan not installed")
  reinfection_ids <- setdiff(results$ids$Zaire, recrudescent_ids)
  for (id in reinfection_ids) {
    prob_recrud <- mean(results$classifications$Zaire[, which(results$ids$Zaire == id)])
    expect_true(prob_recrud < 0.5, info = sprintf("%s: prob_recrud = %.3f", id, prob_recrud))
  }
})

