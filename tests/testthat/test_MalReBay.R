library(testthat)
library(MalReBay)

imported <- readRDS(system.file("extdata", "imported_data.rds", package = "MalReBay"))
marker_file    <- system.file("extdata", "makers_details.xlsx", package = "MalReBay")

zaire_file <- tempfile(fileext = ".xlsx")
writexl::write_xlsx(
  list(Sheet1 = imported$late_failures, Sheet2 = imported$additional),
  zaire_file
)

mcmc_config <- list(
  n_chains     = 2,
  iter         = 1000,
  burn_in_frac = 0.5,
  random_seed  = 42,
  adapt_delta  = 0.8
)

cmdstan_ok <- tryCatch({
  path      <- cmdstanr::cmdstan_path()
  model_ok  <- !inherits(
    try(instantiate::stan_package_model(name = "malrebay_model", package = "MalReBay"), silent = TRUE),
    "try-error"
  )
  nzchar(path) && file.exists(path) && model_ok
}, error = function(e) FALSE)

if (cmdstan_ok) {
  tmp_out <- tempfile()
  result  <- suppressMessages(suppressWarnings(
    MalReBay(
      filepath        = zaire_file,
      marker_filepath = marker_file,
      mcmc_config     = mcmc_config,
      output_folder   = tmp_out,
      n_workers       = 1,
      verbose         = FALSE
    )
  ))
}

test_that("MalReBay returns correct list structure for Zaire", {
  skip_if_not(cmdstan_ok, "CmdStan not installed")
  expect_named(result, c("posterior_probabilities", "comparison", "convergence", "mcmc_loglikelihoods"))
})

test_that("MalReBay probabilities are between 0 and 1", {
  skip_if_not(cmdstan_ok, "CmdStan not installed")
  probs <- result$posterior_probabilities$Probability
  expect_true(all(probs >= 0 & probs <= 1, na.rm = TRUE))
})

test_that("MalReBay saves all expected output files", {
  skip_if_not(cmdstan_ok, "CmdStan not installed")
  expect_true(file.exists(file.path(tmp_out, "posterior_probabilities.csv")))
  expect_true(file.exists(file.path(tmp_out, "bayesian_match_counting_comparison.csv")))
  expect_true(file.exists(file.path(tmp_out, "diversity_length_polymorphic_comparison.png")))
  expect_gt(length(list.files(tmp_out, pattern = "^moi_per_marker_.*\\.png$")), 0)
  expect_true(file.exists(file.path(tmp_out, "recrudescence_probability_histogram.png")))
  expect_true(dir.exists(file.path(tmp_out, "convergence_diagnosis")))
})

test_that("MalReBay non-existent filepath throws error", {
  expect_error(MalReBay(filepath = "nonexistent.xlsx"), regexp = "Cannot read file")
})

recrudescent_ids <- c(
  "ZL21-203", "ZL21-233", "ZL21-245", "ZL21-260", "ZL21-262", "ZL21-263",
  "ZL21-269", "ZL21-287", "ZL21-292", "ZL21-304", "ZQ21-030", "ZQ21-042",
  "ZQ21-054", "ZQ21-077", "ZQ21-085", "ZQ21-103"
)

test_that("MalReBay posterior_probabilities classify the known Zaire recrudescence cases correctly", {
  skip_if_not(cmdstan_ok, "CmdStan not installed")
  for (id in recrudescent_ids) {
    prob <- result$posterior_probabilities$Probability[result$posterior_probabilities$Sample.ID == id]
    expect_true(prob > 0.5, info = sprintf("%s: prob = %.3f", id, prob))
  }
})

test_that("summarise_results NULL mcmc_results throws error", {
  expect_error(summarise_results(NULL, zaire_imported), regexp = "mcmc_results.*is empty")
})

test_that("save_results invalid summary_results throws error", {
  expect_error(save_results(list(wrong = "structure")), regexp = "summary_results.*must be a valid list")
})

test_that("MalReBay convergence is NULL or a data.frame", {
  skip_if_not(cmdstan_ok, "CmdStan not installed")
  expect_true(is.null(result$convergence) || is.data.frame(result$convergence))
})

test_that("MalReBay saves the WHO comparison table", {
  skip_if_not(cmdstan_ok, "CmdStan not installed")
  expect_true(file.exists(file.path(tmp_out, "who_comparison_table.csv")))
})

test_that("saved CSVs contain correct columns", {
  skip_if_not(cmdstan_ok, "CmdStan not installed")
  pp <- utils::read.csv(file.path(tmp_out, "posterior_probabilities.csv"))
  ct <- utils::read.csv(file.path(tmp_out, "bayesian_match_counting_comparison.csv"))
  expect_true(all(c("Sample.ID", "Probability") %in% colnames(pp)))
  expect_true("Probability" %in% colnames(ct))
})

test_that("save_results returns named paths and creates a missing output folder", {
  skip_if_not(cmdstan_ok, "CmdStan not installed")
  fresh_dir <- file.path(tempdir(), "fresh_save_results_test")
  on.exit(unlink(fresh_dir, recursive = TRUE))
  
  paths <- save_results(result, imported_data = imported, output_folder = fresh_dir, verbose = FALSE)
  
  expect_true(dir.exists(fresh_dir))
  expect_type(paths, "character")
  expect_true(all(c("posterior_probabilities", "comparison") %in% names(paths)))
  expect_true(all(file.exists(paths)))
})
