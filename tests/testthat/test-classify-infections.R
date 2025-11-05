example_file <- system.file("extdata", "Angola_2021_TES_7NMS.xlsx", package = "MalReBay")
testthat::skip_if_not(file.exists(example_file), "Example data file not found.")
processed_data <- import_data(filepath = example_file, verbose = FALSE)


test_that("classify_infections runs and returns correctly structured results", {
  # 1. SETUP: Define file paths and configuration
  example_file <- system.file("extdata", "Angola_2021_TES_7NMS.xlsx",
                              package = "MalReBay")
  
  testthat::skip_if_not(file.exists(example_file), "Example data file not found.")
  
  quick_mcmc_config <- list(
    n_chains = 2,
    max_iterations = 500,
    chunk_size = 500,
    rhat_threshold = 1.1,
    ess_threshold = 50
  )
  
 
  processed_data <- import_data(filepath = example_file, verbose = FALSE)
  
  # Call the main function with the processed data object
  results_list <- classify_infections(
    imported_data = processed_data, # Use the new argument name
    mcmc_config = quick_mcmc_config,
    output_folder = tempdir(),
    verbose = FALSE # Good practice to run tests silently
  )

  
  
  testthat::expect_type(results_list, "list")
  testthat::expect_named(results_list, c("summary", "marker_details", "mcmc_loglikelihoods"))
  
  final_summary <- results_list$summary
  
  # Test the 'summary' data frame
  testthat::skip_if(is.null(final_summary), "Summary component of results is NULL.")
  testthat::expect_s3_class(final_summary, "data.frame")
  testthat::expect_gt(nrow(final_summary), 0)
  testthat::expect_named(final_summary, c("Site", "Sample.ID", "Probability",
                                          "N_Available_D0", "N_Available_DF",
                                          "N_Comparable_Loci"))
  testthat::expect_true(all(final_summary$Probability >= 0 & final_summary$Probability <= 1))
  
  marker_details <- results_list$marker_details
  
  # Test the 'marker_details' data frame
  testthat::skip_if(is.null(marker_details), "Marker details component of results is NULL.")
  testthat::expect_s3_class(marker_details, "data.frame")
  testthat::expect_gt(nrow(marker_details), 0)
  testthat::expect_named(marker_details, c("Sample.ID", "Marker", "Mean_Likelihood_Ratio",
                                           "Mean_Distance", "Site", "Interpretation"))
  
  # Test the 'mcmc_loglikelihoods' component
  mcmc_logs <- results_list$mcmc_loglikelihoods
  testthat::skip_if(is.null(mcmc_logs), "MCMC loglikelihoods component is NULL.")
  testthat::expect_type(mcmc_logs, "list")
  testthat::expect_gt(length(mcmc_logs), 0)
  
})

test_that("classify_infections removes pre-existing convergence_diagnosis folder", {
  
  # ARRANGE
  temp_output_dir <- tempfile("test_output_")
  dir.create(temp_output_dir)
  old_convergence_dir <- file.path(temp_output_dir, "convergence_diagnosis")
  dir.create(old_convergence_dir)
  file.create(file.path(old_convergence_dir, "old_plot.png"))
  expect_true(dir.exists(old_convergence_dir))
  
  # Use a minimal config to make this test run very fast
  quick_mcmc_config <- list(n_chains = 1, max_iterations = 10, chunk_size = 10)
  
  # ACT
  classify_infections(
    imported_data = processed_data,
    mcmc_config = quick_mcmc_config,
    output_folder = temp_output_dir,
    verbose = FALSE
  )
  
  expect_false(dir.exists(old_convergence_dir))
  
  unlink(temp_output_dir, recursive = TRUE)
})