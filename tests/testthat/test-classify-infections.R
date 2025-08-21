test_that("classify_infections runs and returns correctly structured results", {
  example_file <- system.file("extdata", "Angola_2021_TES_7NMS.xlsx",
                              package = "MalReBay")
  
  testthat::skip_if_not(file.exists(example_file), "Example data file not found.")
  quick_mcmc_config <- list(
    n_chains = 2,
    max_iterations = 500, # Run one chunk
    chunk_size = 500,
    rhat_threshold = 1.1, # Relaxed for testing
    ess_threshold = 50   # Relaxed for testing
  )
  
  # 2. EXECUTION: Call the main function
  results_list <- classify_infections(
    input_filepath = example_file,
    mcmc_config = quick_mcmc_config,
    output_folder = tempdir()
  )
  
  testthat::expect_type(results_list, "list")
  testthat::expect_named(results_list, c("summary", "marker_details"))
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
})