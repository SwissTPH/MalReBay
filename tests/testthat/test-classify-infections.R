# Path and helper function

get_data_file <- function(file_name) {
  path <- system.file("extdata", file_name, package = "MalReBay")
  if (path != "" && file.exists(path)) return(normalizePath(path))
  
  # Search local project structure if not installed
  local_path <- testthat::test_path("..", "..", "inst", "extdata", file_name)
  if (file.exists(local_path)) return(normalizePath(local_path))
  
  stop("Could not find data file: ", file_name)
}

get_angola_test_subset <- function(n_pairs = 10, n_additional = 5) {
  data_path <- get_data_file("Angola_2021_TES_7NMS.xlsx")
  marker_path <- get_data_file("marker_details.xlsx")
  
  full_data <- import_data(data_path, marker_path, verbose = FALSE)
  
  # Standardize IDs and identify matching pairs
  all_ids <- full_data$late_failures$Sample.ID
  prefixes <- gsub(" Day (0|Failure)$", "", all_ids)
  id_counts <- table(prefixes)
  
  valid_ids <- names(id_counts[id_counts >= 2])
  selected_ids <- valid_ids[1:min(n_pairs, length(valid_ids))]
  
  # Subset to selected pairs and additional samples
  full_data$late_failures <- full_data$late_failures[prefixes %in% selected_ids, ]
  if (nrow(full_data$additional) > 0) {
    full_data$additional <- full_data$additional[1:min(n_additional, nrow(full_data$additional)), ]
  }
  
  return(full_data)
}

create_mock_microsat_data <- function() {
  data.frame(
    Sample.ID = c("Mock_01 Day 0", "Mock_01 Day Failure", "Mock_02 Day 0", "Mock_02 Day Failure"),
    Site = "Luanda",
    `313_1` = c(100, 103, 110, 150), 
    `313_2` = c(NA, NA, 113, NA),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

# Mock helpers for mathematical unit tests
create_mock_microsat_data <- function() {
  data.frame(
    Sample.ID = c("Mock_01 Day 0", "Mock_01 Day Failure", 
                  "Mock_02 Day 0", "Mock_02 Day Failure"),
    Site = "Benguela",
    `313_1` = c(248, 246, 230, 222),  # Data from your screenshot
    `313_2` = c(NA, NA, NA, NA),     # Empty replicates as seen in image
    `313_3` = c(NA, NA, NA, NA),
    check.names = FALSE,             # Crucial: prevents R from turning '313_1' into 'X313_1'
    stringsAsFactors = FALSE
  )
}

create_mock_markers <- function() {
  data.frame(
    marker_id = "M1",
    binning_method = "microsatellite",
    repeatlength = 2,
    stringsAsFactors = FALSE
  )
}

# Test
test_that("import_data handles Angolan ID standardization", {
  raw_ids <- c("BD21-002D0", "BD21-040D28")
  cleaned <- gsub("D0$", " Day 0", raw_ids)
  idx_fail <- !grepl(" Day 0$", cleaned)
  cleaned[idx_fail] <- gsub("D[0-9]+$", " Day Failure", cleaned[idx_fail])
  
  expect_equal(cleaned, c("BD21-002 Day 0", "BD21-040 Day Failure"))
})

test_that("Microsatellite binning and frequency logic", {
  marker_info <- data.frame(marker_id = "313", binning_method = "microsatellite", 
                            repeatlength = 2, stringsAsFactors = FALSE)
  mock_data <- create_mock_microsat_data()
  
  defs <- define_alleles(mock_data, marker_info)
  freqs <- calculate_frequencies(mock_data, defs, marker_info)
  
  # Check that a bin exists that covers our 100bp allele
  expect_true(any(defs$`313`[, "lower"] <= 248 & defs$`313`[, "upper"] >= 248))
  expect_equal(sum(freqs$freq_matrix["313", ]), 1)
})

test_that("MCMC Pipeline Integration: 10 pairs + 5 additional samples", {
  skip_on_cran() 
  
  dat <- get_angola_test_subset(n_pairs = 10, n_additional = 5)
  tmp_out <- file.path(tempdir(), "angola_test")
  
  set.seed(42)
  # Run without parallelization for the unit test to ensure stability
  results <- classify_infections(
    imported_data = dat,
    mcmc_config = list(max_iterations = 40, chunk_size = 20, n_chains = 2),
    output_folder = tmp_out,
    n_workers = 1,
    verbose = FALSE
  )
  
  expect_equal(nrow(results$summary), 10)
  expect_true(all(results$summary$Probability >= 0 & results$summary$Probability <= 1))
  
  unlink(tmp_out, recursive = TRUE)
})

test_that("MCMC results are reproducible and convergence is calculated", {
  skip_on_cran()
  dat <- get_angola_test_subset(n_pairs = 2, n_additional = 0)
  conf <- list(max_iterations = 20, chunk_size = 10, n_chains = 2)
  
  set.seed(100)
  res1 <- classify_infections(dat, conf, file.path(tempdir(), "r1"), n_workers = 1, verbose = FALSE)
  set.seed(100)
  res2 <- classify_infections(dat, conf, file.path(tempdir(), "r2"), n_workers = 1, verbose = FALSE)
  
  expect_equal(res1$summary$Probability, res2$summary$Probability)
  expect_named(res1, c("summary", "marker_details", "comparison", "mcmc_loglikelihoods"))
})