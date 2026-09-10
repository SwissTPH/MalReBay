library(testthat)
library(MalReBay)

imported <- readRDS(system.file("extdata", "imported_data.rds", package = "MalReBay"))

late      <- imported$late_failures[, setdiff(colnames(imported$late_failures), "Site")]
add       <- imported$additional[, setdiff(colnames(imported$additional), "Site")]
markers   <- imported$marker_info
locinames <- markers$marker_id
ids       <- unique(gsub(" Day 0$| recurrence$", "", late$Sample.ID))

allele_cols <- setdiff(colnames(late), "Sample.ID")
maxMOI      <- max(as.integer(gsub(".*(_allele_|_)(\\d+)$", "\\2", allele_cols)))
allele_definitions <- suppressMessages(define_alleles(rbind(late, add), markers))
comparability <- compute_locus_comparability(late, ids, locinames)

stan_data <- prepare_stan_data(
  late_failures_site  = late,
  additional_site     = add,
  allele_definitions  = allele_definitions,
  marker_info         = markers,
  ids                 = ids,
  locinames           = locinames,
  maxMOI              = maxMOI,
  is_locus_comparable = comparability$is_locus_comparable
)

late_free_stan_data <- prepare_stan_data(
  late_failures_site  = late[0, ],
  additional_site     = add,
  allele_definitions  = allele_definitions,
  marker_info         = markers,
  ids                 = ids,
  locinames           = locinames,
  maxMOI              = maxMOI,
  is_locus_comparable = comparability$is_locus_comparable
)

test_that("additional_counts tallies every observed allele from late failures and additional data", {
  for (j in seq_along(locinames)) {
    cols <- grep(paste0("^", locinames[j], "_"), colnames(late), value = TRUE)
    
    trial_vals      <- as.numeric(unlist(late[, cols]))
    trial_vals      <- trial_vals[!is.na(trial_vals)]
    additional_vals <- as.numeric(unlist(add[, cols])) 
    additional_vals <- additional_vals[!is.na(additional_vals)]
    
    expect_equal(sum(stan_data$additional_counts[j, ]), length(trial_vals) + length(additional_vals))
  }
})

test_that("late failures data feeds both classification matrices and allele frequency counts", {
  expect_true(all(late_free_stan_data$recoded0 == 0))
  expect_true(all(late_free_stan_data$recodedf == 0))
  expect_gt(sum(stan_data$additional_counts), sum(late_free_stan_data$additional_counts))
})


test_that("MOI0 and MOIf reflect the true number of alleles observed per patient", {
  for (i in seq_along(ids)) {
    day0  <- late[late$Sample.ID == paste(ids[i], "Day 0"), ]
    recur <- late[late$Sample.ID == paste(ids[i], "recurrence"), ]
    
    n0 <- sapply(locinames, function(locus) {
      cols <- grep(paste0("^", locus, "_"), colnames(late), value = TRUE)
      sum(!is.na(day0[, cols]))
    })
    nf <- sapply(locinames, function(locus) {
      cols <- grep(paste0("^", locus, "_"), colnames(late), value = TRUE)
      sum(!is.na(recur[, cols]))
    })
    
    expect_equal(stan_data$MOI0[i], max(1L, n0))
    expect_equal(stan_data$MOIf[i], max(1L, nf))
  }
})

count_distinct_bins <- function(locus) {
  cols <- grep(paste0("^", locus, "_"), colnames(late), value = TRUE)
  vals <- c(as.numeric(unlist(late[, cols])), as.numeric(unlist(add[, cols])))
  vals <- vals[!is.na(vals)]
  bins <- sapply(vals, function(v) recodeallele(allele_definitions[[locus]], v))
  length(unique(bins[!is.na(bins)]))
}

expected_K <- setNames(sapply(locinames, count_distinct_bins), locinames)

test_that("K matches the real distinct allele count for every marker", {
  expect_equal(setNames(as.vector(stan_data$K), locinames), expected_K)
})

test_that("comparable correctly reflects which patient-locus pairs have real data at both timepoints", {
  for (i in seq_along(ids)) {
    day0  <- late[late$Sample.ID == paste(ids[i], "Day 0"), ]
    recur <- late[late$Sample.ID == paste(ids[i], "recurrence"), ]
    
    for (j in seq_along(locinames)) {
      cols     <- grep(paste0("^", locinames[j], "_"), colnames(late), value = TRUE)
      expected <- as.integer(any(!is.na(day0[, cols])) && any(!is.na(recur[, cols])))
      expect_equal(stan_data$comparable[i, j], expected)
    }
  }
})


test_that("prepare_stan_data output passes the package's own validation gate", {
  expect_true(suppressMessages(validate_stan_data(stan_data)))
})

test_that("prepare_stan_data returns non-empty outputs", {
  expect_true(any(stan_data$recoded0 != 0))
  expect_true(any(stan_data$recodedf != 0))
  expect_true(all(stan_data$K > 0))
  expect_true(any(stan_data$dist_array != 0))
  expect_true(any(stan_data$comparable == 1))
})

test_that("prepare_stan_data hidden alleles match recoded zeros", {
  expect_false(all(stan_data$hidden0 == 1))
  expect_false(all(stan_data$hiddenf == 1))
  expect_true(all(stan_data$hidden0 == (stan_data$recoded0 == 0)))
  expect_true(all(stan_data$hiddenf == (stan_data$recodedf == 0)))
})

test_that("prepare_stan_data MOI, method, and threshold are correct", {
  expect_true(all(stan_data$MOI0 >= 1 & stan_data$MOI0 <= maxMOI))
  expect_true(all(stan_data$MOIf >= 1 & stan_data$MOIf <= maxMOI))
  expect_true(all(stan_data$method_int == 1))
  expect_equal(as.numeric(stan_data$threshold), markers$repeatlength[match(locinames, markers$marker_id)])
})

