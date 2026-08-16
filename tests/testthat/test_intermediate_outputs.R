library(testthat)
library(MalReBay)

zaire_imported <- readRDS(system.file("extdata", "imported_data.rds", package = "MalReBay"))

late      <- zaire_imported$late_failures[, setdiff(colnames(zaire_imported$late_failures), "Site")]
add       <- zaire_imported$additional[, setdiff(colnames(zaire_imported$additional), "Site")]
markers   <- zaire_imported$marker_info
locinames <- markers$marker_id
ids       <- unique(gsub(" Day 0$| recurrence$", "", late$Sample.ID))

allele_cols <- setdiff(colnames(late), "Sample.ID")
maxMOI      <- max(as.integer(gsub(".*(_allele_|_)(\\d+)$", "\\2", allele_cols)))

test_that("import_data returns non-empty Zaire data", {
  expect_gt(nrow(late), 0)
  expect_gt(nrow(add), 0)
  expect_gt(nrow(markers), 0)
  expect_true(all(colSums(!is.na(late[, allele_cols])) > 0))
})

comparability <- compute_locus_comparability(late, ids, locinames)

test_that("compute_locus_comparability returns non-empty Zaire results", {
  expect_equal(nrow(comparability$locus_summary), length(ids))
  expect_true(any(comparability$is_locus_comparable))
})

allele_definitions <- suppressMessages(define_alleles(rbind(late, add), markers))

test_that("define_alleles returns non-empty bins for every Zaire locus", {
  for (locus in locinames) expect_gt(nrow(allele_definitions[[locus]]), 0)
})

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

test_that("prepare_stan_data output passes the package's own validation gate", {
  expect_true(suppressMessages(validate_stan_data(stan_data)))
})

test_that("prepare_stan_data returns non-empty Zaire outputs", {
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

test_that("observed alleles belong to the defined allele population", {
  for (j in seq_along(locinames)) {
    cols <- (maxMOI * (j - 1) + 1):(maxMOI * (j - 1) + maxMOI)
    observed0 <- stan_data$recoded0[, cols][stan_data$hidden0[, cols] == 0]
    observedf <- stan_data$recodedf[, cols][stan_data$hiddenf[, cols] == 0]
    expect_true(all(observed0 >= 1 & observed0 <= stan_data$K[j]))
    expect_true(all(observedf >= 1 & observedf <= stan_data$K[j]))
  }
})

test_that("prepare_stan_data MOI, method, and threshold are correct for Zaire", {
  expect_true(all(stan_data$MOI0 >= 1 & stan_data$MOI0 <= maxMOI))
  expect_true(all(stan_data$MOIf >= 1 & stan_data$MOIf <= maxMOI))
  expect_true(all(stan_data$method_int == 1))
  expect_equal(as.numeric(stan_data$threshold), markers$repeatlength[match(locinames, markers$marker_id)])
})
