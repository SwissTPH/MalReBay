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

add_free_stan_data <- prepare_stan_data(
  late_failures_site  = late,
  additional_site     = add[0, ],
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

test_that("additional data feeds allele frequency counts, not classification matrices", {
  expect_equal(stan_data$recoded0, add_free_stan_data$recoded0)
  expect_equal(stan_data$recodedf, add_free_stan_data$recodedf)
  expect_gt(sum(stan_data$additional_counts), sum(add_free_stan_data$additional_counts))
})

test_that("late failures data feeds both classification matrices and allele frequency counts", {
  expect_true(all(late_free_stan_data$recoded0 == 0))
  expect_true(all(late_free_stan_data$recodedf == 0))
  expect_gt(sum(stan_data$additional_counts), sum(late_free_stan_data$additional_counts))
})
