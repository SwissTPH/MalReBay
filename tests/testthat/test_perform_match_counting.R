library(testthat)
library(MalReBay)

zaire_imported <- readRDS(system.file("extdata", "imported_data.rds", package = "MalReBay"))

late      <- zaire_imported$late_failures
markers   <- zaire_imported$marker_info
locinames <- markers$marker_id
ids       <- unique(gsub(" Day 0$| recurrence$", "", late$Sample.ID))

result <- perform_match_counting(late, markers)

test_that("perform_match_counting returns correct columns and one row per patient", {
  expect_true(all(c("Sample.ID", "Number_Matches", "Number_Loci_Compared", locinames) %in% colnames(result)))
  expect_equal(nrow(result), length(ids))
  expect_true(all(ids %in% result$Sample.ID))
})

test_that("perform_match_counting matches independently-derived expectations for every Zaire patient and locus", {
  for (id in ids) {
    day0  <- late[late$Sample.ID == paste(id, "Day 0"), ]
    recur <- late[late$Sample.ID == paste(id, "recurrence"), ]
    row   <- result[result$Sample.ID == id, ]
    
    for (locus in locinames) {
      cols <- grep(paste0("^", locus, "_"), colnames(late), value = TRUE)
      d0   <- as.numeric(day0[, cols]);  d0 <- d0[!is.na(d0)]
      df   <- as.numeric(recur[, cols]); df <- df[!is.na(df)]
      repeat_length <- markers$repeatlength[markers$marker_id == locus]
      
      expected <- if (length(d0) == 0 || length(df) == 0) {
        "IND"
      } else if (any(sapply(df, function(x) any(abs(x - d0) <= repeat_length)))) {
        "R"
      } else {
        "NI"
      }
      
      expect_equal(row[[locus]], expected, info = sprintf("%s locus %s", id, locus))
    }
  }
})

test_that("perform_match_counting returns known outcomes for ZL21-218", {
  zl21_218 <- result[result$Sample.ID == "ZL21-218", ]
  expected <- c("313" = "IND", "383" = "NI", "TA1" = "IND", "POLYA" = "NI",
                "PFPK2" = "NI", "2490" = "IND", "TA109" = "NI")
  
  for (locus in names(expected)) {
    expect_equal(zl21_218[[locus]], unname(expected[locus]), info = locus)
  }
})

