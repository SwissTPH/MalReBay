# tests/testthat

Test files for the MalReBay package.

## Structure
- `helper.R`          — shared fixtures loaded automatically by testthat
- `test_1_import_data.R`           — import_data() tests
- `test_2_define_alleles.R`        — define_alleles(), recodeallele(), calculate_frequencies()
- `test_3_classify_infections.R`   — classify_infections() tests
- `test_4_summarise_and_save.R`    — summarise_results() and save_results()
- `test_5_perform_match_counting.R`— perform_match_counting() tests
- `test_6_MalReBay.R`              — full pipeline integration tests

## Mock data
All tests use in-memory mock data defined in `helper.R`.
Real Angola data is loaded from `inst/extdata/imported_data.rds`
which was pre-saved via import_data() to avoid file path issues
during devtools::test() before the package is installed.

## Notes
- All file output tests use tempfile() — nothing is written permanently here.
- The _snaps/ folder is not used (no snapshot tests).
- The results/ folder is not used (all outputs go to tempfile()).
