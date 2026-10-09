# Check dataset viability per site

For each site, checks whether there are enough paired samples (Day 0 +
recurrence) with at least one comparable locus to be viable for
analysis.

## Usage

``` r
data_quality_check(
  imported_data,
  min_paired_samples = 10,
  min_comparable_loci = 1,
  verbose = TRUE
)
```

## Arguments

- imported_data:

  Output list from import_data().

- min_paired_samples:

  Minimum paired samples (each with at least min_comparable_loci
  comparable loci) required per site. Default 10.

- min_comparable_loci:

  Minimum comparable loci a paired sample must have to count toward
  min_paired_samples. Default 1.

- verbose:

  Logical. Warn for sites that fail the check.

## Value

A data frame: Site, n_paired, n_viable_pairs, viable.

## Examples

``` r
if (FALSE) { # \dontrun{
  data_file <- system.file("extdata", 
                           "Dataset_microsatellite_panel.xlsx", 
                           package = "MalReBay")
  marker_file <- system.file("extdata", 
                             "makers_details.xlsx", 
                             package = "MalReBay")
  imported <- import_data(filepath = data_file, 
                          marker_filepath = marker_file)
  quality  <- data_quality_check(imported)
} # }
```
