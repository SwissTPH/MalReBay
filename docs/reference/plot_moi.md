# Plot Multiplicity of Infection (MOI)

This function calculates the MOI (number of distinct alleles per marker
per sample) and generates one violin plot per site for improved
readability, with a narrow boxplot (median/IQR) overlaid inside each
violin and jittered points for the individual samples.

## Usage

``` r
plot_moi(
  genotypedata,
  marker_pattern = "(_allele_\\d+|_\\d+)$",
  output_folder = NULL,
  filename_prefix = "moi_per_marker"
)
```

## Arguments

- genotypedata:

  A data frame containing `Sample.ID`, `Site`, and marker columns.

- marker_pattern:

  A regex pattern to identify marker columns. Defaults to standard
  allele suffix patterns.

- output_folder:

  Path to the directory where plots will be saved. If NULL, plots are
  not saved to disk.

- filename_prefix:

  Prefix for output PNG filenames. Each file will be named
  `<prefix>_<site>.png`.

## Value

A named list of ggplot objects, one per site. The underlying MOI data is
attached to each plot as the attribute `"moi_data"`.

## Examples

``` r
if (FALSE) { # \dontrun{
  gdata <- data.frame(
    Sample.ID = c("P1 Day 0", "P1 recurrence", "P2 Day 0", "P2 recurrence"),
    Site      = "SiteA",
    TA1_1     = c(174, 177, 162, 162),
    TA1_2     = c(177,  NA, 171,  NA)
  )
  plot_moi(genotypedata = gdata, output_folder = NULL)
} # }
```
