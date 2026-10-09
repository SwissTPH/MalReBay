# Generate and Save Diversity Pie Charts

Creates pie charts visualising allele or haplotype diversity, pooling
Day 0 and recurrence samples: one figure for all sites combined and one
figure per site (when a `Site` column is present). Pies are arranged
with at most `max_cols` per row.

## Usage

``` r
plot_markers_diversity(
  genotypedata,
  data_type,
  marker_info = NULL,
  output_folder = NULL,
  filename_prefix = "diversity",
  max_cols = 4
)
```

## Arguments

- genotypedata:

  A dataframe containing the genotyping data.

- data_type:

  A string: `"length_polymorphic"` or `"ampseq"`.

- marker_info:

  A dataframe with marker definitions; required when
  `data_type = "length_polymorphic"`.

- output_folder:

  Path to the directory where the output PNGs will be saved:
  `<prefix>_<data_type>_comparison.png` for all sites and
  `<prefix>_<data_type>_<site>.png` per site. If `NULL`, nothing is
  saved.

- filename_prefix:

  A string prefix for the output filenames.

- max_cols:

  Maximum number of pie charts per row. Defaults to 4.

## Value

Invisibly returns a list with `all_sites` (one figure pooling all sites)
and `by_site` (a list of figures named by site; `NULL` if there is no
`Site` column), or `NULL` if no data.

## Examples

``` r
if (FALSE) { # \dontrun{
  gdata <- data.frame(
    Sample.ID     = c("P1 Day 0", "P1 recurrence", "P2 Day 0", "P2 recurrence"),
    Site          = "SiteA",
    cpmp_allele_1 = c("HAPL_A", "HAPL_A", "HAPL_A", "HAPL_B"),
    cpmp_allele_2 = c(NA, NA, NA, NA)
  )
  plot_markers_diversity(genotypedata = gdata, data_type = "ampseq")
} # }
```
