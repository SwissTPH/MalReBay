# Plot Match-Counting vs MalReBay Comparison Heatmap

Visualizes the WHO match-counting rules (loose and strict) alongside
MalReBay's posterior probability, one heatmap per site, from the
bayesian_match_counting_comparison table after it's been passed through
build_who_table().

## Usage

``` r
plot_comparison_heatmap(
  summary_results,
  marker_info,
  output_folder = NULL,
  title_prefix = "",
  verbose = TRUE
)
```

## Arguments

- summary_results:

  A list returned by
  [`summarise_results`](https://swisstph.github.io/MalReBay/reference/summarise_results.md),
  containing a `comparison` data frame (Bayesian probabilities plus
  traditional match-counting results).

- marker_info:

  Marker metadata from
  [`import_data`](https://swisstph.github.io/MalReBay/reference/import_data.md),
  passed on to `build_who_table()` to derive the WHO_loose/WHO_strict
  columns.

- output_folder:

  A string specifying the directory where one PNG per site will be
  saved. If `NULL`, heatmaps are drawn on the current graphics device
  instead.

- title_prefix:

  Optional prefix before the site name in each plot's title.

- verbose:

  Logical. If `TRUE`, prints a message when each file is saved.

## Value

Invisibly NULL; draws one heatmap per site as a side effect.
