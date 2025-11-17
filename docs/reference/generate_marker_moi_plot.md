# Generate and Save an MOI Plot

Generate and Save an MOI Plot

## Usage

``` r
generate_marker_moi_plot(
  moi_data,
  output_folder,
  filename = "moi_per_marker_by_site.pdf"
)
```

## Arguments

- moi_data:

  A data frame of MOI values, typically from `calculate_marker_moi`.

- output_folder:

  Path to the directory for saving the plot.

- filename:

  The name for the output PDF file.

## Value

The ggplot object of the MOI plot.
