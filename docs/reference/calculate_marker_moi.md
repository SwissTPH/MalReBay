# Calculate Marker-Specific MOI

Calculate Marker-Specific MOI

## Usage

``` r
calculate_marker_moi(genotypedata, marker_pattern = "_\\d+$")
```

## Arguments

- genotypedata:

  A data frame with `Sample.ID`, `Site`, and marker columns.

- marker_pattern:

  A regex pattern to identify marker columns.

## Value

A data frame summarizing MOI for each sample and marker.
