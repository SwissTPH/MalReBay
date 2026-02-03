# Import Genotyping Data from Excel

Reads data from a specified Excel file, automatically detecting the data
type and separating sheets into a structured list.

## Usage

``` r
import_data(filepath, marker_filepath = NULL, verbose = TRUE)
```

## Arguments

- filepath:

  The full path to the input Excel file.

- marker_filepath:

  Path to Excel file containing marker metadata (optional if marker_info
  sheet is present)

- verbose:

  Logical. If TRUE, prints progress and data-cleaning messages.

## Value

A list containing the imported data.
