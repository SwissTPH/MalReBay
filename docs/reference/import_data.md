# Import Genotyping Data from Excel

Reads data from a specified Excel file, automatically detecting the data
type and separating sheets into a structured list.

## Usage

``` r
import_data(
  filepath = system.file("extdata", "Dataset_microsatellite_panel.xlsx", package = "MalReBay"),
  marker_filepath = system.file("extdata", "makers_details.xlsx", package = "MalReBay"),
  additional_filepath = NULL,
  verbose = TRUE
)
```

## Arguments

- filepath:

  The full path to the input Excel file.

- marker_filepath:

  Path to Excel file containing marker metadata

- additional_filepath:

  Optional path to a separate additional/background data file (csv or
  xlsx). Ignored if the main file already has 2 sheets. (optional if
  marker_info sheet is present)

- verbose:

  Logical. If TRUE, prints progress and data-cleaning messages.

## Value

A list containing the imported data.

## Examples

``` r
if (FALSE) { # \dontrun{
  data_file <- system.file("extdata", 
                           "Dataset_microsatellite_panel.xlsx", 
                           package = "MalReBay")
  marker_file <- system.file("extdata", 
                             "makers_details.xlsx", 
                             package = "MalReBay")
  imported    <- import_data(filepath = data_file, 
                             marker_filepath = marker_file)
} # }
```
