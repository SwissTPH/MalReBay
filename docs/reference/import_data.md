# Import Genotyping Data from Excel

Reads data from a specified Excel file, automatically detecting the data
type and separating sheets into a structured list.

## Usage

``` r
import_data(filepath, verbose = TRUE)
```

## Arguments

- filepath:

  The full path to the input Excel file.

- verbose:

  Logical. If TRUE, prints progress and data-cleaning messages.

## Value

A list containing the imported data.
