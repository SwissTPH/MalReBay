# Generate and Save Diversity Pie Charts

This function creates pie charts to visualize allele or haplotype
diversity at Day 0 and at the day of recurrence. It saves the combined
plot as a PDF file.

## Usage

``` r
generate_diversity_plots(
  genotypedata,
  data_type,
  marker_info = NULL,
  output_folder,
  filename_prefix = "diversity"
)
```

## Arguments

- genotypedata:

  A dataframe containing the genotyping data.

- data_type:

  A string, either "length_polymorphic" or "ampseq".

- marker_info:

  A dataframe with information about the markers, required if
  `data_type` is "length_polymorphic".

- output_folder:

  The path to the directory where the output PDF will be saved.

- filename_prefix:

  A string prefix for the output PDF filename.

## Value

Invisibly returns the final combined ggplot object.
