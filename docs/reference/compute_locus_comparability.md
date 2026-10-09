# Compute locus comparability per patient

For each patient, determines which loci have non-missing data at both
Day 0 and recurrence timepoints ("comparable"), and tallies per-patient
availability. A locus is comparable only if at least one non-NA allele
is present at that locus for both timepoints.

## Usage

``` r
compute_locus_comparability(late_site, ids, locinames)
```

## Arguments

- late_site:

  The late failures data for a single site (Site column already
  removed), with Sample.ID and per-locus allele columns.

- ids:

  Character vector of patient IDs (without " Day 0"/" recurrence"
  suffixes) to check.

- locinames:

  Character vector of locus names to check.

## Value

A list with two elements:

- locus_summary:

  A data frame with one row per patient: patient_id, n_available_d0,
  n_available_df, n_comparable_loci.

- is_locus_comparable:

  A logical matrix (patients x loci) indicating which patient-locus
  pairs are comparable.
