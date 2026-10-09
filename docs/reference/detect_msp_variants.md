# Identify which locinames belong to the MSP1/MSP2 families

Cross-checks locinames against the standard MSP1 (K1, MAD20, RO33) and
MSP2 (3D7, FC27, IC) allelic family names to identify which markers in a
dataset belong to which family. Digit "0" and letter "O" are treated as
equivalent (e.g. "RO33" vs "R033"), matching ignores case/whitespace.

## Usage

``` r
detect_msp_variants(locinames)
```

## Arguments

- locinames:

  Character vector of marker/locus names present in the dataset.

## Value

A list with msp1 and msp2 elements: the subset of locinames (original
spelling, not normalized) matching each family's known variants.
