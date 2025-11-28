# Harmonise LD matrix against summary data

LD matrix returns with rsid_ea_oa identifiers. Make sure that they are
oriented to the same effect allele as the summary dataset. Summary
dataset can be exposure dataset or harmonised dataset.

## Usage

``` r
harmonise_ld_dat(x, ld)
```

## Arguments

- x:

  Exposure dataset or harmonised dataset

- ld:

  Output from
  [`ld_matrix()`](https://mrcieu.github.io/TwoSampleMR/reference/ld_matrix.md)

## Value

List of exposure dataset and harmonised LD matrix
