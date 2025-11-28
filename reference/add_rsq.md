# Estimate r-squared of each association

Can be applied to exposure_dat, outcome_dat or harmonised_data. Note
that it will be beneficial in some circumstances to add the meta data to
the data object using
[`add_metadata()`](https://mrcieu.github.io/TwoSampleMR/reference/add_metadata.md)
before running this function. Also adds effective sample size for case
control data.

## Usage

``` r
add_rsq(dat)
```

## Arguments

- dat:

  exposure_dat, outcome_dat or harmonised_data

## Value

data frame
