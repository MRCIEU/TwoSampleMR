# Perform MRMix analysis on harmonised dat object

See <https://github.com/gqi/MRMix> for more details.

## Usage

``` r
run_mrmix(dat)
```

## Arguments

- dat:

  Output from
  [`harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.md).
  Ensures that no eaf.exposure values are missing.

## Value

List of results, with one list item for every exposure/outcome pair in
dat object
