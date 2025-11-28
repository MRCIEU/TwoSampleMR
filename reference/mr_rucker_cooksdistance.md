# MR Rucker with outliers automatically detected and removed

Uses Cook's distance D \> 4/nsnp to iteratively remove outliers.

## Usage

``` r
mr_rucker_cooksdistance(dat, parameters = default_parameters())
```

## Arguments

- dat:

  Output from
  [`harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.md).

- parameters:

  List of parameters. The default is
  [`default_parameters()`](https://mrcieu.github.io/TwoSampleMR/reference/default_parameters.md).

## Value

List
