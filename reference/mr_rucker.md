# MR Rucker framework

MR Rucker framework.

## Usage

``` r
mr_rucker(dat, parameters = default_parameters())
```

## Arguments

- dat:

  Output from
  [`harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.md).

- parameters:

  List of Qthresh for determining transition between models, and alpha
  values for calculating confidence intervals. Defaults to 0.05 for both
  in
  [`default_parameters()`](https://mrcieu.github.io/TwoSampleMR/reference/default_parameters.md).

## Value

list
