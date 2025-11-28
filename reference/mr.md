# Perform all Mendelian randomization tests

Perform all Mendelian randomization tests

## Usage

``` r
mr(
  dat,
  parameters = default_parameters(),
  method_list = subset(mr_method_list(), use_by_default)$obj
)
```

## Arguments

- dat:

  Harmonised exposure and outcome data. Output from
  [`harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.md).

- parameters:

  Parameters to be used for various MR methods. Default is output from
  [`default_parameters()`](https://mrcieu.github.io/TwoSampleMR/reference/default_parameters.md).

- method_list:

  List of methods to use in analysis. See
  [`mr_method_list()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_method_list.md)
  for details.

## Value

List with the following elements:

- mr:

  Table of MR results

- extra:

  Table of extra results
