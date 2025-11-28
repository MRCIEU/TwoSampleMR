# Perform enrichment analysis

Perform enrichment analysis

## Usage

``` r
enrichment(dat, method_list = enrichment_method_list()$obj)
```

## Arguments

- dat:

  Harmonised exposure and outcome data. Output from
  [`harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.md).

- method_list:

  List of methods to use in analysis. Default is
  `enrichment_method_list()$obj`. See
  [`enrichment_method_list()`](https://mrcieu.github.io/TwoSampleMR/reference/enrichment_method_list.md)
  for details.

## Value

data frame
