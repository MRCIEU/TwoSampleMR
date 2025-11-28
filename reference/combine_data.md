# Combine data

Taking exposure or outcome data (returned from
[`format_data()`](https://mrcieu.github.io/TwoSampleMR/reference/format_data.md))
combine multiple datasets together so they can be analysed in one batch.
Removes duplicate SNPs, preferentially keeping those usable in MR
analysis.

## Usage

``` r
combine_data(x)
```

## Arguments

- x:

  List of data frames returned from
  [`format_data()`](https://mrcieu.github.io/TwoSampleMR/reference/format_data.md).

## Value

data frame
