# Perform 2 sample MR on each SNP individually

Perform 2 sample MR on each SNP individually

## Usage

``` r
mr_singlesnp(
  dat,
  parameters = default_parameters(),
  single_method = "mr_wald_ratio",
  all_method = c("mr_ivw", "mr_egger_regression")
)
```

## Arguments

- dat:

  Output from
  [`harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.md).

- parameters:

  List of parameters. The default is
  [`default_parameters()`](https://mrcieu.github.io/TwoSampleMR/reference/default_parameters.md).

- single_method:

  Function to use for MR analysis. The default is `"mr_wald_ratio"`.

- all_method:

  Functions to use for MR analysis. The default is
  `c("mr_ivw", "mr_egger_regression")`.

## Value

List of data frames
