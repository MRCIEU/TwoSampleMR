# Weighted median method

Perform MR using summary statistics. Bootstraps used to calculate
standard error.

## Usage

``` r
mr_weighted_median(
  b_exp,
  b_out,
  se_exp,
  se_out,
  parameters = default_parameters()
)
```

## Arguments

- b_exp:

  Vector of genetic effects on exposure.

- b_out:

  Vector of genetic effects on outcome.

- se_exp:

  Standard errors of genetic effects on exposure.

- se_out:

  Standard errors of genetic effects on outcome.

- parameters:

  The default is
  [`default_parameters()`](https://mrcieu.github.io/TwoSampleMR/reference/default_parameters.md).
  Specify the number of bootstrap replications to calculate the SE with
  `nboot`. The default is `list(nboot=1000)`.

## Value

List with the following elements:

- b:

  MR estimate

- se:

  Standard error

- pval:

  p-value
