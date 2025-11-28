# MR weighted mode estimator

MR weighted mode estimator

## Usage

``` r
mr_weighted_mode(
  b_exp,
  b_out,
  se_exp,
  se_out,
  parameters = default_parameters()
)
```

## Arguments

- b_exp:

  Vector of genetic effects on exposure

- b_out:

  Vector of genetic effects on outcome

- se_exp:

  Standard errors of genetic effects on exposure

- se_out:

  Standard errors of genetic effects on outcome

- parameters:

  List containing `phi` - Bandwidth parameter, and `nboot` - number of
  bootstraps to calculate SE.
  [`default_parameters()`](https://mrcieu.github.io/TwoSampleMR/reference/default_parameters.md)
  sets `list(phi=1, nboot=1000)`.

## Value

List with the following elements:

- b:

  MR estimate

- se:

  Standard error

- pval:

  p-value
