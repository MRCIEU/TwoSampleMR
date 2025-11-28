# Penalised weighted median MR

Modification to standard weighted median MR Updated based on Burgess
2016 "Robust instrumental variable methods using multiple candidate
instruments with application to Mendelian randomization"

## Usage

``` r
mr_penalised_weighted_median(
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

  List containing `penk` - Constant term in penalisation, and `nboot` -
  number of bootstrap replications to calculate SE.
  [`default_parameters()`](https://mrcieu.github.io/TwoSampleMR/reference/default_parameters.md)
  sets `parameters=list(penk=20, nboot=1000)`.

## Value

List with the following elements:

- b:

  MR estimate

- se:

  Standard error

- pval:

  p-value
