# Simple median method

Perform MR using summary statistics. Bootstraps used to calculate
standard error.

## Usage

``` r
mr_simple_median(
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

  The number of bootstrap replications used to calculate the SE can be
  set through `parameters=list(nboot = 1000)`. The default is
  `list(nboot=1000)`.

## Value

List with the following elements:

- b:

  MR estimate

- se:

  Standard error

- pval:

  p-value

- nsnp:

  The number of SNPs
