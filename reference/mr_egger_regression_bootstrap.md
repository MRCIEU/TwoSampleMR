# Run bootstrap to generate standard errors for MR

Run bootstrap to generate standard errors for MR

## Usage

``` r
mr_egger_regression_bootstrap(b_exp, b_out, se_exp, se_out, parameters)
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

  List of parameters. Specifically, the `nboot` parameter can be
  specified for the number of bootstrap replications. The default is
  `parameters=list(nboot=1000)`.

## Value

List of with the following elements:

- b:

  MR estimate

- se:

  Standard error of MR estimate

- pval:

  p-value of MR estimate

- b_i:

  Estimate of horizontal pleiotropy (intercept)

- se_i:

  Standard error of intercept

- pval_i:

  p-value of intercept

- mod:

  Summary of regression

- dat:

  Original data used for MR Egger regression
