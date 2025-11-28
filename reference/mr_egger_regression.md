# Egger's regression for Mendelian randomization

Egger's regression for Mendelian randomization

## Usage

``` r
mr_egger_regression(b_exp, b_out, se_exp, se_out, parameters)
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

  List of parameters.

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

- Q, Q_df, Q_pval:

  Heterogeneity stats

- mod:

  Summary of regression

- dat:

  Original data used for MR Egger regression
