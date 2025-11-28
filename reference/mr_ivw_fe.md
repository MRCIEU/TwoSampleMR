# Inverse variance weighted regression (fixed effects)

Inverse variance weighted regression (fixed effects)

## Usage

``` r
mr_ivw_fe(b_exp, b_out, se_exp, se_out, parameters = default_parameters())
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

List with the following elements:

- b:

  MR estimate

- se:

  Standard error

- pval:

  p-value

- Q, Q_df, Q_pval:

  Heterogeneity stats
