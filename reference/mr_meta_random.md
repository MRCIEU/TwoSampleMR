# Perform 2 sample IV using random effects meta analysis and delta method for standard errors

Perform 2 sample IV using random effects meta analysis and delta method
for standard errors

## Usage

``` r
mr_meta_random(b_exp, b_out, se_exp, se_out, parameters)
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

  causal effect estimate

- se:

  standard error

- pval:

  p-value

- Q, Q_df, Q_pval:

  Heterogeneity stats
