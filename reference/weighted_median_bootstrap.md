# Calculate standard errors for weighted median method using bootstrap

Based on new script for weighted median confidence interval, update 31
July 2015.

## Usage

``` r
weighted_median_bootstrap(b_exp, b_out, se_exp, se_out, weights, nboot)
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

- weights:

  Weights to apply to each SNP.

- nboot:

  Number of bootstrap replications. The default is `1000`.

## Value

Empirical standard error
