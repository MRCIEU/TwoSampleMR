# Maximum likelihood MR method

Maximum likelihood MR method

## Usage

``` r
mr_two_sample_ml(b_exp, b_out, se_exp, se_out, parameters)
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
