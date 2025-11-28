# Perform 2 sample IV using Wald ratio.

Perform 2 sample IV using Wald ratio.

## Usage

``` r
mr_wald_ratio(b_exp, b_out, se_exp, se_out, parameters)
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

- nsnp:

  1
