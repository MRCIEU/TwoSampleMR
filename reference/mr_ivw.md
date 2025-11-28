# Inverse variance weighted regression

The default multiplicative random effects IVW estimate. The standard
error is corrected for under dispersion Use the
[`mr_ivw_mre()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_ivw_mre.md)
function for estimates that don't correct for under dispersion.

## Usage

``` r
mr_ivw(b_exp, b_out, se_exp, se_out, parameters = default_parameters())
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
