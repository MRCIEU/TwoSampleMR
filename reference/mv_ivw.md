# Perform IVW multivariable MR

Performs modified multivariable MR analysis. For each exposure the
instruments are selected then all exposures for those SNPs are regressed
against the outcome together, weighting for the inverse variance of the
outcome.

## Usage

``` r
mv_ivw(mvdat, pval_threshold = 5e-08)
```

## Arguments

- mvdat:

  Output from
  [`mv_harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/mv_harmonise_data.md).

- pval_threshold:

  P-value threshold to include instruments. The default is `5e-8`.

## Value

List of results
