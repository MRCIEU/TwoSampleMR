# Perform basic multivariable MR

Performs initial multivariable MR analysis from Burgess et al 2015. For
each exposure the outcome is residualised for all the other exposures,
then unweighted regression is applied.

## Usage

``` r
mv_basic(mvdat, pval_threshold = 5e-08)
```

## Arguments

- mvdat:

  Output from
  [`mv_harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/mv_harmonise_data.md).

- pval_threshold:

  P-value threshold to include instruments. The default is `5e-8`.

## Value

List of results
