# Perform IVW multivariable MR

Performs modified multivariable MR analysis. For each exposure the
instruments are selected then all exposures for those SNPs are regressed
against the outcome together, weighting for the inverse variance of the
outcome.

## Usage

``` r
mv_multiple(
  mvdat,
  intercept = FALSE,
  instrument_specific = FALSE,
  pval_threshold = 5e-08,
  plots = FALSE
)
```

## Arguments

- mvdat:

  Output from
  [`mv_harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/mv_harmonise_data.md).

- intercept:

  Should the intercept by estimated (`TRUE`) or force line through the
  origin (`FALSE`, default).

- instrument_specific:

  Should the estimate for each exposure be obtained by using all
  instruments from all exposures (`FALSE`, default) or by using only the
  instruments specific to each exposure (`TRUE`).

- pval_threshold:

  P-value threshold to include instruments. The default is `5e-8`.

- plots:

  Create plots? The default is `FALSE`.

## Value

List of results
