# Perform multivariable MR on subset of features

The function proceeds as follows:

1.  Select features (by default this is done using LASSO feature
    selection).

2.  Subset the mvdat to only retain relevant features and instruments.

3.  Perform MVMR on remaining data.

## Usage

``` r
mv_subset(
  mvdat,
  features = mv_lasso_feature_selection(mvdat),
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

- features:

  Dataframe of features to retain, must have column with name 'exposure'
  that has list of exposures to retain from mvdat. The default is
  `mvdat_lasso_feature_selection(mvdat)`.

- intercept:

  Should the intercept by estimated (`TRUE`) or force line through the
  origin (`FALSE`, the default).

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
