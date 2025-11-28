# Combine all mr results

This function combines results of
[`mr()`](https://mrcieu.github.io/TwoSampleMR/reference/mr.md),
[`mr_heterogeneity()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_heterogeneity.md),
[`mr_pleiotropy_test()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_pleiotropy_test.md)
and
[`mr_singlesnp()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_singlesnp.md)
into a single data frame. It also merges the results with outcome study
level characteristics in
[`available_outcomes()`](https://mrcieu.github.io/TwoSampleMR/reference/available_outcomes.md).
If desired it also exponentiates results (e.g. if the user wants log
odds ratio converted into odds ratios with 95 percent confidence
intervals). The exposure and outcome columns from the output from
[`mr()`](https://mrcieu.github.io/TwoSampleMR/reference/mr.md) contain
both the trait names and trait ids. The `combine_all_mrresults()`
function splits these into separate columns by default.

## Usage

``` r
combine_all_mrresults(
  res,
  het,
  plt,
  sin,
  ao_slc = TRUE,
  Exp = FALSE,
  split.exposure = FALSE,
  split.outcome = FALSE
)
```

## Arguments

- res:

  Results from
  [`mr()`](https://mrcieu.github.io/TwoSampleMR/reference/mr.md).

- het:

  Results from
  [`mr_heterogeneity()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_heterogeneity.md).

- plt:

  Results from
  [`mr_pleiotropy_test()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_pleiotropy_test.md).

- sin:

  Results from
  [`mr_singlesnp()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_singlesnp.md).

- ao_slc:

  Logical; if set to `TRUE` then outcome study level characteristics are
  retrieved from
  [`available_outcomes()`](https://mrcieu.github.io/TwoSampleMR/reference/available_outcomes.md).
  Default is `TRUE`.

- Exp:

  Logical; if set to `TRUE` results are exponentiated. Useful if user
  wants log odds ratios expressed as odds ratios. Default is `FALSE`.

- split.exposure:

  Logical; if set to `TRUE` the exposure column is split into separate
  columns for the exposure name and exposure ID. Default is `FALSE`.

- split.outcome:

  Logical; if set to `TRUE` the outcome column is split into separate
  columns for the outcome name and outcome ID. Default is `FALSE`.

## Value

data frame
