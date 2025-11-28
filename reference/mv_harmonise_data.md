# Harmonise exposure and outcome for multivariable MR

Harmonise exposure and outcome for multivariable MR

## Usage

``` r
mv_harmonise_data(exposure_dat, outcome_dat, harmonise_strictness = 2)
```

## Arguments

- exposure_dat:

  Output from
  [`mv_extract_exposures()`](https://mrcieu.github.io/TwoSampleMR/reference/mv_extract_exposures.md).

- outcome_dat:

  Output from `extract_outcome_data(exposure_dat$SNP, id_output)`.

- harmonise_strictness:

  See the `action` option of
  [`harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.md).
  The default is `2`.

## Value

List of vectors and matrices required for mv analysis.

- exposure_beta:

  a matrix of beta coefficients, in which rows correspond to SNPs and
  columns correspond to exposures.

- exposure_se:

  is the same as `exposure_beta`, but for standard errors.

- exposure_pval:

  the same as `exposure_beta`, but for p-values.

- expname:

  A data frame with two variables, `id.exposure` and `exposure` which
  are character strings.

- outcome_beta:

  an array of effects for the outcome, corresponding to the SNPs in
  `exposure_beta`.

- outcome_se:

  an array of standard errors for the outcome.

- outcome_pval:

  an array of p-values for the outcome.

- outname:

  A data frame with two variables, `id.outcome` and `outcome` which are
  character strings.
