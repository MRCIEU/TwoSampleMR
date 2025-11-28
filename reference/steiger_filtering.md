# Steiger filtering function

This function takes an object from
[`harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.md)
and does the following: If there is no rsq.exposure or rsq.outcome
column it will try to estimate it. This is done differently for traits
that have "log odds" units. To estimate rsq for quantitative traits
there must be either p-values and sample sizes for each SNP, or effect
sizes and standard errors AND the units are in SD units (the column must
contain "SD"). To estimate rsq for binary traits the units must be
called "log odds" and there must be beta.exposure, eaf.exposure,
ncase.exposure, ncontrol.exposure, prevalence.exposure. The same
principles apply for calculating the rsq for the outcome trait, except
column names are beta.outcome etc. If prevalence isn't supplied then it
uses 0.1 by default.

## Usage

``` r
steiger_filtering(dat)
```

## Arguments

- dat:

  Output from
  [`harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.md).

## Value

[`harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.md)
style data frame with additional columns rsq.exposure, rsq.outcome,
steiger_dir (which is `TRUE` if the rsq.exposure is larger than
rsq.outcome) and steiger_pval which is a test to see if rsq.exposure is
significantly larger than rsq.outcome.

## Details

Once rsq is calculated for the exposure and outcome, it will then
perform the Steiger test for each SNP to see if the rsq of the exposure
is larger than the rsq of the outcome.

Note that Steiger filtering, while useful, does have its own pitfalls.
Try to use replication effect estimates for the exposure (which are not
biased by winner's curse), and note that if there is strong antagonistic
confounding or differential measurement error between the exposure and
outcome then the causal directions could be inferred incorrectly.
