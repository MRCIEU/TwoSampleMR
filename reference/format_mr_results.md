# Format MR results for forest plot

This function takes the results from
[`mr()`](https://mrcieu.github.io/TwoSampleMR/reference/mr.md) and is
particularly useful if the MR has been applied using multiple exposures
and multiple outcomes. It creates a new data frame with the following:

- Variables: exposure, outcome, category, outcome sample size, effect,
  upper ci, lower ci, pval, nsnp

- only one estimate for each exposure-outcome

- exponentiated effects if required

## Usage

``` r
format_mr_results(
  mr_res,
  exponentiate = FALSE,
  single_snp_method = "Wald ratio",
  multi_snp_method = "Inverse variance weighted",
  ao_slc = TRUE,
  priority = "Cardiometabolic"
)
```

## Arguments

- mr_res:

  Results from
  [`mr()`](https://mrcieu.github.io/TwoSampleMR/reference/mr.md).

- exponentiate:

  Convert effects to OR? The default is `FALSE`.

- single_snp_method:

  Which of the single SNP methods to use when only 1 SNP was used to
  estimate the causal effect? The default is `"Wald ratio"`.

- multi_snp_method:

  Which of the multi-SNP methods to use when there was more than 1 SNPs
  used to estimate the causal effect? The default is
  `"Inverse variance weighted"`.

- ao_slc:

  Logical; retrieve sample size and subcategory using
  [`available_outcomes()`](https://mrcieu.github.io/TwoSampleMR/reference/available_outcomes.md).
  If set to `FALSE` `mr_res` must contain the following additional
  columns: `subcategory` and `sample_size`.

- priority:

  Name of category to prioritise at the top of the forest plot. The
  default is `"Cardiometabolic"`.

## Value

data frame.

## Details

By default it uses the
[`available_outcomes()`](https://mrcieu.github.io/TwoSampleMR/reference/available_outcomes.md)
function to retrieve the study level characteristics for the outcome
trait, including sample size and outcome category. This assumes the MR
analysis was performed using outcome GWAS(s) contained in OpenGWAS.

If `ao_slc` is set to `TRUE` then the user must supply their own study
level characteristics. This is useful when the user has supplied their
own outcome GWAS results (i.e. they are not in OpenGWAS).
