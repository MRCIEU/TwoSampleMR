# Forest plot for multiple exposures and multiple outcomes

Perform MR of multiple exposures and multiple outcomes. This plots the
results.

## Usage

``` r
forest_plot(
  mr_res,
  exponentiate = FALSE,
  single_snp_method = "Wald ratio",
  multi_snp_method = "Inverse variance weighted",
  group_single_categories = TRUE,
  by_category = TRUE,
  in_columns = FALSE,
  threshold = NULL,
  xlab = "",
  xlim = NULL,
  trans = "identity",
  ao_slc = TRUE,
  priority = "Cardiometabolic"
)
```

## Arguments

- mr_res:

  Results from
  [`mr()`](https://mrcieu.github.io/TwoSampleMR/reference/mr.md).

- exponentiate:

  Convert effects to OR? Default is `FALSE`.

- single_snp_method:

  Which of the single SNP methods to use when only 1 SNP was used to
  estimate the causal effect? The default is `"Wald ratio"`.

- multi_snp_method:

  Which of the multi-SNP methods to use when there was more than 1 SNPs
  used to estimate the causal effect? The default is
  `"Inverse variance weighted"`.

- group_single_categories:

  If there are categories with only one outcome, group them together
  into an "Other" group. The default is `TRUE`.

- by_category:

  Separate the results into sections by category? The default is `TRUE`.

- in_columns:

  Separate the exposures into different columns. The default is `FALSE`.

- threshold:

  p-value threshold to use for colouring points by significance level.
  If `NULL` (default) then colour layer won't be applied.

- xlab:

  x-axis label. If `in_columns=TRUE` then the exposure values are
  appended to the end of `xlab`. e.g. if `xlab="Effect of"` then
  x-labels will read `"Effect of exposure1"`, `"Effect of exposure2"`
  etc. Otherwise will be printed as is.

- xlim:

  limit x-axis range. Provide vector of length 2, with lower and upper
  bounds. The default is `NULL`.

- trans:

  Transformation to apply to x-axis. e.g. `"identity"`, `"log2"`, etc.
  The default is `"identity"`.

- ao_slc:

  retrieve sample size and subcategory from
  [`available_outcomes()`](https://mrcieu.github.io/TwoSampleMR/reference/available_outcomes.md).
  If set to `FALSE` then `mr_res` must contain the following additional
  columns: `sample_size` and `subcategory`. The default behaviour is to
  use
  [`available_outcomes()`](https://mrcieu.github.io/TwoSampleMR/reference/available_outcomes.md)
  to retrieve sample size and subcategory.

- priority:

  Name of category to prioritise at the top of the forest plot. The
  default is `"Cardiometabolic"`.

## Value

grid plot object
