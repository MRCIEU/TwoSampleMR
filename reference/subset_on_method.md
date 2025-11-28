# Subset MR-results on method

This function takes MR results from
[`mr()`](https://mrcieu.github.io/TwoSampleMR/reference/mr.md) and
restricts to a single method per exposure x disease combination.

## Usage

``` r
subset_on_method(
  mr_res,
  single_snp_method = "Wald ratio",
  multi_snp_method = "Inverse variance weighted"
)
```

## Arguments

- mr_res:

  Results from
  [`mr()`](https://mrcieu.github.io/TwoSampleMR/reference/mr.md).

- single_snp_method:

  Which of the single SNP methods to use when only 1 SNP was used to
  estimate the causal effect? The default is `"Wald ratio"`.

- multi_snp_method:

  Which of the multi-SNP methods to use when there was more than 1 SNPs
  used to estimate the causal effect? The default is
  `"Inverse variance weighted"`.

## Value

data frame.
