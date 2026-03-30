# Forest plot

If the data frame contains a `category` column, SNPs will be coloured
and grouped by category in the forest plot (e.g., cluster assignments
from MR-Clust). See Figure 4 of Vabistsevits et al. (2024) for examples.

## Usage

``` r
mr_forest_plot(
  singlesnp_results,
  exponentiate = FALSE,
  category_colours = NULL
)
```

## Arguments

- singlesnp_results:

  from
  [`mr_singlesnp()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_singlesnp.md).

- exponentiate:

  Plot on exponential scale. The default is `FALSE`.

- category_colours:

  Named character vector of colours for categories. Names should match
  the values in the `category` column. If `NULL` (the default), the
  Okabe-Ito colourblind-friendly palette is used.

## Value

List of plots

## References

Vabistsevits, M., Davey Smith, G., Richardson, T.G. et al. Mammographic
density mediates the protective effect of early-life body size on breast
cancer risk. *Nature Communications*, **15**, 4021 (2024).
[doi:10.1038/s41467-024-48105-7](https://doi.org/10.1038/s41467-024-48105-7)

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic forest plot
bmi_exp_dat <- extract_instruments(outcomes = "ieu-a-2")
chd_out_dat <- extract_outcome_data(
  snps = bmi_exp_dat$SNP, outcomes = "ieu-a-7"
)
dat <- harmonise_data(bmi_exp_dat, chd_out_dat)
res <- mr_singlesnp(dat)
mr_forest_plot(res)

# Forest plot with RadialMR outliers
radial_dat <- dat_to_RadialMR(dat)
radial_res <- RadialMR::ivw_radial(radial_dat[[1]], alpha = 0.05, weights = 3)
outlier_snps <- radial_res$outliers$SNP
snp_rows <- !grepl("^All", res$SNP)
res$category <- NA_character_
res$category[snp_rows] <- ifelse(
  res$SNP[snp_rows] %in% outlier_snps, "Outlier", "Main"
)
mr_forest_plot(res)

# With custom colours
mr_forest_plot(res, category_colours = c(Main = "grey50", Outlier = "red"))
} # }
```
