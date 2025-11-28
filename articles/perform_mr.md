# Perform MR

``` r
library(TwoSampleMR)
library(ggplot2)
```

## Introduction

Let’s continue with the example of body mass index (BMI) on coronary
heart disease (CHD):

``` r
bmi_exp_dat <- extract_instruments(outcomes = 'ieu-a-2')
chd_out_dat <- extract_outcome_data(snps = bmi_exp_dat$SNP, outcomes = 'ieu-a-7')
```

``` r
dat <- harmonise_data(bmi_exp_dat, chd_out_dat)
#> Harmonising Body mass index || id:ieu-a-2 (ieu-a-2) and Coronary heart disease || id:ieu-a-7 (ieu-a-7)
```

Once the exposure and outcome data are harmonised, we have effects and
standard errors for each instrument SNP available for the exposure and
outcome traits. We can use this information to perform a Mendelian
randomization (MR) analysis. To do this, simply run:

``` r
res <- mr(dat)
#> Analysing 'ieu-a-2' on 'ieu-a-7'
res
#>   id.exposure id.outcome                              outcome
#> 1     ieu-a-2    ieu-a-7 Coronary heart disease || id:ieu-a-7
#> 2     ieu-a-2    ieu-a-7 Coronary heart disease || id:ieu-a-7
#> 3     ieu-a-2    ieu-a-7 Coronary heart disease || id:ieu-a-7
#> 4     ieu-a-2    ieu-a-7 Coronary heart disease || id:ieu-a-7
#> 5     ieu-a-2    ieu-a-7 Coronary heart disease || id:ieu-a-7
#>                        exposure                    method nsnp         b
#> 1 Body mass index || id:ieu-a-2                  MR Egger   79 0.5024935
#> 2 Body mass index || id:ieu-a-2           Weighted median   79 0.3870065
#> 3 Body mass index || id:ieu-a-2 Inverse variance weighted   79 0.4459091
#> 4 Body mass index || id:ieu-a-2               Simple mode   79 0.3401554
#> 5 Body mass index || id:ieu-a-2             Weighted mode   79 0.3790910
#>           se         pval
#> 1 0.14396056 8.012590e-04
#> 2 0.07639889 4.071091e-07
#> 3 0.05898302 4.032020e-14
#> 4 0.15001315 2.612742e-02
#> 5 0.10221374 3.881092e-04
```

This returns a data frame of estimates of the causal effect of the
exposure on the outcome for a range of different MR methods.

If there were multiple exposures against multiple outcomes in `dat`, the
[`mr()`](https://mrcieu.github.io/TwoSampleMR/reference/mr.md) function
will perform each MR method for each combination of exposure-outcome
traits.

## MR methods

The list of available MR methods can be obtained:

``` r
mr_method_list()
#>                              obj
#> 1                  mr_wald_ratio
#> 2               mr_two_sample_ml
#> 3            mr_egger_regression
#> 4  mr_egger_regression_bootstrap
#> 5               mr_simple_median
#> 6             mr_weighted_median
#> 7   mr_penalised_weighted_median
#> 8                         mr_ivw
#> 9                  mr_ivw_radial
#> 10                    mr_ivw_mre
#> 11                     mr_ivw_fe
#> 12                mr_simple_mode
#> 13              mr_weighted_mode
#> 14         mr_weighted_mode_nome
#> 15           mr_simple_mode_nome
#> 16                       mr_raps
#> 17                       mr_sign
#> 18                        mr_uwr
#> 19                       mr_grip
#>                                                         name PubmedID
#> 1                                                 Wald ratio         
#> 2                                         Maximum likelihood         
#> 3                                                   MR Egger 26050253
#> 4                                       MR Egger (bootstrap) 26050253
#> 5                                              Simple median 27061298
#> 6                                            Weighted median 27061298
#> 7                                  Penalised weighted median 27061298
#> 8                                  Inverse variance weighted 24114802
#> 9                                                 IVW radial 29961852
#> 10 Inverse variance weighted (multiplicative random effects)         
#> 11                 Inverse variance weighted (fixed effects)         
#> 12                                               Simple mode 29040600
#> 13                                             Weighted mode 29040600
#> 14                                      Weighted mode (NOME) 29040600
#> 15                                        Simple mode (NOME) 29040600
#> 16                      Robust adjusted profile score (RAPS)         
#> 17                                     Sign concordance test         
#> 18                                     Unweighted regression         
#> 19                                                   MR GRIP         
#>                                                    Description use_by_default
#> 1                                                                        TRUE
#> 2                                                                       FALSE
#> 3                                                                        TRUE
#> 4                                                                       FALSE
#> 5                                                                       FALSE
#> 6                                                                        TRUE
#> 7                                                                       FALSE
#> 8                                                                        TRUE
#> 9                                                                       FALSE
#> 10                                                                      FALSE
#> 11                                                                      FALSE
#> 12                                                                       TRUE
#> 13                                                                       TRUE
#> 14                                                                      FALSE
#> 15                                                                      FALSE
#> 16                                                                      FALSE
#> 17 Tests for concordance of signs between exposure and outcome          FALSE
#> 18                                     Doesn't use any weights          FALSE
#> 19                 Allele coding invariant MR-Egger regression          FALSE
#>    heterogeneity_test
#> 1               FALSE
#> 2                TRUE
#> 3                TRUE
#> 4               FALSE
#> 5               FALSE
#> 6               FALSE
#> 7               FALSE
#> 8                TRUE
#> 9                TRUE
#> 10              FALSE
#> 11              FALSE
#> 12              FALSE
#> 13              FALSE
#> 14              FALSE
#> 15              FALSE
#> 16              FALSE
#> 17              FALSE
#> 18               TRUE
#> 19              FALSE
```

To perform them, they can be specified in the
[`mr()`](https://mrcieu.github.io/TwoSampleMR/reference/mr.md) function,
e.g. to only perform MR Egger regression and Inverse variance weighted
methods,

``` r
mr(dat, method_list = c("mr_egger_regression", "mr_ivw"))
#> Analysing 'ieu-a-2' on 'ieu-a-7'
#>   id.exposure id.outcome                              outcome
#> 1     ieu-a-2    ieu-a-7 Coronary heart disease || id:ieu-a-7
#> 2     ieu-a-2    ieu-a-7 Coronary heart disease || id:ieu-a-7
#>                        exposure                    method nsnp         b
#> 1 Body mass index || id:ieu-a-2                  MR Egger   79 0.5024935
#> 2 Body mass index || id:ieu-a-2 Inverse variance weighted   79 0.4459091
#>           se        pval
#> 1 0.14396056 8.01259e-04
#> 2 0.05898302 4.03202e-14
```

By default, all the methods that are labelled `TRUE` in the
`use_by_default` column are used by the
[`mr()`](https://mrcieu.github.io/TwoSampleMR/reference/mr.md) function.

## Sensitivity analyses

### Heterogeneity statistics

Some of the MR methods can also perform tests for heterogeneity. To
obtain those statistics:

``` r
mr_heterogeneity(dat)
#>   id.exposure id.outcome                              outcome
#> 1     ieu-a-2    ieu-a-7 Coronary heart disease || id:ieu-a-7
#> 2     ieu-a-2    ieu-a-7 Coronary heart disease || id:ieu-a-7
#>                        exposure                    method        Q Q_df
#> 1 Body mass index || id:ieu-a-2                  MR Egger 143.3046   77
#> 2 Body mass index || id:ieu-a-2 Inverse variance weighted 143.6508   78
#>         Q_pval
#> 1 6.841585e-06
#> 2 8.728420e-06
```

As with the
[`mr()`](https://mrcieu.github.io/TwoSampleMR/reference/mr.md) function,
the
[`mr_heterogeneity()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_heterogeneity.md)
function can take an argument to only perform heterogeneity tests using
specified methods, e.g.

``` r
mr_heterogeneity(dat, method_list = c("mr_egger_regression", "mr_ivw"))
#>   id.exposure id.outcome                              outcome
#> 1     ieu-a-2    ieu-a-7 Coronary heart disease || id:ieu-a-7
#> 2     ieu-a-2    ieu-a-7 Coronary heart disease || id:ieu-a-7
#>                        exposure                    method        Q Q_df
#> 1 Body mass index || id:ieu-a-2                  MR Egger 143.3046   77
#> 2 Body mass index || id:ieu-a-2 Inverse variance weighted 143.6508   78
#>         Q_pval
#> 1 6.841585e-06
#> 2 8.728420e-06
```

### Horizontal pleiotropy

The intercept term in MR Egger regression can be a useful indication of
whether directional horizontal pleiotropy is driving the results of an
MR analysis. This can be obtained as follows:

``` r
mr_pleiotropy_test(dat)
#>   id.exposure id.outcome                              outcome
#> 1     ieu-a-2    ieu-a-7 Coronary heart disease || id:ieu-a-7
#>                        exposure egger_intercept          se      pval
#> 1 Body mass index || id:ieu-a-2    -0.001719304 0.003985962 0.6674266
```

### Single SNP analysis

To obtain the MR estimates using each of the SNPs singly we can do the
following:

``` r
res_single <- mr_singlesnp(dat)
```

This returns a data.frame of results that is similar to the output from
[`mr()`](https://mrcieu.github.io/TwoSampleMR/reference/mr.md) except it
performs the analysis multiple times for each exposure-outcome
combination - each time using a different single SNP to perform the
analysis.

The method used to perform the single SNP MR is the Wald ratio by
default, though this can be changed, e.g. to use the fixed effects meta
analysis method instead:

``` r
res_single <- mr_singlesnp(dat, single_method = "mr_meta_fixed")
```

The
[`mr_singlesnp()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_singlesnp.md)
function calculates the full MR using all available SNPs as well, and by
default it uses the IVW and MR Egger methods. This can be specified as
so:

``` r
res_single <- mr_singlesnp(dat, all_method = "mr_two_sample_ml")
```

will perform only the maximum likelihood method for the combined test.

### Leave-one-out analysis

It is possible to perform a leave-one-out analysis, where the MR is
performed again but leaving out each SNP in turn, to identify if a
single SNP is driving the association.

``` r
res_loo <- mr_leaveoneout(dat)
head(res_loo)
#>                        exposure                              outcome
#> 1 Body mass index || id:ieu-a-2 Coronary heart disease || id:ieu-a-7
#> 2 Body mass index || id:ieu-a-2 Coronary heart disease || id:ieu-a-7
#> 3 Body mass index || id:ieu-a-2 Coronary heart disease || id:ieu-a-7
#> 4 Body mass index || id:ieu-a-2 Coronary heart disease || id:ieu-a-7
#> 5 Body mass index || id:ieu-a-2 Coronary heart disease || id:ieu-a-7
#> 6 Body mass index || id:ieu-a-2 Coronary heart disease || id:ieu-a-7
#>   id.exposure id.outcome samplesize        SNP         b         se
#> 1     ieu-a-2    ieu-a-7     184305  rs1000940 0.4541284 0.05843628
#> 2     ieu-a-2    ieu-a-7     184305 rs10132280 0.4450525 0.05959723
#> 3     ieu-a-2    ieu-a-7     184305  rs1016287 0.4443855 0.05959844
#> 4     ieu-a-2    ieu-a-7     184305 rs10182181 0.4427817 0.05994598
#> 5     ieu-a-2    ieu-a-7     184305 rs10733682 0.4474244 0.05955045
#> 6     ieu-a-2    ieu-a-7     184305 rs10840100 0.4435657 0.05956120
#>              p
#> 1 7.765841e-15
#> 2 8.162633e-14
#> 3 8.896652e-14
#> 4 1.509205e-13
#> 5 5.762610e-14
#> 6 9.532322e-14
```

By default the method used is the inverse variance weighted method, but
this can be changed by using the `method` argument.

## Plots

There are a few ways to visualise the results, listed below

### Scatter plot

We can depict the relationship of the SNP effects on the exposure
against the SNP effects on the outcome using a scatter plot.

``` r
res <- mr(dat)
#> Analysing 'ieu-a-2' on 'ieu-a-7'
p1 <- mr_scatter_plot(res, dat)
```

A scatter plot is created for each exposure-outcome test, and stored in
`p1` as a list of plots. For example, to plot the first scatter plot:

``` r
p1[[1]]
```

![A scatter plot visualising the two-sample data points and the
following fitted models; Inverse Variance Weighted, MR-Egger, Simple
mode, Weighted median, and Weighted
mode.](perform_mr_files/figure-html/unnamed-chunk-18-1.png)

And to see how many plots there are:

``` r
length(p1)
#> [1] 1
```

Lines are drawn for each method used in `mr(dat)`, the slope of the line
corresponding to the estimated causal effect. To limit which lines are
drawn, simply specify the desired methods, e.g., to only draw the lines
for the MR Egger and IVW causal effect estimates:

``` r
res <- mr(dat, method_list = c("mr_egger_regression", "mr_ivw"))
#> Analysing 'ieu-a-2' on 'ieu-a-7'
p1 <- mr_scatter_plot(res, dat)
p1[[1]]
```

![A scatter plot visualising the two-sample data points and the
following fitted models; Inverse Variance Weighted and
MR-Egger.](perform_mr_files/figure-html/unnamed-chunk-20-1.png)

It is possible to save this plot using the
[`ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
function from the `ggplot2` package, e.g. to save as a pdf

``` r
ggsave(p1[[1]], file = "filename.pdf", width = 7, height = 7)
```

Or to save as a png.

``` r
ggsave(p1[[1]], file = "filename.png", width = 7, height = 7)
```

See `?ggplot2::ggsave()` for more info.

### Forest plot

Use the
[`mr_forest_plot()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_forest_plot.md)
function to compare the MR estimates using the different MR methods
against the single SNP tests.

``` r
p2 <- mr_forest_plot(res_single)
p2[[1]]
```

![A forest plot showing the estimated causal effects using each SNP
separately, and the Inverse Variance Weighted and MR-Egger estimates
using all the
SNPs.](perform_mr_files/figure-html/unnamed-chunk-23-1.png)

Here, the plot shows the causal effect as estimated using each of the
SNPs individually, and comparing against the causal effect as estimated
using the methods that use all the SNPs.

To obtain plots that use different methods, specify them in the
[`mr_singlesnp()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_singlesnp.md)
function:

``` r
res_single <- mr_singlesnp(dat, all_method = c("mr_ivw", "mr_two_sample_ml"))
p2 <- mr_forest_plot(res_single)
p2[[1]]
```

![An alternative forest plot showing the estimated causal effects using
each SNP separately, and the Inverse Variance Weighted and Maximum
Likelihood estimates using all the
SNPs.](perform_mr_files/figure-html/unnamed-chunk-24-1.png)

### Leave-one-out plot

Use the
[`mr_leaveoneout_plot()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_leaveoneout_plot.md)
function to visualise the leave-one-out analysis:

``` r
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]
```

![A leave one out plot showing the Inverse Variance Weighted estimate
with each SNP
omitted.](perform_mr_files/figure-html/unnamed-chunk-25-1.png)

We can also specify the method to use in the function call, e.g., to use
MR-Egger regression in the leave-one-out analysis specify

``` r
mr_leaveoneout(dat, method = mr_egger_regression)
```

### Funnel plot

Asymmetry in a funnel plot is useful for gauging the reliability of a
particular MR analysis. Funnel plots can be produced using the single
SNP results as follows:

``` r
p4 <- mr_funnel_plot(res_single)
p4[[1]]
```

![A funnel plot showing the causal effect for each SNP and the inverse
variance weighted and MR-Egger estimates using all the
SNPs.](perform_mr_files/figure-html/unnamed-chunk-26-1.png)

## 1-to-many forest plot

A 1-to-many MR analysis interrogates the effect of a single exposure on
multiple outcomes or multiple exposures on a single outcome. The results
of this analysis can be visualised using the 1-to-many forest plot, with
or without stratification on a categorical variable. From a visual point
of view, the function works best for 50 or fewer results and is not
really designed to handle more than a 100 results. If your number of
results is much greater than 50, it may be better to split these across
two separate plots. For example, if you have 100 sets of results you
could divide these equally across two plots and then combine the two
plots together in another program like Powerpoint. The function assumes
the results are already in the right order for plotting. As such, users
are advised to sort their results according to how they would like them
to appear in the plot. Users can use their own code to do this or they
can use the
[`sort_1_to_many()`](https://mrcieu.github.io/TwoSampleMR/reference/sort_1_to_many.md)
function.

### Step 1: generate 1-to-many MR results

``` r
exp_dat <- extract_instruments(outcomes = c("ieu-a-2", "ieu-a-100", "ieu-a-1032", "ieu-a-104", "ieu-a-1", "ieu-a-72", "ieu-a-999"))
table(exp_dat$exposure)
chd_out_dat <- extract_outcome_data(
  snps = exp_dat$SNP,
  outcomes = "ieu-a-7"
)

dat2 <- harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = chd_out_dat
)
res <- mr(dat2)
```

### Step 2. Make the 1-to-many forest plot

#### Example 1. Effect of multiple risk factors on coronary heart disease

In this example we wish to plot results from an MR analysis of the
effect of multiple exposures on coronary heart disease, with results
sorted by decreasing effect size (largest effect at the top of the plot)
and with one MR method for each unique exposure-outcome combination. We
will also make the size of each point estimate proportional to its
inverse variance. This is a useful way to draw attention towards the
most reliable results and away from results with very wide confidence
intervals. To specify the size of the point estimate, set the weight
argument to the name of the column in the data with the weight
information.

``` r
# default is to subset on either the IVW method (>1 instrumental SNP) or Wald ratio method (1 instrumental SNP).
res <- subset_on_method(res)
# this sorts results by decreasing effect size (largest effect at top of the plot)
res <- sort_1_to_many(res, b = "b", sort_action = 4)
# to keep the Y axis label clean we exclude the exposure ID labels from the exposure column
res <- split_exposure(res)
res$weight <- 1 / res$se
```

``` r
# identify value for 'lo' in forest_plot_1_to_many
min(exp(res$b - 1.96 * res$se))
#> [1] 0.3233985
# identify value for 'up' in forest_plot_1_to_many
max(exp(res$b + 1.96 * res$se))
#> [1] 2.489726

forest_plot_1_to_many(
  res,
  b = "b",
  se = "se",
  exponentiate = TRUE,
  ao_slc = FALSE,
  lo = 0.3,
  up = 2.5,
  TraitM = "exposure",
  col1_width = 2,
  by = NULL,
  trans = "log2",
  xlab = "OR for CHD per SD increase in risk factor (95% confidence interval)",
  weight = "weight"
)
#> Warning: Removed 6 rows containing missing values or values outside the scale range
#> (`geom_vline()`).
```

![A forest plot showing the estimated causal effects for each risk
factor on the
outcome.](perform_mr_files/figure-html/unnamed-chunk-30-1.png)

It is also possible to add additional columns and column titles and to
choose the size of the text in the columns:

``` r
res$pval <- formatC(res$pval, format = "e", digits = 2)

forest_plot_1_to_many(
  res,
  b = "b",
  se = "se",
  exponentiate = TRUE,
  ao_slc = FALSE,
  lo = 0.3,
  up = 2.5,
  TraitM = "exposure",
  by = NULL,
  trans = "log2",
  xlab = "OR for CHD per SD increase in risk factor (95% CI)",
  weight = "weight",
  subheading_size = 11,
  col1_title = "Risk factor",
  col1_width = 2.5,
  col_text_size = 4,
  addcols = c("nsnp", "pval"),
  addcol_widths = c(1.0, 1.0),
  addcol_titles = c("No. SNPs", "P-val")
)
#> [1] "nsnp"
#> [1] "pval"
#> Warning: Removed 6 rows containing missing values or values outside the scale range
#> (`geom_vline()`).
```

![A forest plot showing the estimated causal effects for each risk
factor on the outcome including columns for the no. SNPs and
p-value.](perform_mr_files/figure-html/unnamed-chunk-31-1.png)

In my own workflow I prefer to to keep the plot free of axis and column
titles and to add them separately in a program like Powerpoint:

``` r
forest_plot_1_to_many(
  res,
  b = "b",
  se = "se",
  exponentiate = TRUE,
  ao_slc = FALSE,
  lo = 0.3,
  up = 3.0,
  TraitM = "exposure",
  col1_width = 2.0,
  by = NULL,
  trans = "log2",
  xlab = "",
  addcols = c("nsnp", "pval"),
  weight = "weight",
  col_text_size = 4,
  addcol_widths = c(0.5, 1.0),
  addcol_titles = c("", "")
)
#> [1] "nsnp"
#> [1] "pval"
```

![A forest plot showing the estimated causal effects for each risk
factor on the outcome including columns for the no. SNPs and p-value
without column
headers.](perform_mr_files/figure-html/unnamed-chunk-32-1.png)

#### Example 2. MR results for multiple MR methods grouped by multiple exposures

In this next example we plot the results from an analysis of the effect
of multiple exposures on coronary heart disease using multiple methods,
with results grouped by exposure. We also want the result for the IVW
method to be given priority and to go above the other methods. We also
want the exposure with the largest IVW effect size to go the top of the
plot. We also set the `TraitM` argument to the column describing the MR
method. This is because we are grouping the results on the exposures.
Normally the row labels would correspond to the exposures but in this
example we want the row names to correspond to the MR method.

``` r
res <- mr(dat2)
#> Analysing 'ieu-a-1' on 'ieu-a-7'
#> Analysing 'ieu-a-100' on 'ieu-a-7'
#> Analysing 'ieu-a-104' on 'ieu-a-7'
#> Analysing 'ieu-a-2' on 'ieu-a-7'
#> Analysing 'ieu-a-72' on 'ieu-a-7'
#> Analysing 'ieu-a-999' on 'ieu-a-7'
# to keep the Y axis label clean we exclude the exposure ID labels from the exposure column
res <- split_exposure(res)

res <-
  sort_1_to_many(
    res,
    group = "exposure",
    sort_action = 3,
    priority = "Inverse variance weighted",
    trait_m = "method"
  )
```

``` r
forest_plot_1_to_many(
  res,
  b = "b",
  se = "se",
  exponentiate = TRUE,
  trans = "log2",
  ao_slc = FALSE,
  lo = 0.03,
  up = 22,
  col1_width = 2,
  by = "exposure",
  TraitM = "method",
  xlab = "OR for CHD per SD increase in risk factor (95% confidence interval)",
  subheading_size = 12,
  col_text_size = 4
)
```

![A forest plot showing the estimated causal effects for each method for
each risk factor on the
outcome.](perform_mr_files/figure-html/unnamed-chunk-34-1.png)

#### Example 3. Stratify results on a grouping variable

In this next example we plot the same results as above but with results
stratified by a grouping variable. We also select one MR method for each
unique exposure-outcome combination and sort the results by decreasing
effect size within each group (i.e. largest effect at the top).

``` r
res <- mr(dat2)
#> Analysing 'ieu-a-1' on 'ieu-a-7'
#> Analysing 'ieu-a-100' on 'ieu-a-7'
#> Analysing 'ieu-a-104' on 'ieu-a-7'
#> Analysing 'ieu-a-2' on 'ieu-a-7'
#> Analysing 'ieu-a-72' on 'ieu-a-7'
#> Analysing 'ieu-a-999' on 'ieu-a-7'
res <- split_exposure(res)
res <- subset_on_method(res)
res$subcategory[res$exposure %in% c("Adiponectin", "Hip circumference", "Waist circumference")] <- "Group 1"
res$subcategory[is.na(res$subcategory)] <- "Group 2"
res$weight <- 1 / res$se
res <- sort_1_to_many(res, sort_action = 1, group = "subcategory")

forest_plot_1_to_many(
  res,
  b = "b",
  se = "se",
  exponentiate = TRUE,
  trans = "log2",
  ao_slc = FALSE,
  lo = 0.3,
  up = 2.5,
  TraitM = "exposure",
  col_text_size = 4,
  col1_width = 1.5,
  by = "subcategory",
  xlab = "OR for CHD per SD increase in risk factor (95% confidence interval)",
  subheading_size = 14,
  weight = "weight"
)
#> Warning: Removed 3 rows containing missing values or values outside the scale range
#> (`geom_vline()`).
#> Removed 3 rows containing missing values or values outside the scale range
#> (`geom_vline()`).
```

![Forest plot split by a subcategory
variable.](perform_mr_files/figure-html/unnamed-chunk-35-1.png)

In the above example we made up an arbitrary grouping variable called
“subcategory” with values “Group 1” and “Group 2”. Typically, however,
the grouping variable might correspond to something like a trait
ontology (e.g. anthropometric and glycemic traits) or study design
(e.g. MR and observational studies).

#### Example 4. Effect of BMI on 103 diseases

The plot function works best with 50 or fewer rows and is not really
designed to handle more than a 100. Visualising a single-column forest
plot with 100 results is also quite difficult. If your number of results
is much greater than 50, it is advisable to split the results across two
different plots. In the example below we select BMI as the exposure and
test this against 103 diseases in the IEU GWAS database:

``` r
# extract instruments for BMI
exp_dat <- extract_instruments(outcomes = "ieu-a-2")
ao <- available_outcomes()
# identify diseases
ao <- ao[ao$category == "Disease", ]
ao <- ao[which(ao$ncase > 100), ]

dis_dat <- extract_outcome_data(
  snps = exp_dat$SNP,
  outcomes = ao$id
)

dat3 <- harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = dis_dat
)

res <- mr(dat3, method_list = c("mr_wald_ratio", "mr_ivw"))
# to keep the Y axis label clean we exclude the exposure ID labels from the exposure column
res <- split_outcome(res)
# this sorts results by decreasing effect size (largest effect at top of the plot)
res <- sort_1_to_many(res, b = "b", sort_action = 4)
```

MR results for 103 diseases can be difficult to visualise in a
single-column forest plot. In my own workflow I would split these across
two plots and then join them together in a separate program, such as
Powerpoint, and do further refinements there. I typically save my plots
using the [`pdf()`](https://rdrr.io/r/grDevices/pdf.html) graphics
device. In this particular example the disease labels probably require
some cleaning up (some are a bit long) or alternatively the column text
size could be made smaller. It is also possible to change the colour of
the plot and the shape of the point estimates. Type
[`?forest_plot_1_to_many`](https://mrcieu.github.io/TwoSampleMR/reference/forest_plot_1_to_many.md)
for further details.

``` r
res1 <- res[1:52, ]
res2 <- res[53:103, ]

plot1 <- forest_plot_1_to_many(
  res1,
  b = "b",
  se = "se",
  exponentiate = TRUE,
  trans = "log2",
  ao_slc = FALSE,
  lo = 0.004,
  up = 461,
  col1_width = 2,
  TraitM = "outcome",
  col_text_size = 3,
  xlab = ""
)

plot2 <- forest_plot_1_to_many(
  res2,
  b = "b",
  se = "se",
  exponentiate = TRUE,
  trans = "log2",
  ao_slc = FALSE,
  lo = 0.004,
  up = 461,
  col1_width = 2,
  TraitM = "outcome",
  subheading_size = 11,
  col_text_size = 3,
  xlab = ""
)

plot1
plot2
```

``` r
pdf("plot1.pdf", height = 10, width = 8)
plot1
dev.off()
```

## MR-RAPS: Many weak instruments analysis

MR-RAPS (Robust Adjusted Profile Score) is a recently proposed method
that considers the measurement error in SNP-exposure effects, is
unbiased when there are many (e.g. hundreds of) weak instruments, and is
robust to systematic and idiosyncratic pleiotropy. See Zhao et al.
(2020) for more detail about the statistical methodology.

MR-RAPS is implemented in the R package **mr.raps** that is available on
[GitHub](https://github.com/qingyuanzhao/mr.raps) which is installed
when you install TwoSampleMR. It can be directly called from TwoSampleMR
by

    #> Harmonising Body mass index || id:ieu-a-2 (ieu-a-2) and Coronary heart disease || id:ieu-a-7 (ieu-a-7)

``` r
res <- mr(dat, method_list = c("mr_raps"))
#> Analysing 'ieu-a-2' on 'ieu-a-7'
res
#>   id.exposure id.outcome                              outcome
#> 1     ieu-a-2    ieu-a-7 Coronary heart disease || id:ieu-a-7
#>                        exposure                               method nsnp
#> 1 Body mass index || id:ieu-a-2 Robust adjusted profile score (RAPS)   79
#>           b        se         pval
#> 1 0.4646893 0.0542478 1.070969e-17
```

MR-RAPS comes with two main options: `over.dispersion` (whether the
method should consider systematic pleiotropy) and `loss.function`
(either `"l2"`, `"huber"`, or `"tukey"`). The latter two loss functions
are robust to idiosyncratic pleiotropy. The default option is
`over.dispersion = TRUE`, `loss.function = "huber"`, and
`shrinkage = FALSE`. To change these options, modify the `parameters`
argument of
[`mr()`](https://mrcieu.github.io/TwoSampleMR/reference/mr.md) by, for
example,

``` r
res <-
  mr(
    dat,
    method_list = c("mr_raps"),
    parameters = list(over.dispersion = FALSE, loss.function = "l2", shrinkage = TRUE)
  )
#> Analysing 'ieu-a-2' on 'ieu-a-7'
res
#>   id.exposure id.outcome                              outcome
#> 1     ieu-a-2    ieu-a-7 Coronary heart disease || id:ieu-a-7
#>                        exposure                               method nsnp
#> 1 Body mass index || id:ieu-a-2 Robust adjusted profile score (RAPS)   79
#>           b         se        pval
#> 1 0.4580697 0.04515624 3.51986e-24
```

## MR-GRIP

MR-GRIP is a recently proposed method due to Dudbridge and Bowden et
al. (2025) which is a modification to the MR-Egger method which has the
Genotype Recoding Invariance Property (GRIP), it can be implemented
using the code below.

``` r
res <- mr(dat, method_list = c("mr_grip"))
#> Analysing 'ieu-a-2' on 'ieu-a-7'
res
#>   id.exposure id.outcome                              outcome
#> 1     ieu-a-2    ieu-a-7 Coronary heart disease || id:ieu-a-7
#>                        exposure  method nsnp         b        se         pval
#> 1 Body mass index || id:ieu-a-2 MR GRIP   79 0.4899371 0.0946796 1.765179e-06
```

Or include it with the default set of methods as follows.

``` r
mr(dat, method_list = c(subset(mr_method_list(), use_by_default)$obj, "mr_grip"))
```

## Reports

A report can be generated that performs all MR analyses, sensitivity
analyses, and plots, and presents them in a single self-contained html
web page, word document, or pdf document.

``` r
mr_report(dat)
```

By default this produces a html file in the current working directory,
but see the help pages on how to modify this.

This function will create a separate report file for every
exposure-outcome combination that is present in the `dat` object.

## MR Steiger directionality test

This is an implementation of the method described here:

[Hemani G, Tilling K, Davey Smith G. **Orienting the causal relationship
between imprecisely measured traits using GWAS summary data.** PLoS
Genetics. 2017. 13(11):
e1007081.](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007081)

In MR it is assumed that the instruments influence the exposure first
and then the outcome through the exposure. But sometimes this is
difficult to evaluate, for example is a cis-acting SNP influencing gene
expression levels or DNA methylation levels first? The causal direction
between the hypothesised exposure and outcomes can be tested using the
Steiger test (Hemani, Tilling, and Davey Smith 2017). For example:

``` r
out <- directionality_test(dat)
#> r.exposure and/or r.outcome not present.
#> Calculating approximate SNP-exposure and/or SNP-outcome correlations, assuming all are quantitative traits. Please pre-calculate r.exposure and/or r.outcome using get_r_from_lor() for any binary traits
knitr::kable(out)
```

| id.exposure | id.outcome | exposure                        | outcome                                | snp_r2.exposure | snp_r2.outcome | correct_causal_direction | steiger_pval |
|:------------|:-----------|:--------------------------------|:---------------------------------------|----------------:|---------------:|:-------------------------|-------------:|
| ieu-a-2     | ieu-a-7    | Body mass index \|\| id:ieu-a-2 | Coronary heart disease \|\| id:ieu-a-7 |       0.0158082 |      0.0013505 | TRUE                     |            0 |

It calculates the variance explained in the exposure and the outcome by
the instrumenting SNPs, and tests if the variance in the outcome is less
than the exposure.

This test is, like many others, liable to give inaccurate causal
directions under some measurement error parameters in the exposure and
the outcome (e.g. if the outcome has much lower measurement precision
then its proportion of variance explained will be underestimated).
Sensitivity can be applied to evaluate the extent to which the inferred
causal direction is liable to measurement error, in two ways.

1.  Provide estimates of measurement error for the exposure and the
    outcome, and obtain an adjusted estimate of the causal direction
2.  For all possible values of measurement error, identify the
    proportion of the parameter space which supports the inferred causal
    direction

These tests are obtained using:

``` r
mr_steiger(
  p_exp = dat$pval.exposure,
  p_out = dat$pval.outcome,
  n_exp = dat$samplesize.exposure,
  n_out = dat$samplesize.outcome,
  r_exp = rep(NA, length(dat$pval.exposure)),
  r_out = rep(NA, length(dat$pval.exposure))
)
```

## Multivariable MR

When SNPs instrument multiple potential exposures, for example in the
case of different lipid fractions, one method for overcoming this
problem is to estimate the influence of each lipid conditioning on the
effects of the SNPs on the other lipids. Multivariable MR can be
performed using the R package as follows. Practically speaking, this is
the process that needs to occur from the perspective of generating the
data in the correct format:

1.  Get instruments for each exposure
2.  Combine these into a set of all instruments
3.  Clump these to avoid the possibility that e.g. a variant for
    exposure 1 is in LD with a variant for exposure 2
4.  Re-extract all the final clumped SNPs from (3) from all of the
    exposures
5.  Harmonise them all to be on the same effect allele
6.  Use the multivariable MR method against these harmonised data

Example - The GWAS IDs for HDL, LDL and total cholesterol are
`ieu-a-299`, `ieu-a-300`, and `ieu-a-302`. The GWAS ID for coronary
heart disease (CHD) is `ieu-a-7`. In this example we will estimate the
multivariable effects of HDL, LDL and total cholesterol on CHD.

``` r
id_exposure <- c("ieu-a-299", "ieu-a-300", "ieu-a-302")
id_outcome <- "ieu-a-7"
```

First obtain the instruments for each lipid fraction. This entails
obtaining a combined set of SNPs including all instruments, and getting
those SNPs for each lipid fraction. Therefore, if there are e.g. 20
instruments for each of 3 lipid fractions, but combined there are 30
unique SNPs, then we need to extract each of the 30 SNPs from each lipid
fraction (exposure).

``` r
mv_exposure_dat <- mv_extract_exposures(id_exposure)
```

Next, also extract those SNPs from the outcome.

``` r
mv_outcome_dat <- extract_outcome_data(exposure_dat$SNP, id_outcome)
```

Once the data has been obtained, harmonise so that all are on the same
reference allele.

``` r
mvdat <- mv_harmonise_data(mv_exposure_dat, mv_outcome_dat)
#> Harmonising HDL cholesterol || id:ieu-a-299 (ieu-a-299) and Coronary heart disease || id:ieu-a-7 (ieu-a-7)
```

Finally, perform the multivariable MR analysis, which generates a table
of results.

``` r
res <- mv_multiple(mvdat)
res
#> $result
#>   id.exposure                        exposure id.outcome
#> 1   ieu-a-299 HDL cholesterol || id:ieu-a-299    ieu-a-7
#> 2   ieu-a-300 LDL cholesterol || id:ieu-a-300    ieu-a-7
#> 3   ieu-a-302   Triglycerides || id:ieu-a-302    ieu-a-7
#>                                outcome nsnp           b         se         pval
#> 1 Coronary heart disease || id:ieu-a-7   79 -0.08919724 0.05970552 1.351879e-01
#> 2 Coronary heart disease || id:ieu-a-7   68  0.37853543 0.04976846 2.828614e-14
#> 3 Coronary heart disease || id:ieu-a-7   42  0.13584165 0.06738291 4.380354e-02
```

### Note about multivariable methods

There are several different ways in which this analysis can be
formulated, e.g., consider 3 exposures against one outcome, one could:

1.  Fit all exposures together or fit one exposure at a time against the
    residuals of the outcome that has been adjusted for the other
    outcomes. The former is recommended by default in this R package
    through the
    [`mv_multiple()`](https://mrcieu.github.io/TwoSampleMR/reference/mv_multiple.md)
    function but the latter was how MV MR was originally described by
    Burgess and Thompson (2015) and can be done with
    [`mv_residual()`](https://mrcieu.github.io/TwoSampleMR/reference/mv_residual.md);
2.  Fitting all instruments for all exposures (default) or only fitting
    the instruments for each exposure sequentially;
3.  Forcing the slopes through the origin (default) or allowing an
    intercept term.

With these three different parameters there are eight different ways to
do multivariable (MV) analysis. We recommend the default settings as
described above.

### Note about visualisation

Plots can be generated using the `plots = TRUE` argument for
[`mv_multiple()`](https://mrcieu.github.io/TwoSampleMR/reference/mv_multiple.md)
and
[`mv_residual()`](https://mrcieu.github.io/TwoSampleMR/reference/mv_residual.md).

The current plots being generated are not necessarily adequate because
while they show the slope through the raw points, they do not
demonstrate that the raw points might be effectively different between
plots because they are conditional on the other exposures.

### Using your own summary data

If you want to perform analysis with your local summary data (i.e. not
in the OpenGWAS database) then use then look up the
[`mv_extract_exposures_local()`](https://mrcieu.github.io/TwoSampleMR/reference/mv_extract_exposures_local.md)
function instead of the
[`mv_extract_exposures()`](https://mrcieu.github.io/TwoSampleMR/reference/mv_extract_exposures.md)
function.

## MR estimates when instruments are correlated

In the examples shown so far it is assumed that instruments are
independent (i.e. they are not in linkage disequilibrium, LD). This is
to avoid ‘double counting’ effects. An alternative approach is to
estimate the MR effects accounting for the correlation between variants.

The TwoSampleMR package has not implemented this yet, but the
[MendelianRandomization](https://CRAN.R-project.org/package=MendelianRandomization)
R package by Yavorska and Burgess does have this functionality. We can
use the TwoSampleMR package to extract, format and harmonise data, and
then convert to the format required by the MendelianRandomization
package. The IEU GWAS database server has the individual level genetic
data for ~500 Europeans in 1000 genomes data, and can obtain the LD
matrix for a set of SNPs using these data. For example:

``` r
snplist <- c("rs234", "rs1205")
ld_mat <- ld_matrix(snplist)
```

``` r
ld_mat
#>            rs234_A_G rs1205_T_C
#> rs234_A_G  1.0000000  0.0797023
#> rs1205_T_C 0.0797023  1.0000000
```

Here
[`ld_matrix()`](https://mrcieu.github.io/TwoSampleMR/reference/ld_matrix.md)
returns the LD correlation values (not R²) for each pair of variants
present in the 1000 genomes data set.

``` r
dat <- harmonise_data(
  exposure_dat = bmi_exp_dat,
  outcome_dat = chd_out_dat
)
#> Harmonising Body mass index || id:ieu-a-2 (ieu-a-2) and Coronary heart disease || id:ieu-a-7 (ieu-a-7)
```

Convert to the `MRInput` format for the MendelianRandomization package:

``` r
dat2 <- dat_to_MRInput(dat)
#> Converting:
#>  - exposure: Body mass index || id:ieu-a-2
#>  - outcome: Coronary heart disease || id:ieu-a-7
```

This produces a list of `MRInput` objects that can be used with the
MendelianRandomization functions, e.g.

``` r
MendelianRandomization::mr_ivw(dat2[[1]])
#> 
#> Inverse-variance weighted method
#> (variants uncorrelated, random-effect model)
#> 
#> Number of Variants : 79 
#> 
#> ------------------------------------------------------------------
#>  Method Estimate Std Error 95% CI       p-value
#>     IVW    0.446     0.059 0.330, 0.562   0.000
#> ------------------------------------------------------------------
#> Residual standard error =  1.357 
#> Heterogeneity test statistic (Cochran's Q) = 143.6508 on 78 degrees of freedom, (p-value = 0.0000). I^2 = 45.7%. 
#> F statistic = 65.6.
```

Alternatively, convert to the `MRInput` format but also obtaining the LD
matrix for the instruments

``` r
dat2 <- try(dat_to_MRInput(dat, get_correlation = TRUE))
#> Converting:
#>  - exposure: Body mass index || id:ieu-a-2
#>  - outcome: Coronary heart disease || id:ieu-a-7
#>  - obtaining LD matrix
#> Please look at vignettes for options on running this locally if you need to run many instances of this command.
#> Error in get_query_content(.) : 
#> Status code from OpenGWAS API: 401
#> 
#> Message: ERROR - Go to https://api.opengwas.io/ - From 1st May 2024 you must provide a token (JWT) alongside most of your requests. Read more at https://api.opengwas.io/ and also check for the latest version at https://mrcieu.github.io/ieugwasr/
if (class(dat2) != "try-error") MendelianRandomization::mr_ivw(dat2[[1]], correl = TRUE)
```

## MR-MoE: Using a mixture of experts machine learning approach

We recently developed MR-MoE, a method to choose the most appropriate
amongst several MR tests using a machine learning algorithm. Note that
the method is still under review, but full details are described here:
<https://doi.org/10.1101/173682>.

MR-MoE operates by taking a set of harmonised data, inferring some
characteristics about the dataset, and using those characteristics to
predict how well each of the different MR methods will perform on the
dataset, in terms of maximising power while minimising false discovery
rates.

In order to run the analysis you must download an RData object that
contains the trained random forests that are used to predict the
efficacy of each method. This can be downloaded from here:

[dropbox.com/s/5la7y38od95swcf](https://mrcieu.github.io/TwoSampleMR/articles/dropbox.com/s/5la7y38od95swcf)

**Caution: this is a large file (approx 167Mb)**

Once downloaded, read in the object and use the
[`mr_moe()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_moe.md)
function to perform the analysis. An example is shown here, estimating
the causal effect of BMI on coronary heart disease:

``` r
# Extact instruments for BMI
exposure_dat <- extract_instruments("ieu-a-2")

# Get corresponding effects for CHD
outcome_dat <- extract_outcome_data(exposure_dat$SNP, "ieu-a-7")

# Harmonise
dat <- harmonise_data(exposure_dat, outcome_dat)

# Load the downloaded RData object. This loads the rf object
load("rf.rdata")

# Obtain estimates from all methods, and generate data metrics
res_all <- mr_wrapper(dat)

# MR-MoE - predict the performance of each method
res_moe <- mr_moe(res_all, rf)
```

The function does the following:

1.  Performs MR using each of 11 MR methods
2.  Applies Steiger filtering or heterogeneity filtering or both to
    remove SNPs that do not have substantially larger R² with the
    exposure than the outcome. Note - for binary traits ensure number of
    cases, number of controls, and allele frequencies are available for
    each SNP. For continuous traits make sure the p-value and sample
    size is available. The function infers if a trait is binary or
    continuous based on the units.exposure and units.outcome columns -
    binary traits must have those values set to ‘log odds’
3.  Performs the 14 MR methods again but using the subset of SNPs that
    survive Steiger filtering
4.  Generates meta data about the summary data to predict the most
    reliable of the 28 methods applied.

For every exposure / outcome combination in the `dat` object, the MR-MoE
method is applied. The function returns a list which is as long as the
number of exposure / outcome combinations. In this case, it will be of
length 1, containing the result for BMI on CHD.

The result object itself is a list with the following elements:

- `estimates` (results from each MR)
- `heterogeneity` (results from heterogeneity for different filtering
  approaches)
- `directional_pleiotropy` (egger intercepts)
- `info` (metrics used to generate MOE)

Looking at the `estimates`, we see that there is a column called `MOE`
which is the predicted AUROC curve performance of each method.

## Post MR results management

The TwoSampleMR package also provides the following functions for
managing or editing MR results.

### Split outcome names

The outcome column in the output of mr() combines the original outcome
name with the outcome trait ID.

``` r
head(res)
#> $result
#>   id.exposure                        exposure id.outcome
#> 1   ieu-a-299 HDL cholesterol || id:ieu-a-299    ieu-a-7
#> 2   ieu-a-300 LDL cholesterol || id:ieu-a-300    ieu-a-7
#> 3   ieu-a-302   Triglycerides || id:ieu-a-302    ieu-a-7
#>                                outcome nsnp           b         se         pval
#> 1 Coronary heart disease || id:ieu-a-7   79 -0.08919724 0.05970552 1.351879e-01
#> 2 Coronary heart disease || id:ieu-a-7   68  0.37853543 0.04976846 2.828614e-14
#> 3 Coronary heart disease || id:ieu-a-7   42  0.13584165 0.06738291 4.380354e-02
```

The outcome column can be split into separate columns for the id and
outcome name using the split_outcome function:

``` r
res <- mr(dat)
#> Analysing 'ieu-a-2' on 'ieu-a-7'
split_outcome(res)
#>   id.exposure id.outcome                outcome                      exposure
#> 1     ieu-a-2    ieu-a-7 Coronary heart disease Body mass index || id:ieu-a-2
#> 2     ieu-a-2    ieu-a-7 Coronary heart disease Body mass index || id:ieu-a-2
#> 3     ieu-a-2    ieu-a-7 Coronary heart disease Body mass index || id:ieu-a-2
#> 4     ieu-a-2    ieu-a-7 Coronary heart disease Body mass index || id:ieu-a-2
#> 5     ieu-a-2    ieu-a-7 Coronary heart disease Body mass index || id:ieu-a-2
#>                      method nsnp         b         se         pval
#> 1                  MR Egger   79 0.5024935 0.14396056 8.012590e-04
#> 2           Weighted median   79 0.3870065 0.07671416 4.540309e-07
#> 3 Inverse variance weighted   79 0.4459091 0.05898302 4.032020e-14
#> 4               Simple mode   79 0.3401554 0.15196883 2.805152e-02
#> 5             Weighted mode   79 0.3790910 0.10242200 3.980416e-04
```

### Split exposure names

Similarly to the outcome column, the exposure column in the output of
[`mr()`](https://mrcieu.github.io/TwoSampleMR/reference/mr.md) combines
the original exposure name with the exposure trait ID. This can be split
into separate columns for the id and exposure name using the
split_exposure function.

### Generate odds ratios with 95% confidence intervals

Users can convert log odds ratios into odds ratios with 95% confidence
intervals using:

``` r
generate_odds_ratios(res)
#>   id.exposure id.outcome                              outcome
#> 1     ieu-a-2    ieu-a-7 Coronary heart disease || id:ieu-a-7
#> 2     ieu-a-2    ieu-a-7 Coronary heart disease || id:ieu-a-7
#> 3     ieu-a-2    ieu-a-7 Coronary heart disease || id:ieu-a-7
#> 4     ieu-a-2    ieu-a-7 Coronary heart disease || id:ieu-a-7
#> 5     ieu-a-2    ieu-a-7 Coronary heart disease || id:ieu-a-7
#>                        exposure                    method nsnp         b
#> 1 Body mass index || id:ieu-a-2                  MR Egger   79 0.5024935
#> 2 Body mass index || id:ieu-a-2           Weighted median   79 0.3870065
#> 3 Body mass index || id:ieu-a-2 Inverse variance weighted   79 0.4459091
#> 4 Body mass index || id:ieu-a-2               Simple mode   79 0.3401554
#> 5 Body mass index || id:ieu-a-2             Weighted mode   79 0.3790910
#>           se         pval      lo_ci     up_ci       or or_lci95 or_uci95
#> 1 0.14396056 8.012590e-04 0.22033081 0.7846562 1.652838 1.246489 2.191653
#> 2 0.07671416 4.540309e-07 0.23664672 0.5373662 1.472566 1.266993 1.711493
#> 3 0.05898302 4.032020e-14 0.33030238 0.5615158 1.561909 1.391389 1.753328
#> 4 0.15196883 2.805152e-02 0.04229652 0.6380143 1.405166 1.043204 1.892719
#> 5 0.10242200 3.980416e-04 0.17834387 0.5798381 1.460956 1.195236 1.785749
```

### Subset on method

It is sometimes useful to subset results on MR method, so that there is
one unique result for each exposure-outcome combination:

``` r
subset_on_method(res)
#>   id.exposure id.outcome                              outcome
#> 3     ieu-a-2    ieu-a-7 Coronary heart disease || id:ieu-a-7
#>                        exposure                    method nsnp         b
#> 3 Body mass index || id:ieu-a-2 Inverse variance weighted   79 0.4459091
#>           se        pval
#> 3 0.05898302 4.03202e-14
```

The default is to subset on the IVW method when \>1 SNP is available and
to use the Wald ratio method when a single SNP is available. Users can
specify which multi-SNP method to subset on.

### Combine all results

It is often useful to combine all results and study level
characteristics into a single dataframe or table, e.g. for sharing
results with collaborators or when the user wishes to present all
results in a single table or figure. This can be done using the
[`combine_all_mrresults()`](https://mrcieu.github.io/TwoSampleMR/reference/combine_all_mrresults.md)
function:

``` r
res <- mr(dat)
het <- mr_heterogeneity(dat)
plt <- mr_pleiotropy_test(dat)
sin <- mr_singlesnp(dat)
all_res <-
  combine_all_mrresults(
    res,
    het,
    plt,
    sin,
    ao_slc = TRUE,
    Exp = TRUE,
    split.exposure = FALSE,
    split.outcome = TRUE
  )
head(all_res[, c(
  "Method",
  "outcome",
  "exposure",
  "nsnp",
  "b",
  "se",
  "pval",
  "intercept",
  "intercept_se",
  "intercept_pval",
  "Q",
  "Q_df",
  "Q_pval",
  "consortium",
  "ncase",
  "ncontrol",
  "pmid",
  "population"
)])
```

This combines all results from
[`mr()`](https://mrcieu.github.io/TwoSampleMR/reference/mr.md),
[`mr_heterogeneity()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_heterogeneity.md),
[`mr_pleiotropy_test()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_pleiotropy_test.md)
and
[`mr_singlesnp()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_singlesnp.md)
into a single dataframe. It also merges the results with outcome study
level characteristics from the
[`available_outcomes()`](https://mrcieu.github.io/TwoSampleMR/reference/available_outcomes.md)
function, including sample size characteristics. If requested, it also
exponentiates results (e.g. if the user wants log odds ratio converted
into odds ratios with 95 percent confidence intervals).

## References

Burgess, S, and S Thompson. 2015. “Multivariable Mendelian
Randomization: The Use of Pleiotropic Genetic Variants to Estimate
Causal Effects.” *American Journal of Epidemiology* 181: 251–60.
<https://doi.org/10.1093/aje/kwu283>.

Hemani, Gibran, Kate Tilling, and George Davey Smith. 2017. “Orienting
the Causal Relationship Between Imprecisely Measured Traits Using GWAS
Summary Data.” *PLOS Genetics* 13 (11): e1007081.
<https://doi.org/10.1371/journal.pgen.1007081>.

Zhao, Qingyuan, Jingshu Wang, Gibran Hemani, Jack Bowden, and Dylan S
Small. 2020. “Statistical inference in two-sample summary-data Mendelian
randomization using robust adjusted profile score” 48: 1742–69.
<https://doi.org/10.1214/19-AOS1866>.
