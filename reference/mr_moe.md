# Mixture of experts

Based on the method described here
<https://www.biorxiv.org/content/10.1101/173682v2>. Once all MR methods
have been applied to a summary set, you can then use the mixture of
experts to predict the method most likely to be the most accurate.

## Usage

``` r
mr_moe(res, rf)
```

## Arguments

- res:

  Output from
  [`mr_wrapper()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_wrapper.md).

- rf:

  The trained random forest for the methods. This is available to
  download at <https://www.dropbox.com/s/5la7y38od95swcf/rf.rdata?dl=0>.

## Value

List

## Details

The `mr_moe()` function modifies the `estimates` item in the list of
results from the
[`mr_wrapper()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_wrapper.md)
function. It does three things:

1.  Adds the MOE column, which is a predictor for each method for how
    well it performs in terms of high power and low type 1 error (scaled
    0-1, where 1 is best performance).

2.  It renames the methods to be the estimating method + the instrument
    selection method. There are 4 instrument selection methods: Tophits
    (i.e. no filtering), directional filtering (DF, an unthresholded
    version of Steiger filtering), heterogeneity filtering (HF, removing
    instruments that make substantial (p \< 0.05) contributions to
    Cochran's Q statistic), and DF + HF which is where DF is applied and
    the HF applied on top of that.

3.  It orders the table to be in order of best performing method.

Note that the mixture of experts has only been trained on datasets with
at least 5 SNPs. If your dataset has fewer than 5 SNPs this function
might return errors.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example of body mass index on coronary heart disease
# Extract and harmonise data
a <- extract_instruments("ieu-a-2")
b <- extract_outcome_data(a$SNP, 7)
dat <- harmonise_data(a, b)

# Apply all MR methods
r <- mr_wrapper(dat)

# Load the rf object containing the trained models
load("rf.rdata")
# Update the results with mixture of experts
r <- mr_moe(r, rf)

# Now you can view the estimates, and see that they have
# been sorted in order from most likely to least likely to
# be accurate, based on MOE prediction
r[[1]]$estimates
} # }
```
