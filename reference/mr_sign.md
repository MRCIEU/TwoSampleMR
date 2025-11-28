# MR sign test

Tests how often the SNP-exposure and SNP-outcome signs are concordant.
This is to avoid the problem of averaging over all SNPs, which can
suffer bias due to outliers with strong effects; and to avoid excluding
SNPs which is implicit in median and mode based estimators. The effect
estimate here is not to be interpreted as the effect size - it is the
proportion of SNP-exposure and SNP-outcome effects that have concordant
signs. e.g. +1 means all have the same sign, -1 means all have opposite
signs, and 0 means that there is an equal number of concordant and
discordant signs. Restricted to only work if there are 6 or more valid
SNPs.

## Usage

``` r
mr_sign(b_exp, b_out, se_exp = NULL, se_out = NULL, parameters = NULL)
```

## Arguments

- b_exp:

  Vector of genetic effects on exposure

- b_out:

  Vector of genetic effects on outcome

- se_exp:

  Not required

- se_out:

  Not required

- parameters:

  Not required

## Value

List with the following elements:

- b:

  Concordance (see description)

- se:

  NA

- pval:

  p-value

- nsnp:

  Number of SNPs (excludes NAs and effect estimates that are 0)
