# Estimate trait SD by obtaining beta estimates from z-scores and finding the ratio with original beta values

Assumes that sample size and allele frequency is correct for each SNP,
and that allele frequency gives a reasonable estimate of the variance of
the SNP.

## Usage

``` r
estimate_trait_sd(b, se, n, p)
```

## Arguments

- b:

  vector of effect sizes.

- se:

  vector of standard errors.

- n:

  vector of sample sizes.

- p:

  vector of allele frequencies.

## Value

Vector of sd estimates for each association.
