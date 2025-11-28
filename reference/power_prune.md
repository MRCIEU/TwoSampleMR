# Power prune

When there are duplicate summary sets for a particular exposure-outcome
combination, this function keeps the exposure-outcome summary set with
the highest expected statistical power. This can be done by dropping the
duplicate summary sets with the smaller sample sizes. Alternatively, the
pruning procedure can take into account instrument strength and outcome
sample size. The latter is useful, for example, when there is
considerable variation in SNP coverage between duplicate summary sets
(e.g. because some studies have used targeted or fine mapping arrays).
If there are a large number of SNPs available to instrument an exposure,
the outcome GWAS with the better SNP coverage may provide better power
than the outcome GWAS with the larger sample size.

## Usage

``` r
power_prune(dat, method = 1, dist.outcome = "binary")
```

## Arguments

- dat:

  Results from
  [`harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.md).

- method:

  Should the duplicate summary sets be pruned on the basis of sample
  size alone (`method = 1`) or a combination of instrument strength and
  sample size (`method = 2`)? Default set to `1`. When set to 1, the
  duplicate summary sets are first dropped on the basis of the outcome
  sample size (smaller duplicates dropped). If duplicates are still
  present, remaining duplicates are dropped on the basis of the exposure
  sample size (smaller duplicates dropped). When method is set to `2`,
  duplicates are dropped on the basis of instrument strength (amount of
  variation explained in the exposure by the instrumental SNPs) and
  sample size, and assumes that the SNP-exposure effects correspond to a
  continuous trait with a normal distribution (i.e. exposure cannot be
  binary). The SNP-outcome effects can correspond to either a binary or
  continuous trait. If the exposure is binary then `method=1` should be
  used.

- dist.outcome:

  The distribution of the outcome. Can either be `"binary"` or
  `"continuous"`. Default set to `"binary"`.

## Value

data.frame with duplicate summary sets removed
