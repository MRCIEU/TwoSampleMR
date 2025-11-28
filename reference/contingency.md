# Obtain 2x2 contingency table from marginal parameters and odds ratio

Columns are the case and control frequencies. Rows are the frequencies
for allele 1 and allele 2.

## Usage

``` r
contingency(af, prop, odds_ratio, eps = 1e-15)
```

## Arguments

- af:

  Allele frequency of effect allele.

- prop:

  Proportion of cases.

- odds_ratio:

  Odds ratio.

- eps:

  tolerance, default is `1e-15`.

## Value

2x2 contingency table as matrix
