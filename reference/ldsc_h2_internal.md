# Univariate LDSC

Imported here to help estimate sample overlap between studies

## Usage

``` r
ldsc_h2_internal(Z, r2, N, W = NULL)
```

## Arguments

- Z:

  summary Z-statistics for M variants

- r2:

  average reference LD scores for M variants

- N:

  GWAS sample size for each variant (could be different across variants)

- W:

  variant weight

## Value

model fit
