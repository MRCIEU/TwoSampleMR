# Bivariate LDSC

Imported here to help estimate sample overlap between studies

## Usage

``` r
ldsc_rg_internal(Zs, r2, h1, h2, N1, N2, Nc = 0, W = NULL)
```

## Arguments

- Zs:

  Mx2 matrix of summary Z-statistics for M variants from two GWAS

- r2:

  average reference LD scores for M variants

- h1:

  hsq for trait 1

- h2:

  hsq for trait 2

- N1:

  sample size for the 1st GWAS

- N2:

  sample size for the 2nd GWAS

- Nc:

  overlapped sample size between the two GWAS

- W:

  variant weight

## Value

List of models

## References

Bulik-Sullivan,B.K. et al. (2015) An atlas of genetic correlations
across human diseases and traits. Nat. Genet. 47, 1236–1241.

Guo,B. and Wu,B. (2018) Principal component based adaptive association
test of multiple traits using GWAS summary statistics. bioRxiv 269597;
doi: 10.1101/269597

Gua,B. and Wu,B. (2019) Integrate multiple traits to detect novel
trait-gene association using GWAS summary data with an adaptive test
approach. Bioinformatics. 2019 Jul 1;35(13):2251-2257. doi:
10.1093/bioinformatics/bty961.

https://github.com/baolinwu/MTAR
