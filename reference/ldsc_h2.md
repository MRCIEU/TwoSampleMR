# Univariate LDSC

Imported here to help estimate sample overlap between studies

## Usage

``` r
ldsc_h2(id, ancestry = "infer", snpinfo = NULL, splitsize = 20000)
```

## Arguments

- id:

  ID to analyse

- ancestry:

  ancestry of traits 1 and 2 (AFR, AMR, EAS, EUR, SAS) or 'infer'
  (default) in which case it will try to guess based on allele
  frequencies

- snpinfo:

  Output from `ieugwasr::afl2_list("hapmap3")`, or `NULL` for it to be
  done automatically

- splitsize:

  How many SNPs to extract at one time. Default=`20000`

## Value

model fit

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

<https://github.com/baolinwu/MTAR>
