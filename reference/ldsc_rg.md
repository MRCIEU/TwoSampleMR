# Bivariate LDSC

Imported here to help estimate sample overlap between studies

## Usage

``` r
ldsc_rg(id1, id2, ancestry = "infer", snpinfo = NULL, splitsize = 20000)
```

## Arguments

- id1:

  ID 1 to analyse

- id2:

  ID 2 to analyse

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
