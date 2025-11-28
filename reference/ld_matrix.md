# Get LD matrix for list of SNPs

This function takes a list of SNPs and searches for them in a specified
super-population in the 1000 Genomes phase 3 reference panel. It then
creates an LD matrix of r values (signed, and not squared). All LD
values are with respect to the major alleles in the 1000G dataset. You
can specify whether the allele names are displayed.

## Usage

``` r
ld_matrix(snps, with_alleles = TRUE, pop = "EUR")
```

## Arguments

- snps:

  List of SNPs.

- with_alleles:

  Whether to append the allele names to the SNP names. The default is
  `TRUE`.

- pop:

  Super-population to use as reference panel. Default = `"EUR"`. Options
  are `"EUR"`, `"SAS"`, `"EAS"`, `"AFR"`, `"AMR"`. `'legacy'` also
  available - which is a previously used version of the EUR panel with a
  slightly different set of markers.

## Value

Matrix of LD r values

## Details

The data used for generating the LD matrix includes only bi-allelic SNPs
with MAF \> 0.01, so it's quite possible that a variant you want to
include will be absent. If it is absent, it will be automatically
excluded from the results.

You can check if your variants are present in the LD reference panel
using
[`ieugwasr::ld_reflookup()`](https://mrcieu.github.io/ieugwasr/reference/ld_reflookup.html).

This function does put load on the OpenGWAS servers, which makes life
more difficult for other users, and has been limited to analyse only up
to 500 variants at a time. We have implemented a method and made
available the LD reference panels to perform the operation locally, see
[`ieugwasr::ld_matrix()`](https://mrcieu.github.io/ieugwasr/reference/ld_matrix.html)
and related vignettes for details.
