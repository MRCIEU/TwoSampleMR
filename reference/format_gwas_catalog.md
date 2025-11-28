# Get data selected from GWAS catalog into correct format

DEPRECATED. Please use
[`format_data()`](https://mrcieu.github.io/TwoSampleMR/reference/format_data.md)
instead.

## Usage

``` r
format_gwas_catalog(gwas_catalog_subset, type = "exposure")
```

## Arguments

- gwas_catalog_subset:

  The GWAS catalog subset.

- type:

  The default is `"exposure"`.

## Value

Data frame

## Examples

``` r
if (FALSE) { # \dontrun{
require(MRInstruments)
data(gwas_catalog)
bmi <- subset(gwas_catalog, Phenotype=="Body mass index" & Year==2010 & grepl("kg", Units))
bmi <- format_data(bmi)
} # }
```
