# Convert TwoSampleMR format to MendelianRandomization format

The MendelianRandomization package offers MR methods that can be used
with the same data used in the TwoSampleMR package. This function
converts from the TwoSampleMR format to the MRInput class.

## Usage

``` r
dat_to_MRInput(dat, get_correlations = FALSE, pop = "EUR")
```

## Arguments

- dat:

  Output from the
  [`harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.md)
  function.

- get_correlations:

  Default `FALSE`. If `TRUE` then extract the LD matrix for the SNPs
  from the European 1000 genomes data on OpenGWAS.

- pop:

  If `get_correlations` is `TRUE` then use the following

## Value

List of MRInput objects for each exposure/outcome combination
