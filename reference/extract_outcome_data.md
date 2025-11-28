# Supply the output from [`read_exposure_data()`](https://mrcieu.github.io/TwoSampleMR/reference/read_exposure_data.md) and all the SNPs therein will be queried against the requested outcomes in remote database using API.

Supply the output from
[`read_exposure_data()`](https://mrcieu.github.io/TwoSampleMR/reference/read_exposure_data.md)
and all the SNPs therein will be queried against the requested outcomes
in remote database using API.

## Usage

``` r
extract_outcome_data(
  snps,
  outcomes,
  proxies = TRUE,
  rsq = 0.8,
  align_alleles = 1,
  palindromes = 1,
  maf_threshold = 0.3,
  opengwas_jwt = ieugwasr::get_opengwas_jwt(),
  splitsize = 10000,
  proxy_splitsize = 500
)
```

## Arguments

- snps:

  Array of SNP rs IDs.

- outcomes:

  Array of IDs (see `id` column in output from
  [`available_outcomes()`](https://mrcieu.github.io/TwoSampleMR/reference/available_outcomes.md)).

- proxies:

  Look for LD tags? Default is `TRUE`.

- rsq:

  Minimum LD rsq value (if proxies = 1). Default = `0.8`.

- align_alleles:

  Try to align tag alleles to target alleles (if proxies = 1). `1` =
  yes, `0` = no. The default is `1`.

- palindromes:

  Allow palindromic SNPs (if proxies = 1). `1` = yes, `0` = no. The
  default is `1`.

- maf_threshold:

  MAF threshold to try to infer palindromic SNPs. The default is `0.3`.

- opengwas_jwt:

  Used to authenticate protected endpoints. Login to
  <https://api.opengwas.io> to obtain a jwt. Provide the jwt string
  here, or store in .Renviron under the keyname OPENGWAS_JWT.

- splitsize:

  The default is `10000`.

- proxy_splitsize:

  The default is `500`.

## Value

Dataframe of summary statistics for all available outcomes
