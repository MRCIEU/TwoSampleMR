# Find instruments for use in MR from the OpenGWAS database

This function searches for GWAS significant SNPs (for a given p-value)
for a specified set of outcomes. It then performs LD based clumping to
return only independent significant associations.

## Usage

``` r
extract_instruments(
  outcomes,
  p1 = 5e-08,
  clump = TRUE,
  p2 = 5e-08,
  r2 = 0.001,
  kb = 10000,
  opengwas_jwt = ieugwasr::get_opengwas_jwt(),
  force_server = FALSE
)
```

## Arguments

- outcomes:

  Array of outcome IDs (see
  [`available_outcomes()`](https://mrcieu.github.io/TwoSampleMR/reference/available_outcomes.md)).

- p1:

  Significance threshold. The default is `5e-8`.

- clump:

  Logical; whether to clump results. The default is `TRUE`.

- p2:

  Secondary clumping threshold. The default is `5e-8`.

- r2:

  Clumping r2 cut off. The default is `0.001`.

- kb:

  Clumping distance cutoff. The default is `10000`.

- opengwas_jwt:

  Used to authenticate protected endpoints. Login to
  <https://api.opengwas.io> to obtain a jwt. Provide the jwt string
  here, or store in .Renviron under the keyname OPENGWAS_JWT.

- force_server:

  Force the analysis to extract results from the server rather than the
  MRInstruments package.

## Value

data frame
