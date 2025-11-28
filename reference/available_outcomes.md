# Get list of studies with available GWAS summary statistics through API

Get list of studies with available GWAS summary statistics through API

## Usage

``` r
available_outcomes(opengwas_jwt = ieugwasr::get_opengwas_jwt())
```

## Arguments

- opengwas_jwt:

  Used to authenticate protected endpoints. Login to
  <https://api.opengwas.io> to obtain a jwt. Provide the jwt string
  here, or store in .Renviron under the keyname OPENGWAS_JWT.

## Value

Dataframe of details for all available studies
