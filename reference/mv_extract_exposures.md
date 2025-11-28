# Extract exposure variables for multivariable MR

Requires a list of IDs from available_outcomes. For each ID, it extracts
instruments. Then, it gets the full list of all instruments and extracts
those SNPs for every exposure. Finally, it keeps only the SNPs that are
a) independent and b) present in all exposures, and harmonises them to
be all on the same strand.

## Usage

``` r
mv_extract_exposures(
  id_exposure,
  clump_r2 = 0.001,
  clump_kb = 10000,
  harmonise_strictness = 2,
  opengwas_jwt = ieugwasr::get_opengwas_jwt(),
  find_proxies = TRUE,
  force_server = FALSE,
  pval_threshold = 5e-08,
  pop = "EUR",
  plink_bin = NULL,
  bfile = NULL
)
```

## Arguments

- id_exposure:

  Array of IDs (e.g. c(299, 300, 302) for HDL, LDL, trigs)

- clump_r2:

  The default is `0.01`.

- clump_kb:

  The default is `10000`.

- harmonise_strictness:

  See the `action` option of
  [`harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.md).
  The default is `2`.

- opengwas_jwt:

  Used to authenticate protected endpoints. Login to
  <https://api.opengwas.io> to obtain a jwt. Provide the jwt string
  here, or store in .Renviron under the keyname OPENGWAS_JWT.

- find_proxies:

  Look for proxies? This slows everything down but is more accurate. The
  default is `TRUE`.

- force_server:

  Whether to search through pre-clumped dataset or to re-extract and
  clump directly from the server. The default is `FALSE`.

- pval_threshold:

  Instrument detection p-value threshold. Default = `5e-8`

- pop:

  Which 1000 genomes super population to use for clumping when using the
  server

- plink_bin:

  If `NULL` and `bfile` is not `NULL` then will detect packaged plink
  binary for specific OS. Otherwise specify path to plink binary.
  Default = `NULL`

- bfile:

  If this is provided then will use the API. Default = `NULL`

## Value

data frame in `exposure_dat` format
