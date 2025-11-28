# Convenient function to create a harmonised dataset

Convenient function to create a harmonised dataset.

## Usage

``` r
make_dat(
  exposures = c("ieu-a-2", "ieu-a-301"),
  outcomes = c("ieu-a-7", "ieu-a-1001"),
  proxies = TRUE
)
```

## Arguments

- exposures:

  The default is `c("ieu-a-2", "ieu-a-301")` (BMI and LDL).

- outcomes:

  The default is `c("ieu-a-7", "ieu-a-1001")` (CHD and EDU).

- proxies:

  Look for proxies? Default = `TRUE`

## Value

Harmonised data frame
