# MR mode estimators

Perform simple, weighted, penalised modes, as well as versions that use
the NOME assumption.

## Usage

``` r
mr_mode(dat, parameters = default_parameters(), mode_method = "all")
```

## Arguments

- dat:

  Output from
  [`harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.md).

- parameters:

  List of parameters. The default is
  [`default_parameters()`](https://mrcieu.github.io/TwoSampleMR/reference/default_parameters.md).

- mode_method:

  The default is `"all"`. The other choices are `'Simple mode'`,
  `'Weighted mode'`, `'Penalised mode'`, `'Simple mode (NOME)'`,
  `'Weighted mode (NOME)'`.

## Value

data frame
