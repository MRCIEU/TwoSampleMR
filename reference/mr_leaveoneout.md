# Leave one out sensitivity analysis

Leave one out sensitivity analysis

## Usage

``` r
mr_leaveoneout(dat, parameters = default_parameters(), method = mr_ivw)
```

## Arguments

- dat:

  Output from
  [`harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.md).

- parameters:

  List of parameters.

- method:

  Choose which method to use. The default is `mr_ivw`.

## Value

List of data frames
