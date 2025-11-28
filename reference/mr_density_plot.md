# Density plot

Density plot

## Usage

``` r
mr_density_plot(
  singlesnp_results,
  mr_results,
  exponentiate = FALSE,
  bandwidth = "nrd0"
)
```

## Arguments

- singlesnp_results:

  from
  [`mr_singlesnp()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_singlesnp.md).

- mr_results:

  Results from
  [`mr()`](https://mrcieu.github.io/TwoSampleMR/reference/mr.md).

- exponentiate:

  Plot on exponentiated scale. The default is `FALSE`.

- bandwidth:

  Density bandwidth parameter.

## Value

List of plots
