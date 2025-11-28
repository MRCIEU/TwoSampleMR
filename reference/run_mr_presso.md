# Wrapper for MR-PRESSO

See <https://github.com/rondolab/MR-PRESSO> for more details.

## Usage

``` r
run_mr_presso(dat, NbDistribution = 1000, SignifThreshold = 0.05)
```

## Arguments

- dat:

  Output from
  [`harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.md).

- NbDistribution:

  Number of bootstrap replications. The default is `1000`.

- SignifThreshold:

  Outlier significance threshold. The default is `0.05`.

## Value

List of results for every exposure/outcome combination
