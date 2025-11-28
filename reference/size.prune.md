# Size prune

Whens there are duplicate summary sets for a particular exposure-outcome
combination, this function drops the duplicates with the smaller total
sample size (for binary outcomes, the number of cases is used instead of
total sample size).

## Usage

``` r
size.prune(dat)
```

## Arguments

- dat:

  Results from
  [`harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.md).

## Value

data frame
