# Calculate variance explained from p-values and sample size

This method is an approximation, and may be numerically unstable.
Ideally you should estimate r directly from independent replication
samples. Use
[`get_r_from_lor()`](https://mrcieu.github.io/TwoSampleMR/reference/get_r_from_lor.md)
for binary traits.

## Usage

``` r
get_r_from_pn(p, n)
```

## Arguments

- p:

  Array of pvals

- n:

  Array of sample sizes

## Value

Vector of r values (all arbitrarily positive)
