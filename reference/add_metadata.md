# Add meta data to extracted data

Previously the meta data was returned alongside association information.
This is mostly unnecessary as it is needlessly repeating information.
This is a convenience function that reinstates that information. Can be
applied to either exposure data, outcome data, or harmonised data

## Usage

``` r
add_metadata(dat, cols = c("sample_size", "ncase", "ncontrol", "unit", "sd"))
```

## Arguments

- dat:

  Either exposure data, outcome data or harmonised data

- cols:

  Which metadata fields to add. Default =
  `c("sample_size", "ncase", "ncontrol", "unit", "sd")`

## Value

Data frame
