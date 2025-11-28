# Get data from methylation QTL results

See
[`format_data()`](https://mrcieu.github.io/TwoSampleMR/reference/format_data.md).

## Usage

``` r
format_aries_mqtl(aries_mqtl_subset, type = "exposure")
```

## Arguments

- aries_mqtl_subset:

  Selected rows from `aries_mqtl` data loaded from `MRInstruments`
  package.

- type:

  Are these data used as `"exposure"` or `"outcome"`? Default is
  `"exposure"`.

## Value

Data frame
