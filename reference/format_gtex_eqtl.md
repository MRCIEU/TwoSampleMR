# Get data from eQTL catalog into correct format

See
[`format_data()`](https://mrcieu.github.io/TwoSampleMR/reference/format_data.md).

## Usage

``` r
format_gtex_eqtl(gtex_eqtl_subset, type = "exposure")
```

## Arguments

- gtex_eqtl_subset:

  Selected rows from `gtex_eqtl` data loaded from `MRInstruments`
  package.

- type:

  Are these data used as `"exposure"` or `"outcome"`? Default is
  `"exposure"`.

## Value

Data frame
