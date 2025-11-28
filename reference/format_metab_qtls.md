# Get data from metabolomic QTL results

See
[`format_data()`](https://mrcieu.github.io/TwoSampleMR/reference/format_data.md).

## Usage

``` r
format_metab_qtls(metab_qtls_subset, type = "exposure")
```

## Arguments

- metab_qtls_subset:

  Selected rows from `metab_qtls` data loaded from `MRInstruments`
  package.

- type:

  Are these data used as `"exposure"` or `"outcome"`? Default is
  `"exposure"`.

## Value

Data frame
