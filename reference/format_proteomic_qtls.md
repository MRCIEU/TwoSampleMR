# Get data from proteomic QTL results

See
[`format_data()`](https://mrcieu.github.io/TwoSampleMR/reference/format_data.md).

## Usage

``` r
format_proteomic_qtls(proteomic_qtls_subset, type = "exposure")
```

## Arguments

- proteomic_qtls_subset:

  Selected rows from `proteomic_qtls` data loaded from `MRInstruments`
  package.

- type:

  Are these data used as `"exposure"` or `"outcome"`? Default is
  `"exposure"`.

## Value

Data frame
