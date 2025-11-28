# A basic forest plot

This function is used to create a basic forest plot. It requires the
output from
[`format_mr_results()`](https://mrcieu.github.io/TwoSampleMR/reference/format_mr_results.md).

## Usage

``` r
forest_plot_basic(
  dat,
  section = NULL,
  colour_group = NULL,
  colour_group_first = TRUE,
  xlab = NULL,
  bottom = TRUE,
  trans = "identity",
  xlim = NULL,
  threshold = NULL
)
```

## Arguments

- dat:

  Output from
  [`format_mr_results()`](https://mrcieu.github.io/TwoSampleMR/reference/format_mr_results.md).

- section:

  Which category in dat to plot. If `NULL` then prints everything.

- colour_group:

  Which exposure to plot. If `NULL` then prints everything grouping by
  colour.

- colour_group_first:

  The default is `TRUE`.

- xlab:

  x-axis label. Default=`NULL`.

- bottom:

  Show x-axis? Default=`FALSE`.

- trans:

  Transformation of x axis.

- xlim:

  x-axis limits.

- threshold:

  p-value threshold to use for colouring points by significance level.
  If `NULL` (default) then colour layer won't be applied.

## Value

ggplot object
