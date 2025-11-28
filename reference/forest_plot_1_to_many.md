# 1-to-many forest plot

Plot results from an analysis of multiple exposures against a single
outcome or a single exposure against multiple outcomes. Plots effect
estimates and 95 percent confidence intervals. The ordering of results
in the plot is determined by the order supplied by the user. Users may
find
[`sort_1_to_many()`](https://mrcieu.github.io/TwoSampleMR/reference/sort_1_to_many.md)
helpful for sorting their results prior to using the 1-to-many forest
plot. The plot function works best for 50 results and is not designed to
handle more than 100 results.

## Usage

``` r
forest_plot_1_to_many(
  mr_res = "mr_res",
  b = "b",
  se = "se",
  TraitM = "outcome",
  col1_width = 1,
  col1_title = "",
  exponentiate = FALSE,
  trans = "identity",
  ao_slc = TRUE,
  lo = NULL,
  up = NULL,
  by = NULL,
  xlab = "Effect (95% confidence interval)",
  addcols = NULL,
  addcol_widths = NULL,
  addcol_titles = "",
  subheading_size = 6,
  shape_points = 15,
  colour_scheme = "black",
  col_text_size = 5,
  weight = NULL
)
```

## Arguments

- mr_res:

  Data frame of results supplied by the user. The default is `"mr_res"`.

- b:

  Name of the column specifying the effect of the exposure on the
  outcome. The default is `"b"`.

- se:

  Name of the column specifying the standard error for b. The default is
  `"se"`.

- TraitM:

  The column specifying the names of the traits. Corresponds to 'many'
  in the 1-to-many forest plot. The default is `"outcome"`.

- col1_width:

  Width of Y axis label for the column specified by the TraitM argument.
  The default is `1`.

- col1_title:

  Title for the column specified by the TraitM argument. The default is
  `""`.

- exponentiate:

  Convert log odds ratios to odds ratios? Default is `FALSE`.

- trans:

  Specify x-axis scale. e.g. "identity", "log2", etc. If set to
  "identity" an additive scale is used. If set to log2 the x-axis is
  plotted on a multiplicative / doubling scale (preferable when plotting
  odds ratios). Default is `"identity"`.

- ao_slc:

  Logical; retrieve trait subcategory information using
  available_outcomes(). Default is `FALSE`.

- lo:

  Lower limit of X axis to plot.

- up:

  upper limit of X axis to plot.

- by:

  Name of the grouping variable to stratify results on. Default is
  `NULL`.

- xlab:

  X-axis label, default is `"Effect (95% confidence interval)"`.

- addcols:

  Name of additional columns to plot. Character vector. The default is
  `NULL`.

- addcol_widths:

  Widths of Y axis labels for additional columns specified by the
  addcols argument. Numeric vector. The default is `NULL`.

- addcol_titles:

  Titles of additional columns specified by the addcols argument.
  Character vector. The default is `""`.

- subheading_size:

  text size for the subheadings specified in by argument. The default is
  `6`.

- shape_points:

  the shape of the data points to pass to
  [`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html).
  Default is set to `15` (filled square).

- colour_scheme:

  the general colour scheme for the plot. Default is to make all text
  and data points `"black"`.

- col_text_size:

  The default is `5`.

- weight:

  The default is `NULL`.

## Value

grid plot object
