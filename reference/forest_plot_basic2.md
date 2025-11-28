# A basic forest plot

This function is used to create a basic forest plot. It requires the
output from
[`format_1_to_many()`](https://mrcieu.github.io/TwoSampleMR/reference/format_1_to_many.md).

## Usage

``` r
forest_plot_basic2(
  dat,
  section = NULL,
  colour_group = NULL,
  colour_group_first = TRUE,
  xlab = NULL,
  bottom = TRUE,
  trans = "identity",
  xlim = NULL,
  lo = lo,
  up = up,
  subheading_size = subheading_size,
  colour_scheme = "black",
  shape_points = 15
)
```

## Arguments

- dat:

  Output from
  [`format_1_to_many()`](https://mrcieu.github.io/TwoSampleMR/reference/format_1_to_many.md)

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

  x-axis scale.

- xlim:

  x-axis limits.

- lo:

  Lower limit of x axis.

- up:

  Upper limit of x axis.

- subheading_size:

  text size for the subheadings. The subheadings correspond to the
  values of the section argument.

- colour_scheme:

  the general colour scheme for the plot. Default is to make all text
  and data points `"black"`.

- shape_points:

  the shape of the data points to pass to
  [`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html).
  Default is set to `15` (filled square).

## Value

ggplot object
