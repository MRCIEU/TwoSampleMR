# Format MR results for a 1-to-many forest plot

This function formats user-supplied results for the
[`forest_plot_1_to_many()`](https://mrcieu.github.io/TwoSampleMR/reference/forest_plot_1_to_many.md)
function. The user supplies their results in the form of a data frame.
The data frame is assumed to contain at least three columns of data:

1.  effect estimates, from an analysis of the effect of an exposure on
    an outcome;

2.  standard errors for the effect estimates; and

3.  a column of trait names, corresponding to the 'many' in a 1-to-many
    forest plot.

## Usage

``` r
format_1_to_many(
  mr_res,
  b = "b",
  se = "se",
  exponentiate = FALSE,
  ao_slc = FALSE,
  by = NULL,
  TraitM = "outcome",
  addcols = NULL,
  weight = NULL
)
```

## Arguments

- mr_res:

  Data frame of results supplied by the user.

- b:

  Name of the column specifying the effect of the exposure on the
  outcome. Default = `"b"`.

- se:

  Name of the column specifying the standard error for b. Default =
  `"se"`.

- exponentiate:

  Convert log odds ratios to odds ratios? Default=`FALSE`.

- ao_slc:

  Logical; retrieve trait subcategory information using
  [`available_outcomes()`](https://mrcieu.github.io/TwoSampleMR/reference/available_outcomes.md).
  Default=`FALSE`.

- by:

  Name of the column indicating a grouping variable to stratify results
  on. Default=`NULL`.

- TraitM:

  The column specifying the names of the traits. Corresponds to 'many'
  in the 1-to-many forest plot. Default=`"outcome"`.

- addcols:

  Name of any additional columns to add to the plot. Character vector.
  The default is `NULL`.

- weight:

  The default is `NULL`.

## Value

data frame.
