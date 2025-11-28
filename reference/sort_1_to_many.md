# Sort results for 1-to-many forest plot

This function sorts user-supplied results for the
[`forest_plot_1_to_many()`](https://mrcieu.github.io/TwoSampleMR/reference/forest_plot_1_to_many.md)
function. The user supplies their results in the form of a data frame.

## Usage

``` r
sort_1_to_many(
  mr_res,
  b = "b",
  trait_m = "outcome",
  sort_action = 4,
  group = NULL,
  priority = NULL
)
```

## Arguments

- mr_res:

  Data frame of results supplied by the user.

- b:

  Name of the column specifying the effect of the exposure on the
  outcome. The default is `"b"`.

- trait_m:

  The column specifying the names of the traits. Corresponds to 'many'
  in the 1-to-many forest plot. The default is `"outcome"`.

- sort_action:

  Choose how to sort results.

  - `sort_action = 1`: sort results by effect size within groups. Use
    the group order supplied by the user.

  - `sort_action = 2`: sort results by effect size and group. Overrides
    the group ordering supplied by the user.

  - `sort_action = 3`: group results for the same trait together (e.g.
    multiple results for the same trait from different MR methods).

  - `sort_action = 4`: sort by decreasing effect size (largest effect
    size at top and smallest at bottom).

  - `sort_action = 5`: sort by increasing effect size (smallest effect
    size at top and largest at bottom).

- group:

  Name of grouping variable in `mr_res`.

- priority:

  If `sort_action = 3`, choose which value of the `trait_m` variable
  should be given priority and go above the other `trait_m` values. The
  trait with the largest effect size for the prioritised group will go
  to the top of the plot.

## Value

data frame.
