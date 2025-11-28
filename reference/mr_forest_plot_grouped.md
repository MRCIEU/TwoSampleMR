# Grouped forest plot

Grouped forest plot

## Usage

``` r
mr_forest_plot_grouped(
  name,
  eff_Col = "b",
  exposure_Name = "exposure",
  outcome_Name = "outcome",
  forest_Title = "",
  outfile_Name = "annot_FP.pdf",
  left_Col_Names = c("Exposure", "Outcome"),
  left_Col_Titles = NULL,
  right_Col_Names = c("p", "Outcome.n.case", "Outcome.n.control", "Outcome.sample.size"),
  right_Col_Titles = NULL,
  debug = FALSE,
  log_ES = FALSE,
  decrease = TRUE,
  returnRobj = TRUE,
  se_Col = "se"
)
```

## Arguments

- name:

  (character) name of the delimited file containing all of the results
  on the first sheet (needs to have headers), or of the r object.

- eff_Col:

  (character) name of the column in the delimited file that contains the
  effect sizes.

- exposure_Name:

  (character) name of the column in the delimited file containing the
  *types* of studies.

- outcome_Name:

  (character) name of the column in the delimited file containing the
  names of each study.

- forest_Title:

  (character) the title to be used for a forest plot.

- outfile_Name:

  (character) name to be used for output file (*.pdf) or (*.wmf).

- left_Col_Names:

  (character vector) vector containing the names of the left-hand-side
  annotation columns in the delimited file.

- left_Col_Titles:

  (character vector) vector containing the titles for each
  left-hand-side annotation column.

- right_Col_Names:

  (character vector) vector containing the names of the right-hand-side
  annotation columns in the delimited file.

- right_Col_Titles:

  (character vector) vector containing the titles for each
  right-hand-side annotation column.

- debug:

  (logical) show warnings `TRUE`/`FALSE`?

- log_ES:

  (logical) perform natural log transform of effect sizes and confidence
  bounds `TRUE`/`FALSE`?

- decrease:

  (logical) sort the studies by decreasing effect sizes `TRUE`/`FALSE`?

- returnRobj:

  (logical) return the graph as an internal R object `TRUE`/`FALSE`?

- se_Col:

  (character) name of the column giving the standard error of the effect
  sizes.

## Value

grid object giving the forest plot (or plot as pdf)
