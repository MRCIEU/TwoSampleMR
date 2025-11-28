# Generate MR report

Using the output from the
[`mr()`](https://mrcieu.github.io/TwoSampleMR/reference/mr.md) function
this report will generate a report containing tables and graphs
summarising the results. A separate report is produced for each
exposure - outcome pair that was analysed.

## Usage

``` r
mr_report(
  dat,
  output_path = ".",
  output_type = "html",
  author = "Analyst",
  study = "Two Sample MR",
  path = system.file("reports", package = "TwoSampleMR"),
  ...
)
```

## Arguments

- dat:

  Output from
  [`harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.md)

- output_path:

  Directory in which reports should be saved.

- output_type:

  Choose `"html"` or `"md"`. Default is `"html"`. All output files
  including cache and figures will appear in the folder specified in
  `output_path`.

- author:

  Author name.

- study:

  Study title.

- path:

  The filepath to the report template.

- ...:

  Extra options to be passed to
  [`knitr::knit()`](https://rdrr.io/pkg/knitr/man/knit.html).
