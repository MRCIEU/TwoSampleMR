# Knit report using working environment

Warning: It is quite likely that this will be called within an Rmd file
implying a recursive call to `knit()`. This will generate "duplicate
label" errors for unlabelled chunks. To avoid this, all code chunks in
your Rmd file should be named. Supposedly this error can also be avoided
by setting the following option:
`options(knitr.duplicate.label = 'allow')`. I tried this but it didn't
seem to help.

## Usage

``` r
knit_report(input_filename, output_filename, ...)
```

## Arguments

- input_filename:

  Rmd file.

- output_filename:

  Markdown or HTML output file. An HTML file is specified using the
  .htm, .html, .HTM or .HTML file extension. When html is specified, a
  similarly named markdown file is also generated. All output files
  including cache and figures will appear in the same folder as
  `output_filename`.

- ...:

  Arguments to be passed to
  [`knitr::knit()`](https://rdrr.io/pkg/knitr/man/knit.html)
