# Robust adjusted profile score

Robust adjusted profile score

## Usage

``` r
mr_raps(b_exp, b_out, se_exp, se_out, parameters = default_parameters())
```

## Arguments

- b_exp:

  Vector of genetic effects on exposure.

- b_out:

  Vector of genetic effects on outcome.

- se_exp:

  Standard errors of genetic effects on exposure.

- se_out:

  Standard errors of genetic effects on outcome.

- parameters:

  A list of parameters. Specifically, `over.dispersion`,
  `loss.function`, and `shrinkage`:

  - `over.dispersion` is a logical concerning should the model consider
    overdispersion (systematic pleiotropy);

  - `loss.function` allows using either the squared error loss (`"l2"`)
    or robust loss functions/scores (`"huber"` or `"tukey"`);

  - `shrinkage` is a logical specifying whether empirically partially
    Bayes should be used.

  The default is
  `parameters=list(overdispersion = TRUE, loss.function = "huber", shrinkage = FALSE)`.

## Value

List with the following elements:

- b:

  MR estimate

- se:

  Standard error

- pval:

  p-value

- nsnp:

  Number of SNPs

## Details

This function calls the `mr.raps` package. Please refer to the
documentation of that package for more detail.

## References

Qingyuan Zhao, Jingshu Wang, Jack Bowden, Dylan S. Small. Statistical
inference in two-sample summary-data Mendelian randomization using
robust adjusted profile score. Annals of Statistics, 48, 3, 1742–1769,
[doi:10.1214/19-AOS1866](https://doi.org/10.1214/19-AOS1866)
