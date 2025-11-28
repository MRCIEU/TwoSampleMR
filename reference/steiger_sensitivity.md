# Evaluate the Steiger test's sensitivity to measurement error

Evaluate the Steiger test's sensitivity to measurement error

## Usage

``` r
steiger_sensitivity(rgx_o, rgy_o, ...)
```

## Arguments

- rgx_o:

  Observed variance of exposure explained by SNPs

- rgy_o:

  Observed variance of outcome explained by SNPs

- ...:

  Further arguments to be passed to
  [`lattice::wireframe()`](https://rdrr.io/pkg/lattice/man/cloud.html)

## Value

List with the following elements:

- vz:

  Total volume of the error parameter space

- vz0:

  Volume of the parameter space that gives the incorrect answer

- vz1:

  Volume of the parameter space that gives the correct answer

- sensitivity_ratio:

  Ratio of vz1/vz0. Higher means inferred direction is less susceptible
  to measurement error

- pl:

  plot of parameter space
