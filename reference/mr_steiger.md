# MR Steiger test of directionality

A statistical test for whether the assumption that exposure causes
outcome is valid

## Usage

``` r
mr_steiger(p_exp, p_out, n_exp, n_out, r_exp, r_out, r_xxo = 1, r_yyo = 1, ...)
```

## Arguments

- p_exp:

  Vector of p-values of SNP-exposure

- p_out:

  Vector of p-values of SNP-outcome

- n_exp:

  Sample sizes for p_exp

- n_out:

  Sample sizes for p_out

- r_exp:

  Vector of absolute correlations for SNP-exposure

- r_out:

  Vector of absolute correlations for SNP-outcome

- r_xxo:

  Measurement precision of exposure

- r_yyo:

  Measurement precision of outcome

- ...:

  Further arguments to be passed to
  [`lattice::wireframe()`](https://rdrr.io/pkg/lattice/man/cloud.html)

## Value

List with the following elements:

- r2_exp:

  Estimated variance explained in x

- r2_out:

  Estimated variance explained in y

- r2_exp_adj:

  Predicted variance explained in x accounting for estimated measurement
  error

- r2_out_adj:

  Predicted variance explained in y accounting for estimated measurement
  error

- correct_causal_direction:

  `TRUE`/`FALSE`

- steiger_test:

  p-value for inference of direction

- correct_causal_direction_adj:

  `TRUE`/`FALSE`, direction of causality for given measurement error
  parameters

- steiger_test_adj:

  p-value for inference of direction of causality for given measurement
  error parameters

- vz:

  Total volume of the error parameter space

- vz0:

  Volume of the parameter space that gives the incorrect answer

- vz1:

  Volume of the parameter space that gives the correct answer

- sensitivity_ratio:

  Ratio of vz1/vz0. Higher means inferred direction is less susceptible
  to measurement error

- sensitivity_plot:

  Plot of parameter space of causal directions and measurement error
