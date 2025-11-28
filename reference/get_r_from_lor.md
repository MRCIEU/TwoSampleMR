# Estimate proportion of variance of liability explained by SNP in general population

This uses equation 10 in Lee et al. A Better Coefficient of
Determination for Genetic Profile Analysis. Genetic Epidemiology 36:
214–224 (2012)
[doi:10.1002/gepi.21614](https://doi.org/10.1002/gepi.21614) .

## Usage

``` r
get_r_from_lor(
  lor,
  af,
  ncase,
  ncontrol,
  prevalence,
  model = "logit",
  correction = FALSE
)
```

## Arguments

- lor:

  Vector of Log odds ratio.

- af:

  Vector of allele frequencies.

- ncase:

  Vector of Number of cases.

- ncontrol:

  Vector of Number of controls.

- prevalence:

  Vector of Disease prevalence in the population.

- model:

  Is the effect size estimated from the `"logit"` (default) or
  `"probit"` model.

- correction:

  Scale the estimated r by correction value. The default is `FALSE`.

## Value

Vector of signed r values
