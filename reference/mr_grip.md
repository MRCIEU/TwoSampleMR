# MR-GRIP: a modified MR-Egger model with the Genotype Recoding Invariance Property

This implements the modified MR-Egger model with the Genotype Recoding
Invariance Property (MR-GRIP) due to Dudbridge and Bowden et al. (2025).
It is well known that the results of MR-Egger are sensitive to which
alleles are designated as the effect alleles. A pragmatic convention is
to orient all SNPs to have positive effects on the exposure, which has
some advantages in interpretation but also brings some philosophical
limitations. The MR-GRIP model is a modification to the MR-Egger model
in which each term is multiplied by the genotype-phenotype associations.
This makes each term in the model invariant to allele coding.

## Usage

``` r
mr_grip(b_exp, b_out, se_exp, se_out, parameters)
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

  List of parameters.

## Value

List of with the following elements:

- b:

  MR estimate

- se:

  Standard error of MR estimate

- pval:

  p-value of MR estimate

- b_i:

  Intercept

- se_i:

  Standard error of intercept

- pval_i:

  p-value of intercept

- b.adj:

  MR estimate adjusting for weak instruments

- se.adj:

  Standard error adjusting for weak instruments

- pval.adj:

  p-value adjusting for weak instruments

- mod:

  Summary of regression

- dat:

  Original data used for MR-GRIP
