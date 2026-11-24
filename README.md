
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Mendelian randomization with GWAS summary data

<!-- badges: start -->

[![R-CMD-check](https://github.com/MRCIEU/TwoSampleMR/actions/workflows/check-full.yaml/badge.svg)](https://github.com/MRCIEU/TwoSampleMR/actions/workflows/check-full.yaml)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.3298000-blue.svg)](https://doi.org/10.5281/zenodo.3298000)
[![Codecov test
coverage](https://codecov.io/gh/MRCIEU/TwoSampleMR/branch/master/graph/badge.svg)](https://app.codecov.io/gh/MRCIEU/TwoSampleMR?branch=master)
[![TwoSampleMR status
badge](https://mrcieu.r-universe.dev/badges/TwoSampleMR)](https://mrcieu.r-universe.dev/TwoSampleMR)
<!-- badges: end -->

A package for performing Mendelian randomization using GWAS summary
data. It uses the [IEU OpenGWAS database](https://opengwas.io) to obtain
data automatically, and a wide range of methods to run the analysis.

## ALPHA version 0.8.0 - early release of new OpenGWAS infrastructure

The version on this branch replaces the `ieugwasr` package with the new `opengwasr` package, which provides access to an overhauled OpenGWAS data architecture.

## Installation

Users running Windows and macOS, to install the latest version of
TwoSampleMR please install from our MRC IEU r-universe

``` r
remotes::install_github("MRCIEU/TwoSampleMR@0.8.0")
```

## Documentation

**Full documentation available here:**
<https://mrcieu.github.io/TwoSampleMR/>
