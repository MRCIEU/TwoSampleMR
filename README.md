
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Mendelian randomization with GWAS summary data

<!-- badges: start -->

[![Build
Status](https://github.com/MRCIEU/TwoSampleMR/workflows/R-CMD-check/badge.svg)](https://github.com/MRCIEU/TwoSampleMR/actions?workflow=R-CMD-check)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![DOI](https://zenodo.org/badge/49515156.svg)](https://zenodo.org/badge/latestdoi/49515156)
[![Codecov test
coverage](https://codecov.io/gh/MRCIEU/TwoSampleMR/branch/ieugwasr/graph/badge.svg)](https://codecov.io/gh/MRCIEU/TwoSampleMR?branch=ieugwasr)
<!-- badges: end -->

A package for performing Mendelian randomization using GWAS summary
data. It uses the [IEU GWAS database](https://gwas.mrcieu.ac.uk/) to
obtain data automatically, and a wide range of methods to run the
analysis. You can use the [MR-Base web app](http://www.mrbase.org/) to
try out a limited range of the functionality in this package, but for
any serious work we strongly recommend using this R package.

## January 2020 major update

**We have made substantial changes to the package, database and
reference panels.** For full details of the changes, please visit
<https://mrcieu.github.io/TwoSampleMR/articles/gwas2020.html>

## Installation

To install the latest version of TwoSampleMR, perform as normal:

``` r
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")
```

To update the package just run the
`remotes::install_github("MRCIEU/TwoSampleMR")` command again.

We recommend using this new version going forwards but for a limited
time we are enabling backwards compatibility, in case you are in the
middle of analysis or need to reproduce old analysis. In order to use
the legacy version of the package and the database, install using:

``` r
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR@0.4.26")
```

## Docker

A docker image containing R with the TwoSampleMR package pre-installed
is available here: <https://hub.docker.com/r/mrcieu/twosamplemr>

<!-- Additional content -->

**Full documentation available here:**
<https://mrcieu.github.io/TwoSampleMR>
