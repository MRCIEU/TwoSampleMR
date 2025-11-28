# Mendelian randomization with GWAS summary data

A package for performing Mendelian randomization using GWAS summary
data. It uses the [IEU OpenGWAS database](https://gwas.mrcieu.ac.uk/) to
obtain data automatically, and a wide range of methods to run the
analysis.

## January 2020 major update

**We have made substantial changes to the package, database and
reference panels.** For full details of the changes, please visit
<https://mrcieu.github.io/TwoSampleMR/articles/gwas2020.html>

## Installation

Users running Windows and macOS, to install the latest version of
TwoSampleMR please install from our MRC IEU r-universe

``` r
install.packages("TwoSampleMR", repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))
```

Users running Linux or WebR please see the [following
instructions](https://github.com/MRCIEU/mrcieu.r-universe.dev#readme).

To update the package run the same command again.

### Installing from source

``` r
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")
```

To update the package just run the
`remotes::install_github("MRCIEU/TwoSampleMR")` command again.

## Docker

A multi-platform docker image containing R with the TwoSampleMR package
pre-installed (for both x86_64 and ARM computers) is available here:
<https://hub.docker.com/r/mrcieu/twosamplemr>
