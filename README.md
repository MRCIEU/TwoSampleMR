# Two Sample MR

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/MRCIEU/TwoSampleMR.svg?branch=ieugwasr)](https://travis-ci.org/MRCIEU/TwoSampleMR) [![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental) [![DOI](https://zenodo.org/badge/49515156.svg)](https://zenodo.org/badge/latestdoi/49515156)
[![Codecov test coverage](https://codecov.io/gh/MRCIEU/TwoSampleMR/branch/ieugwasr/graph/badge.svg)](https://codecov.io/gh/MRCIEU/TwoSampleMR?branch=ieugwasr)
<!-- badges: end -->

* * * 

This package is performing Mendelian randomization using GWAS summary data. It uses the [IEU GWAS database](https://gwas.mrcieu.ac.uk/) to obtain data automatically, and a wide range of methods to run the analysis. You can use the [MR-Base web app](http://www.mrbase.org/) to try it out a limited range of the functionality in this package, but for any serious work we strongly recommend using this R package.

**Full documentation available here:** [https://mrcieu.github.io/TwoSampleMR](https://mrcieu.github.io/TwoSampleMR/)

## January 2020 major update 

**We have made substantial changes to the package, database and reference panels.** For full details of the changes, please visit https://mrcieu.github.io/gwasglue/articles/gwas2020.html

To install the latest version of TwoSampleMR, perform as normal:

```
install.packages("devtools")
devtools::install_github("MRCIEU/TwoSampleMR")
```

To update the package just run the `install_github("MRCIEU/TwoSampleMR")` command again.

We recommend using this new version going forwards but for a limited time we are enabling backwards compatibility, in case you are in the middle of analysis or need to reproduce old analysis. In order to use the legacy version of the package and the database, install using:

```
install.packages("devtools")
devtools::install_github("MRCIEU/TwoSampleMR@0.4.25")
```

## Citation

If using MR-Base, IEU GWAS database or the TwoSampleMR R package:

[Hemani G, Zheng J, Elsworth B, Wade KH, Baird D, Haberland V, Laurin C, Burgess S, Bowden J, Langdon R, Tan VY, Yarmolinsky J, Shihab HA, Timpson NJ, Evans DM, Relton C, Martin RM, Davey Smith G, Gaunt TR, Haycock PC, The MR-Base Collaboration.</br>
**The MR-Base platform supports systematic causal inference across the human phenome.** <br/>
eLife 2018;7:e34408. doi: 10.7554/eLife.34408](https://elifesciences.org/articles/34408)

If also using the MR-Steiger test:

[Hemani G, Tilling K, Davey Smith G.<br/>
**Orienting the causal relationship between imprecisely measured traits using GWAS summary data.**<br/>
PLoS Genetics. 2017. 13(11).](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007081)

## License

This project is made available under the open source MIT license.
