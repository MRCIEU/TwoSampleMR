# Two Sample MR

Extended documentation is available here:

[https://mrcieu.github.io/TwoSampleMR/](https://mrcieu.github.io/TwoSampleMR/)

* * *

Two sample Mendelian randomisation is a technique that makes causal inference about an exposure on an outcome using only summary statistics from a GWAS. This means you obtain SNPs (the instruments) that are robustly associated with your exposure, obtain a set of GWAS summary associations for the outcome you are interested, extract the instrument SNPs from the outcome GWAS, and by contrasting the effect sizes of the SNPs on the exposure with the effect sizes of the SNPs on the outcome one can estimate the causal effect.

This package provides the following functionality to help perform two sample MR:

- Extraction of user-specified SNP effects from a choice of hundreds of GWAS summary statistics
- LD pruning of exposure SNPs
- Harmonisation of direction of effects between exposure and outcome associations
- Two sample MR analysis methods and diagnostic tools
	- Including the MR Steiger test for orienting causal directions, described here [http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007081](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007081)
- Plots and automatically generated reports


A set of instruments from several sources including GWAS catalogs, metabolite QTLs, etc can be obtained from the `MRInstruments` package [https://github.com/MRCIEU/MRInstruments](https://github.com/MRCIEU/MRInstruments). Users can also specify instruments manually by providing a text file. The package uses the MR Base database, a host to hundreds of GWAS results, that can be queried as a resource for outcome data. Users can alternatively provide their own outcome summary associations. It also provides LD pruning on a remote server for exposure SNP data provided by the user.

## Citation

If using MR-Base or the TwoSampleMR R package:

[Hemani G, Zheng J, Wade KH, Laurin C, Elsworth B, Burgess S, Bowden J, Langdon R, Tan V, Yarmolinsky J, Shihab HA, Timpson NJ, Evans DM, Relton C, Martin RM, Davey Smith G, Gaunt TR, Haycock PC, The MR-Base Collaboration.</br>
MR-Base: a platform for systematic causal inference across the phenome using billions of genetic associations. <br/>
bioRxiv. doi: https://doi.org/10.1101/078972](http://biorxiv.org/content/early/2016/12/16/078972)

If also using the MR-Steiger test:

[Hemani G, Tilling K, Davey Smith G. 
Orienting the causal relationship between imprecisely measured traits using GWAS summary data. 
PLoS Genetics. 2017. 13(11).](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007081)

## License

This project is licensed under GNU GPL v3.

## Installing the TwoSampleMR R package

The package is hosted on github, and this allows installation and update to be very easy. First make sure you have the `biomaRt` and `devtools` packages installed:

    install.packages("devtools")

Then to install:

    library(devtools)
    install_github("MRCIEU/TwoSampleMR")

To update the package just run the `install_github("MRCIEU/TwoSampleMR")` command again.
