# Two Sample MR

This package provides the following functionality to help perform two sample MR:

- Extraction of user-specified SNP effects from a choice of GWAS summary statistics
- LD pruning of exposure SNPs
- Harmonisation of direction of effects between exposure and outcome results
- Two sample MR analysis
- Plots and automatically generated reports 

It uses a database that hosts hundreds of GWAS results that can be queried as a resource for outcome data. It also provides LD pruning on a remote server for exposure SNP data provided by the user.

## Installing the TwoSampleMR R package

The package is hosted on github, and this allows installation and update to be very easy. First make sure you have the `biomaRt` and `devtools` packages installed:

    source("https://bioconductor.org/biocLite.R")
    biocLite("biomaRt")
    install.packages("devtools")

Then to install:

    library(devtools)
    install_github("MRCIEU/TwoSampleMR")

To update the package just run the `install_github("MRCIEU/TwoSampleMR")` command again.


## Using the package

Load library

    library(TwoSampleMR)

Read the data (e.g. telomere length example data)

    tl_file <- system.file("data/telomere_length.txt", package="TwoSampleMR")
    exposure_dat <- read_exposure_data(tl_file, "Telomere length")

Prune the SNPs in LD using clumping on the remote server:

    exposure_dat <- clump_data(exposure_dat)

Get the available outcome studies

    ao <- available_outcomes()

Extract the outcome associations, e.g. for Celiac disease and T2D
    
    outcome_dat <- extract_outcome_data(exposure_dat, c(6, 13))

Harmonise the exposure and outcome data
    
    dat <- harmonise_exposure_outcome(exposure_dat, outcome_dat)

Perform the MR
    
    mr_results <- mr(dat)

Plot the different methods

    p <- mr_scatter_plot(mr_results, dat)
    p[[1]]
    p[[2]]

Leave one out analysis

    l <- mr_leaveoneout(dat)
    p <- mr_leaveoneout_plot(l)
    p[[1]]
    p[[2]]

Forest plot

    s <- mr_singlesnp(dat)
    p <- mr_forest_plot(s)
    p[[1]]
    p[[2]]

More details about each step are outlined in the vignette.


