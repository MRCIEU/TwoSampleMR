# Two Sample MR

Two sample Mendelian randomisation is a technique that makes causal inference about an exposure on an outcome using only summary statistics from a GWAS. This means you obtain SNPs (the instruments) that are robustly associated with your exposure, obtain a set of GWAS summary associations for the outcome you are interested, extract the instrument SNPs from the outcome GWAS, and by contrasting the effect sizes of the SNPs on the exposure with the effect sizes of the SNPs on the outcome one can estimate the causal effect.

This package provides the following functionality to help perform two sample MR:

- Extraction of user-specified SNP effects from a choice of hundreds of GWAS summary statistics
- LD pruning of exposure SNPs
- Harmonisation of direction of effects between exposure and outcome associations
- Two sample MR analysis methods and diagnostic tools
- Plots and automatically generated reports

A set of instruments from several sources including GWAS catalogs, metabolite QTLs, etc can be obtained from the `MRInstruments` package [https://github.com/MRCIEU/MRInstruments](https://github.com/MRCIEU/MRInstruments). Users can also specify instruments manually by providing a text file. The package uses the MR Base database, a host to hundreds of GWAS results, that can be queried as a resource for outcome data. Users can alternatively provide their own outcome summary associations. It also provides LD pruning on a remote server for exposure SNP data provided by the user. 

## Installing the TwoSampleMR R package

The package is hosted on github, and this allows installation and update to be very easy. First make sure you have the `biomaRt` and `devtools` packages installed:

    install.packages("devtools")

Then to install:

    library(devtools)
    install_github("MRCIEU/TwoSampleMR")

To update the package just run the `install_github("MRCIEU/TwoSampleMR")` command again.


## Using the package

Load library

    library(TwoSampleMR)

Define your exposure (i.e. genetic proxies/instruments for body mass index) 

    bmi_file <- system.file("data/bmi.txt", package="TwoSampleMR")
    exposure_dat <- read_exposure_data(bmi_file)

Alternatively, instruments can be identified from various data sources in the  `MRInstruments` package [https://github.com/MRCIEU/MRInstruments](https://github.com/MRCIEU/MRInstruments)., e.g. metabolomic QTLs: 

    devtools::install_github("MRCIEU/MRInstruments")
    library(MRInstruments)
    data(metab_qtl) #to load metabolomic QTLs
    exposure_dat <- format_metab_qtls(metab_qtls) 
    
To use the MR-Base GWAS catalog to define instruments:
  ao<-available_outcomes() 
  exposure_dat <- extract_instruments(ao$id[c(1)]) 
 
Prune the SNPs in LD using clumping on the remote server:

    exposure_dat <- clump_data(exposure_dat)

Get the available outcomes in MR Base

    ao <- available_outcomes()

Extract the outcome associations, e.g. for Celiac disease and T2D
    
    outcome_dat <- extract_outcome_data(exposure_dat$SNP, c(6, 13))

Harmonise the exposure and outcome data. This checks that the effect alleles in the exposure and outcome data align, flips them and the effect size directions when necessary, and drops SNPs when it is impossible to determine the correct orientation.
    
    dat <- harmonise_data(exposure_dat, outcome_dat)

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

Funnel plot

    p <- mr_funnel_plot(s)
    p[[1]]
    p[[2]]

Also possible to perform heterogeneity tests, tests for directional horizontal pleiotropy, tests for the causal direction, enrichment analysis, and more.

More details are outlined in the vignette.
