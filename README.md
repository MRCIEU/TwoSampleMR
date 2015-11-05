# Two Sample MR

This package provides the following functionality to help perform two sample MR:

- Extraction of user-specified SNP effects from a choice of GWAS summary statistics
- Harmonisation of direction of effects between exposure and outcome results
- LD pruning of exposure SNPs
- Two sample MR analysis
- Plots and automatically generated reports 


## Accessing the MySQL server

In order to gain access to the GWAS summary statistics database we have to create a tunnel to epifranklin. In a terminal session run:

    ssh -L 3306:localhost:3306 gh13047@epi-franklin.epi.bris.ac.uk

Replace gh13047 with your own username. Leave this terminal session open and perform all additional steps in a separate terminal window.


## Installing the TwoSampleMR R package

First clone the git repository:

    git clone git@scmv-ieugit.epi.bris.ac.uk:gh13047/mr_base.git

Next, in an R session set your work directory to the `mr_base` directory that you just cloned and using the `R/devtools` library run:

    library(devtools)
    install()

This will install the package from the cloned repository. 


## Using the package

### Overview

    library(TwoSampleMR)
    
    # Read the data
    exposure_dat <- read_exposure_data("inst/data/telomere_length.txt", "Telomere length")
    
    # Get the available outcome studies
    ao <- available_outcomes()
    
    # Extract the outcome associations
    outcome_dat <- extract_outcome_data(exposure_dat, ao$id[1:7])
    
    # Harmonise the exposure and outcome data
    dat <- harmonise_exposure_outcome(exposure_dat, outcome_dat)
    
    # Perform the MR
    mr_results <- mr(dat)
    
    # Make some plots
    mr_scatter_plot(mr_results, dat)
    mr_leaveoneout(dat)
    mr_forest_plot(mr_singlesnp(dat))

More details about each step are outlined below.


### Reading in exposure data

The simplest way to read in exposure data is to read it in from a text file. To see an example of what format the textfile should take, see the file located here:

    mr_base/inst/data/telomere_length.txt

To read in this data:

    library(TwoSampleMR)
    exposure_dat <- read_exposure_data("inst/data/telomere_length.txt", "Telomere length")

This reads in the text file, checks that the data is as expected, and performs a lookup to extract chromosome and position info for each SNP.

Alternatively you can create the data frame in R manually, and then use a helper function to perform the checks. For example, using the GWAS catalog data that is available here:

    https://scmv-ieugit.epi.bris.ac.uk/gh13047/mr_base_shiny/tree/master/data

you can get the Speliotes et al 2010 SNPs by running:

    bmi <- subset(gwas_catalog, 
        Phenotype=="Body mass index" & Year==2010 & grepl("kg", Units), 
        select=c(SNP, Effect, eaf, Allele, other_allele, SE)
    )
    names(bmi) <- c("SNP", "beta", "eaf", "effect_allele", "other_allele", "se")
    exposure_dat <- format_exposure_dat(bmi, "BMI")

**TO DO:** Create a proper helper function for this and integrate GWAS catalog into this repository.


### Extracting outcome data from MySQL database

To extract summary associations for the exposure SNPs first identify the GWASs against which the SNPs will be searched. The list of all available studies can be obtained in a dataframe by:

    ao <- available_outcomes()

To extract the exposure SNPs from GWAS results for study IDs 1-7 you would then do:

    outcome_dat <- extract_outcome_data(exposure_dat, ao$id[1:7])

Alternatively, outcome data can be read in from a file directly:

    outcome_dat <- rbind(
        read_outcome_data("inst/data/cardiogram.txt", "Cardiogram"),
        read_outcome_data("inst/data/bladdercancer.txt", "Bladder cancer")
    )

### Harmonising the exposure and outcome data

To harmonise effects for exposure and outcome do this:

    dat <- harmonise_exposure_outcome(exposure_dat, outcome_dat)


### Perform MR

This performs all MR analyses and sensitivity analyses:

    mr_results <- mr(dat)


### Plots

Some plots can be created too. Leave one out analysis:

    mr_leaveoneout(dat)

Scatter plot:

    mr_scatter_plot(mr_results, dat)

Single SNP sensitivity analysis:

    singlesnp <- mr_singlesnp(dat)
    mr_forest_plot(singlesnp)


### Generate report

A report of the analysis can be generated like this:

    mr_report(mr_results, dat, path="inst/reports", output_path=".")


## To do

- Create a proper helper function for this and integrate GWAS catalog into this repository.
- Provide functionality for the LD pruning based on access to the remote server (it currently requires local access to the 1000 genomes reference data)
- Integrate Charles's forest plots
- Spruce up the report
- Perhaps restructure the MR results so that generation of graphs etc is a bit more uniform
- Improve documentation
- Create a better template structure for all the different MR analyses so that it's very simple to implement a new MR function. Currently it's quite easy but could be better
- Add in additional MR functions that allow for correlated SNPS, using Steve Burgess' likelihood and GLM approaches; these extensions will require a function that generates a correlation matrix; will also require an LD pruning step that ensures no SNP-pairs have a correlation >0.9
- MR in in non-European populations; LD pruning and corrMR functions assume European ancestry 



- Move to do list to issues tracker