## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(eval=TRUE, cache=TRUE)
suppressWarnings(suppressPackageStartupMessages(library(TwoSampleMR)))
suppressWarnings(suppressPackageStartupMessages(library(MRInstruments)))
suppressWarnings(suppressPackageStartupMessages(library(knitr)))

## ----eval=FALSE---------------------------------------------------------------
#  library(devtools)
#  install_github("MRCIEU/TwoSampleMR")

## ----eval=FALSE---------------------------------------------------------------
#  
#  library(TwoSampleMR)
#  
#  # List available GWASs
#  ao <- available_outcomes()
#  
#  # Get instruments
#  exposure_dat <- extract_instruments("ieu-a-2")
#  
#  # Get effects of instruments on outcome
#  outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-a-7")
#  
#  # Harmonise the exposure and outcome data
#  dat <- harmonise_data(exposure_dat, outcome_dat)
#  
#  # Perform MR
#  res <- mr(dat)

