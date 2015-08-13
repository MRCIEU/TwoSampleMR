load("inst/data/tel.RData")


# exposure_dat <- read_exposure_data(filename, "Telomere length")
# exposure_dat <- ld_pruning(exposure_dat)
# outcome_dat <- extract_outcome_data(exposure_dat, c("cardiogram", "bladder cancer"))
# dat <- harmonise_exposure_outcome(exposure_dat, outcome_dat)
res <- mr(dat)
# knit_mr_summary(res)
