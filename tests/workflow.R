# system.file("data/telomere_length.txt", package="meffil")

exposure_dat <- read_exposure_data("inst/data/telomere_length.txt", "Telomere length")
outcome_dat <- read_outcome_data("inst/data/cardiogram.txt", "Cardiogram")

# exposure_dat <- ld_pruning(exposure_dat)
# outcome_dat <- extract_outcome_data(exposure_dat, c("cardiogram", "bladder cancer"))

dat <- harmonise_exposure_outcome(exposure_dat, outcome_dat)
mr_res <- mr(dat)

# knit_mr_summary(res)


