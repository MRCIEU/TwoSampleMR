# system.file("data/telomere_length.txt", package="meffil")

exposure_dat <- read_exposure_data("inst/data/telomere_length.txt", "Telomere length")

# outcome_dat <- rbind(
# 	read_outcome_data("inst/data/cardiogram.txt", "Cardiogram"),
# 	read_outcome_data("inst/data/bladdercancer.txt", "Bladder cancer")
# )

# exposure_dat <- ld_pruning(exposure_dat)

ao <- available_outcomes(user='epxjz', password='wv-92n_YjB')
outcomes <- extract_outcome_data(exposure_dat, ao$filename[1:4], user='epxjz', password='wv-92n_YjB')

dat <- harmonise_exposure_outcome(exposure_dat, outcomes)
mr_res <- mr(dat)

# knit_mr_summary(res)


