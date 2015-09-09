ssh -L 3306:localhost:3306 mruser@epi-franklin.epi.bris.ac.uk


# system.file("data/telomere_length.txt", package="meffil")


exposure_dat <- read_exposure_data("inst/data/telomere_length.txt", "Telomere length")
outcome_dat <- rbind(
	read_outcome_data("inst/data/cardiogram.txt", "Cardiogram"),
	read_outcome_data("inst/data/bladdercancer.txt", "Bladder cancer")
)

# exposure_dat <- ld_pruning(exposure_dat)

ao <- available_outcomes()
outcome_dat <- extract_outcome_data(exposure_dat, ao$filename[1:4])

dat <- harmonise_exposure_outcome(exposure_dat, outcome_dat)
mr_results <- mr(dat)

mr_report(mr_results, dat, path="inst/reports", output_path=".")


