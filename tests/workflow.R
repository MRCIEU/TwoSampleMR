ssh -L 3306:localhost:3306 gh13047@epi-franklin.epi.bris.ac.uk


# system.file("data/telomere_length.txt", package="meffil")


exposure_dat <- read_exposure_data("inst/data/telomere_length.txt", "Telomere length")
outcome_dat <- rbind(
	read_outcome_data("inst/data/cardiogram.txt", "Cardiogram"),
	read_outcome_data("inst/data/bladdercancer.txt", "Bladder cancer")
)

# exposure_dat <- ld_pruning(exposure_dat)

ao <- available_outcomes()
outcome_dat <- extract_outcome_data(exposure_dat, "GIANT_2015_HIPadjBMI_COMBINED_AllAncestries.txt.gz.uniform.af.txt")

outcome_dat <- extract_outcome_data(exposure_dat, ao$filename[30:37])



dat <- harmonise_exposure_outcome(exposure_dat, outcome_dat)

d <- subset(dat, outcome=="CD.gwas_ichip_meta_release.txt.gz.uniform.af.txt")
mr_results <- mr(dat)

m <- mr(d)

mrs <- mr_leaveoneout(dat)

p <- mr_leaveoneout_plot(mrs)


mr_report(mr_results, dat, path="inst/reports", output_path=".")


ggplot(mrs[[3]], aes(x=SNP, y=b)) + 
geom_point() + 
geom_errorbar(aes(ymin=b-se, ymax=b+se),width=0, colour=grey) +
coord_flip()






exposure_dat <- read_exposure_data("inst/data/telomere_length.txt", "Telomere length")

outcome_dat <- rbind(
	read_outcome_data("inst/data/cardiogram.txt", "Cardiogram"),
	read_outcome_data("inst/data/bladdercancer.txt", "Bladder cancer")
)

ao <- available_outcomes()
outcome_dat <- extract_outcome_data(exposure_dat, ao$filename[1:7])
dat <- harmonise_exposure_outcome(exposure_dat, outcome_dat)
mr_results <- mr(dat)
mrs <- mr_leaveoneout(dat)

m <- mr_singlesnp(dat)

p <- mr_leaveoneout_plot(mrs)
p1 <- mr_forest_plot(m)



a <- mr_scatter_plot(mr_results,dat)
mr_leaveoneout_plot
mr_funnel_plot

