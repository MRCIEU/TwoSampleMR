tl_file <- system.file("data/telomere_length.txt", package="TwoSampleMR")
exposure_dat <- read_exposure_data(tl_file, "Telomere length")
exposure_dat <- clump_data(exposure_dat)
ao <- available_outcomes()
outcome_dat <- extract_outcome_data(exposure_dat, c(6, 13))


dat <- harmonise_exposure_outcome(exposure_dat, outcome_dat)
dato <- harmonise_exposure_outcome(exposure_dat, outcome_dat)
mr_results <- mr(dat)
mr_results0 <- mr(dato)

p <- mr_scatter_plot(mr_results, dat)
p <- mr_scatter_plot(mr_results0, dato)
p[[1]]
p[[2]]
l <- mr_leaveoneout(dat)
p <- mr_leaveoneout_plot(l)
p[[1]]
p[[2]]
s <- mr_singlesnp(dat)
p <- mr_forest_plot(s)
p[[1]]
p[[2]]



bmi <- subset(gwas_catalog, Phenotype=="Body mass index" & Year==2010 & grepl("kg", Units))
exposure_dat <- format_gwas_catalog(bmi)





library(TwoSampleMR)

exposure_dat <- read_exposure_data("~/Downloads/BMI.txt", "BMI")
exposure_dat <- clump_data(exposure_dat)

ao <- available_outcomes()
outcome_dat <- extract_outcome_data(exposure_dat, ao$id[33:35])
dim(outcome_dat)



exposure_dat <- read_exposure_data("~/Downloads/HDL.txt", "HDL")
exposure_dat <- clump_data(exposure_dat)

ao <- available_outcomes()
outcome_dat <- extract_outcome_data(exposure_dat, ao$id[33:35])
dim(outcome_dat)



exposure_dat <- read_exposure_data("~/Downloads/smoking_behaviour.txt", "smoking_behaviour")
exposure_dat <- clump_data(exposure_dat)

ao <- available_outcomes()
outcome_dat <- extract_outcome_data(exposure_dat, ao$id[33:35])
dim(outcome_dat)
