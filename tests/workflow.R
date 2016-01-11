ssh -L 3306:localhost:3306 gh13047@epi-franklin.epi.bris.ac.uk


# From GWAS catalog and mysql
load("../mr_base_shiny/data/gwas_catalog.RData")
bmi <- subset(gwas_catalog, Phenotype=="Body mass index" & Year==2010 & grepl("kg", Units), select=c(SNP, Effect, eaf, Allele, other_allele, SE))
bmi <- subset(gwas_catalog, Phenotype=="Body mass index", select=c(SNP, Effect, eaf, Allele, other_allele, SE))
names(bmi) <- c("SNP", "beta", "eaf", "effect_allele", "other_allele", "se")

exposure_dat <- format_exposure_dat(bmi, "BMI")
outcome_dat <- extract_outcome_data(bmi, 1:7)
dat <- harmonise_exposure_outcome(exposure_dat, outcome_dat)
mr_results <- mr(dat)



# From text files
exposure_dat <- read_exposure_data("inst/data/telomere_length.txt", "Telomere length")
outcome_dat <- rbind(
	read_outcome_data("inst/data/cardiogram.txt", "Cardiogram"),
	read_outcome_data("inst/data/bladdercancer.txt", "Bladder cancer")
)
dat <- harmonise_exposure_outcome(exposure_dat, outcome_dat)
mr_results <- mr(dat)





# From text files and mysql
exposure_dat <- read_exposure_data("inst/data/telomere_length.txt", "Telomere length")
ao <- available_outcomes()
outcome_dat <- extract_outcome_data_using_get(exposure_dat, ao$id[1:7])
outcome_dat2 <- extract_outcome_data(exposure_dat, ao$id[1:7])
dat <- harmonise_exposure_outcome(exposure_dat, outcome_dat)
mr_results <- mr(dat)
mrs <- mr_singlesnp(dat)




## Other stuff


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






mrs <- mr_leaveoneout(dat)
m <- mr_singlesnp(dat)
p <- mr_leaveoneout_plot(mrs)
p1 <- mr_forest_plot(m)

a <- mr_scatter_plot(mr_results,dat)
mr_leaveoneout_plot
mr_funnel_plot

# exposure_dat <- ld_pruning(exposure_dat)
