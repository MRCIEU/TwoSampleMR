library(devtools)
load_all()
library(RadialMR)
library(dplyr)

a <- extract_instruments(300)
b <- extract_outcome_data(a$SNP, 7)

dat <- harmonise_data(a, b)

radial$outliers
dim(dat)
mr_heterogeneity(dat)
mr(dat)

b[[1]]$Out$Pval

outlier_scan <- function(dat)
{

	# Get outliers
	radial <- RadialMR(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP, "IVW", "YES", "NO", 0.05/nrow(dat), "NO")

	# Find associations with outliers

	ao <- available_outcomes()
	ids <- subset(ao, priority == 1 & nsnp > 500000 & sample_size > 5000) %>%
		arrange(desc(sample_size)) %>%
		filter(!duplicated(trait), mr == 1)

	out <- extract_outcome_data(radial$outliers$SNP, ids$id, proxies=FALSE)
	out2 <- subset(out, pval.outcome < 5e-8)



	return(outlier_list)

	# b <- run_mr_presso(dat)

}



find_associations <- function(outlier_list)
{
	
}


