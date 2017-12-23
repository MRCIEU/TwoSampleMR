devs <- function()
{

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

}

outlier_scan <- function(dat, use_proxies=FALSE, threshold = 0.01)
{
	# Get outliers
	radial <- RadialMR(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP, "IVW", "YES", "NO", 0.05/nrow(dat), "NO")

	# Find associations with outliers

	ao <- available_outcomes()
	ids <- subset(ao, priority == 1 & nsnp > 500000 & sample_size > 5000) %>%
		arrange(desc(sample_size)) %>%
		filter(!duplicated(trait), mr == 1) %>%
		filter(! id %in% c(dat$id.exposure[1], dat$id.outcome[1]))

	out <- extract_outcome_data(radial$outliers$SNP, ids$id, proxies=use_proxies)
	out2 <- subset(out, pval.outcome < 5e-8)

	inst <- extract_instruments(unique(out2$id.outcome))

	exp1 <- extract_outcome_data(unique(inst$SNP), dat$id.exposure[1], proxies=use_proxies)
	out1 <- extract_outcome_data(unique(inst$SNP), dat$id.outcome[1], proxies=use_proxies)

	dat1 <- harmonise_data(inst, exp1)
	dat2 <- harmonise_data(inst, out1)

	mr1 <- mr(dat1)
	mr11 <- group_by(mr1, id.exposure, id.outcome) %>%
		do({
			x <- .
			print(nrow(x))
			if(x$nsnp[1] == 1)
			{
				x
			} else if(x$nsnp[1] < 5) {
				subset(x, method == "Inverse variance weighted")
			} else {
				subset(x, method == "Weighted mode")
			}
		})
	mr2 <- mr(dat2)
	mr22 <- group_by(mr2, id.exposure, id.outcome) %>%
		do({
			x <- .
			print(nrow(x))
			if(x$nsnp[1] == 1)
			{
				x
			} else if(x$nsnp[1] < 5) {
				subset(x, method == "Inverse variance weighted")
			} else {
				subset(x, method == "Weighted mode")
			}
		})

	# Use threshold to retain causal relationships

	# Make graph of how each thing relates to exposure and outcome

	return(outlier_list)
}



find_associations <- function(outlier_list)
{
	
}


