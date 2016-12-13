# library(TwoSampleMR)


# set.seed(100)

# n <- 500000
# nsnp <- 6
# g <- matrix(0, n, nsnp)
# for(i in 1:nsnp)
# {
# 	g[,1] <- rbinom(n, 2, 0.5)
# 	g[,2] <- rbinom(n, 2, 0.5)
# 	g[,3] <- rbinom(n, 2, 0.5)
# 	g[,4] <- rbinom(n, 2, 0.5)
# 	g[,5] <- rbinom(n, 2, 0.5)
# 	g[,6] <- rbinom(n, 2, 0.5)
# }

# Bx1 <- rnorm(4)
# Bx2 <- rnorm(3)
# Bx3 <- rnorm(3)

# x1 <- g[,1] * Bx1[1] + g[,2] * Bx1[2] + g[,3] * Bx1[3] + g[,4] * Bx1[4] + rnorm(n) * sd(g[,1])
# x2 <- g[,4] * Bx2[1] + g[,5] * Bx2[2] + g[,6] * Bx2[3] + rnorm(n) * sd(g[,1])
# x3 <- g[,1] * Bx3[1] + g[,2] * Bx3[2] + g[,6] * Bx3[3] + rnorm(n) * sd(g[,1])

# y <- x1 * 0.5 + x2 * -5 + rnorm(n) * sd(x1)


# by <- bx1 <- bx2 <- bx3 <- array(0, 6)
# for(i in 1:6)
# {
# 	bx1[i] <- lm(x1 ~ g[,i])$coef[2]
# 	bx2[i] <- lm(x2 ~ g[,i])$coef[2]
# 	bx3[i] <- lm(x3 ~ g[,i])$coef[2]
# 	by[i]  <- lm(y ~ g[,i])$coef[2]
# }

# beta1 = lm(lm(by ~ bx2 + bx3)$res ~ bx1)$coef[2]
# beta2 = lm(lm(by ~ bx1 + bx3)$res ~ bx2)$coef[2]
# beta3 = lm(lm(by ~ bx1 + bx2)$res ~ bx3)$coef[2]
# beta1se = summary(lm(lm(by ~ bx2 + bx3)$res ~ bx1))$coef[2,2]
# beta2se = summary(lm(lm(by ~ bx1 + bx3)$res ~ bx2))$coef[2,2]
# beta3se = summary(lm(lm(by ~ bx1 + bx2)$res ~ bx3))$coef[2,2]


# beta1
# beta2
# beta3

# beta1se
# beta2se
# beta3se


# 1. identify exposures
# 2. get best instruments for each exposure
# 3. get effects of each instrument from each exposure and the outcome
# 4. create matrix of xs


#' Perform multivariable MR
#'
#' WARNING: EXPERIMENTAL, needs more testing
#' 1. identify exposures
#' 2. get best instruments for each exposure
#' 3. get effects of each instrument from each exposure and the outcome
#' 4. create matrix of exposure effects
#' 5. perform multivariable MR
#'
#' @param id_exposure Array of ids from \code{available_outcomes}
#' @param id_outcome Single id from \code{available_outcomes}
#'
#' @export
#' @return List of results table, exposure effects and outcome effects
multivariable_mr <- function(id_exposure, id_outcome)
{
	require(reshape2)
	message("Warning: This analysis is still experimental")
	stopifnot(length(id_exposure) > 1)
	stopifnot(length(id_outcome) == 1)

	# Get best instruments for each exposure
	exposure_dat <- extract_instruments(id_exposure)
	temp <- exposure_dat
	temp$id.exposure <- 1
	temp <- clump_data(temp)
	exposure_dat <- subset(exposure_dat, SNP %in% temp$SNP)


	# Get effects of each instrument from each exposure
	d1 <- extract_outcome_data(exposure_dat$SNP, id_exposure)
	d1 <- subset(d1, mr_keep.outcome)
	d2 <- subset(d1, id.outcome != id_exposure[1])
	d1 <- convert_outcome_to_exposure(subset(d1, id.outcome == id_exposure[1]))

	# Harmonise against the first id
	d <- harmonise_data(d1, d2)

	# Only keep SNPs that are present in all
	tab <- table(d$SNP)
	keepsnps <- names(tab)[tab == length(id_exposure)-1]
	d <- subset(d, SNP %in% keepsnps)
	
	# Reshape exposures
	dh1 <- subset(d, id.outcome == id.outcome[1], select=c(SNP, exposure, id.exposure, effect_allele.exposure, other_allele.exposure, beta.exposure, se.exposure))
	dh2 <- subset(d, select=c(SNP, outcome, id.outcome, effect_allele.outcome, other_allele.outcome, beta.outcome, se.outcome))
	names(dh2) <- gsub("outcome", "exposure", names(dh2))
	dh <- rbind(dh1, dh2)
	exposure_mat <- reshape2::dcast(dh, SNP ~ exposure, value.var="beta.exposure")


	# Get outcome data
	outcome_dat <- extract_outcome_data(keepsnps, id_outcome)
	dat <- harmonise_data(d1, outcome_dat)
	exposure_mat <- subset(exposure_mat, SNP %in% dat$SNP)
	dat$SNP <- as.character(dat$SNP)
	exposure_mat$SNP <- as.character(exposure_mat$SNP)
	index <- match(exposure_mat$SNP, dat$SNP)
	dat <- dat[index, ]
	stopifnot(all(dat$SNP == exposure_mat$SNP))

	exposure_mat <- as.matrix(exposure_mat[,-1])
	rownames(exposure_mat) <- dat$SNP
	effs <- array(1:length(id_exposure))
	se <- array(1:length(id_exposure))
	pval <- array(1:length(id_exposure))
	for(i in 1:length(id_exposure))
	{
		mod <- summary(
			lm(lm(dat$beta.outcome ~ exposure_mat[,-c(i)])$res ~ exposure_mat[,i]))
		effs[i] <- mod$coef[2,1]
		se[i] <- mod$coef[2,2]
		pval[i] <- pnorm(abs(effs[i]) / se[i], lower.tail=FALSE)
	}

	nom <- unique(dh$exposure)
	return(list(
		results = data.frame(
			exposure=nom,
			outcome=unique(dat$outcome),
			b=effs,
			se=se,
			pval=pval,
			stringsAsFactors=FALSE
		),
		exposure_effects = exposure_mat,
		outcome_effects  = dat$beta.outcome
	))
}


convert_outcome_to_exposure <- function(outcome_dat)
{
	exposure_dat <- format_data(
		outcome_dat,
		beta_col = "beta.outcome",
		se_col="se.outcome",
		pval_col="pval.outcome",
		phenotype_col="outcome",
		effect_allele_col="effect_allele.outcome",
		other_allele_col="other_allele.outcome",
		eaf_col="eaf.outcome",
		units_col="units.outcome"
	)
	return(exposure_dat)
}



instruments_for_multivariable_mr <- function(id_exposure, snps=NULL)
{
	require(reshape2)
	message("Warning: This analysis is still experimental")
	stopifnot(length(id_exposure) > 1)
	if(is.null(snps))
	{
		exposure_dat <- extract_instruments(id_exposure)
		temp <- exposure_dat
		temp$id.exposure <- 1
		temp <- clump_data(temp)
		exposure_dat <- subset(exposure_dat, SNP %in% temp$SNP)
	} else {
		exposure_dat <- convert_outcome_to_exposure(extract_outcome_data(snps, id_exposure))
		temp <- exposure_dat
		temp$id.exposure <- 1
		temp <- clump_data(temp)
		exposure_dat <- subset(exposure_dat, SNP %in% temp$SNP)
	}

}




# multivariable_mr <- function(id_exposure, id_outcome)
# {
# 	require(reshape2)
# 	message("Warning: This analysis is still experimental")
# 	stopifnot(length(id_exposure) > 1)
# 	stopifnot(length(id_outcome) == 1)

# 	# Get best instruments for each exposure
# 	exposure_dat <- extract_instruments(id_exposure)
# 	temp <- exposure_dat
# 	temp$id.exposure <- 1
# 	temp <- clump_data(temp)
# 	exposure_dat <- subset(exposure_dat, SNP %in% temp$SNP)


# 	# Get effects of each instrument from each exposure
	d1 <- extract_outcome_data(exposure_dat$SNP, id_exposure)
	d1 <- subset(d1, mr_keep.outcome)
	d2 <- subset(d1, id.outcome != id_exposure[1])
	d1 <- convert_outcome_to_exposure(subset(d1, id.outcome == id_exposure[1]))

# 	# Harmonise against the first id
	d <- harmonise_data(d1, d2)

# 	# Only keep SNPs that are present in all
	tab <- table(d$SNP)
	keepsnps <- names(tab)[tab == length(id_exposure)-1]
	d <- subset(d, SNP %in% keepsnps)
# }



library(MRInstruments)
data(gwas_catalog)
load("~/repo/mr_base_paper/results/lpa_ldl_trigs.RData")


hdl_inst <- format_data(subset(gwas_catalog, grepl("HDL cholest", Phenotype) & Author == "Willer CJ" & Year == 2013))
ldl_inst <- format_data(subset(gwas_catalog, grepl("LDL cholesterol", Phenotype) & Author == "Willer CJ" & Year == 2013))
trig_inst <- format_data(subset(gwas_catalog, grepl("Triglycerides", Phenotype) & Author == "Willer CJ" & Year == 2013))

hdl_inst$exposure <- "HDL cholesterol"
ldl_inst$exposure <- "LDL cholesterol"

inst <- rbind(hdl_inst, ldl_inst, trig_inst)

snps <- inst$SNP
exposure_dat 


toggle_dev("test")
ao <- available_outcomes()

temp <- ao[order(ao$sample_size, decreasing=TRUE), ]
temp <- subset(temp, !duplicated(trait))
sum(temp$sample_size, na.rm=TRUE)

