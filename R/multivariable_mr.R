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



extract_multivariable_instruments <- function(id)
{
	stopifnot(length(id) > 1)
	# Get best instruments for each exposure
	l <- list()
	for(i in 1:length(id))
	{
		message("Extracting instruments for ", id[i])
		l[[i]] <- extract_instruments(id[i])
	}

	exposure_dat <- rbind.fill(l)
	exposure_dat <- clump_data(exposure_dat)


	# Get effects of each instrument from each exposure

	d1 <- extract_outcome_data(exposure_dat$SNP, id[1])
	d2 <- extract_outcome_data(exposure_dat$SNP, id[-1])
	d1 <- convert_outcome_to_exposure(d1)

	d <- harmonise_data(d1, d2)

	


	od <- table(od$SNP)

	#


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

