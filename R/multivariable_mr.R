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
multivariable_mr <- function(id_exposure, id_outcome, harmonise_strictness=2)
{
	require(reshape2)
	message("Warning: This analysis is still experimental")
	message("Testing the regression based multivariable MR")
	stopifnot(length(id_exposure) > 1)
	stopifnot(length(id_outcome) == 1)

	# Get best instruments for each exposure
	exposure_dat <- extract_instruments(id_exposure, r2 = 0.0000001, kb=10000)
	temp <- exposure_dat
	temp$id.exposure <- 1
	temp <- clump_data(temp, clump_r2=0.0000001)
	exposure_dat <- subset(exposure_dat, SNP %in% temp$SNP)


	# Get effects of each instrument from each exposure
	d1 <- extract_outcome_data(exposure_dat$SNP, id_exposure)
	d1 <- subset(d1, mr_keep.outcome)
	d2 <- subset(d1, id.outcome != id_exposure[1])
	d1 <- convert_outcome_to_exposure(subset(d1, id.outcome == id_exposure[1]))

	# Harmonise against the first id
	d <- harmonise_data(d1, d2, action=harmonise_strictness)

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
	dat <- harmonise_data(d1, outcome_dat, action=harmonise_strictness)
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
		mod <- summary(lm(lm(dat$beta.outcome ~ exposure_mat[,-c(i)])$res ~ exposure_mat[,i]))
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


#' Convert outcome format to exposure format
#'
#' @param outcome_dat Output from \code{format_data(type="outcome")}
#'
#' @export
#' @return Data frame
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

#' Extract exposure variables for multivariable MR
#'
#' Requires a list of IDs from \code{available_outcomes()}. For each ID, it extracts instruments. Then, it gets the full list of all instruments and extracts those SNPs for every exposure. Finally, it keeps only the SNPs that are a) independent and b) present in all exposures, and harmonises them to be all on the same strand. 
#'
#' @param id_exposure Array of IDs (e.g. c(299, 300, 302) for HDL, LDL, trigs)
#' @param clump_r2=0.01 Once a full list of
#' @param clump_kb=10000 <what param does>
#'
#' @export
#' @return data frame in exposure_dat format
mv_extract_exposures <- function(id_exposure, clump_r2=0.01, clump_kb=10000, harmonise_strictness=2)
{
	require(reshape2)
	message("Warning: This analysis is still experimental")
	message("Testing the regression based multivariable MR")
	stopifnot(length(id_exposure) > 1)

	# Get best instruments for each exposure
	exposure_dat <- extract_instruments(id_exposure, r2 = clump_r2, kb=clump_kb)
	temp <- exposure_dat
	temp$id.exposure <- 1
	temp <- clump_data(temp, clump_r2=clump_r2, clump_kb=clump_kb)
	exposure_dat <- subset(exposure_dat, SNP %in% temp$SNP)


	# Get effects of each instrument from each exposure
	d1 <- extract_outcome_data(exposure_dat$SNP, id_exposure)
	d1 <- subset(d1, mr_keep.outcome)
	d2 <- subset(d1, id.outcome != id_exposure[1])
	d1 <- convert_outcome_to_exposure(subset(d1, id.outcome == id_exposure[1]))

	# Harmonise against the first id
	d <- harmonise_data(d1, d2, action=harmonise_strictness)

	# Only keep SNPs that are present in all
	tab <- table(d$SNP)
	keepsnps <- names(tab)[tab == length(id_exposure)-1]
	d <- subset(d, SNP %in% keepsnps)
	
	# Reshape exposures
	dh1 <- subset(d, id.outcome == id.outcome[1], select=c(SNP, exposure, id.exposure, effect_allele.exposure, other_allele.exposure, eaf.exposure, beta.exposure, se.exposure, pval.exposure))
	dh2 <- subset(d, select=c(SNP, outcome, id.outcome, effect_allele.outcome, other_allele.outcome, eaf.outcome, beta.outcome, se.outcome, pval.outcome))
	names(dh2) <- gsub("outcome", "exposure", names(dh2))
	dh <- rbind(dh1, dh2)
	return(dh)
}


#' Harmonise exposure and outcome for multivariable MR
#'
#'
#' @param exposure_dat Output from \code{mv_extract_exposures}
#' @param outcome_dat Output from \code{extract_outcome_data(exposure_dat$SNP, id_output)}
#'
#' @export
#' @return List of vectors and matrices required for mv analysis. exposure_beta is a matrix of beta coefficients, rows correspond to SNPs and columns correspond to exposures. exposure_pval is the same as exposure_beta, but for p-values. exposure_se is the same as exposure_beta, but for standard errors. outcome_beta is an array of effects for the outcome, corresponding to the SNPs in exposure_beta. outcome_se and outcome_pval are as in outcome_beta.
mv_harmonise_data <- function(exposure_dat, outcome_dat, harmonise_strictness=2)
{

	stopifnot(all(c("SNP", "exposure", "effect_allele.exposure", "beta.exposure", "se.exposure", "pval.exposure") %in% names(exposure_dat)))
	nexp <- length(unique(exposure_dat$exposure))
	stopifnot(nexp > 1)
	tab <- table(exposure_dat$SNP)
	keepsnp <- names(tab)[tab == nexp]
	exposure_dat <- subset(exposure_dat, SNP %in% keepsnp)


	exposure_mat <- reshape2::dcast(exposure_dat, SNP ~ exposure, value.var="beta.exposure")


	# Get outcome data
	dat <- harmonise_data(subset(exposure_dat, exposure == exposure_dat$exposure[1]), outcome_dat, action=harmonise_strictness)
	dat <- subset(dat, mr_keep)
	dat$SNP <- as.character(dat$SNP)

	exposure_beta <- reshape2::dcast(exposure_dat, SNP ~ exposure, value.var="beta.exposure")
	exposure_beta <- subset(exposure_beta, SNP %in% dat$SNP)
	exposure_beta$SNP <- as.character(exposure_beta$SNP)

	exposure_pval <- reshape2::dcast(exposure_dat, SNP ~ exposure, value.var="pval.exposure")
	exposure_pval <- subset(exposure_pval, SNP %in% dat$SNP)
	exposure_pval$SNP <- as.character(exposure_pval$SNP)

	exposure_se <- reshape2::dcast(exposure_dat, SNP ~ exposure, value.var="se.exposure")
	exposure_se <- subset(exposure_se, SNP %in% dat$SNP)
	exposure_se$SNP <- as.character(exposure_se$SNP)

	index <- match(exposure_beta$SNP, dat$SNP)
	dat <- dat[index, ]
	stopifnot(all(dat$SNP == exposure_beta$SNP))

	exposure_beta <- as.matrix(exposure_beta[,-1])
	exposure_pval <- as.matrix(exposure_pval[,-1])
	exposure_se <- as.matrix(exposure_se[,-1])

	rownames(exposure_beta) <- dat$SNP
	rownames(exposure_pval) <- dat$SNP
	rownames(exposure_se) <- dat$SNP

	outcome_beta <- dat$beta.outcome
	outcome_se <- dat$se.outcome
	outcome_pval <- dat$pval.outcome

	return(list(exposure_beta=exposure_beta, exposure_pval=exposure_pval, exposure_se=exposure_se, outcome_beta=outcome_beta, outcome_pval=outcome_pval, outcome_se=outcome_se))
}


#' Perform basic multivariable MR
#'
#' Performs initial multivariable MR analysis from Burgess et al 2015.
#'
#' @param mvdat Output from \code{mv_harmonise_data}
#' @param pval_threshold=5e-8 P-value threshold to include instruments
#'
#' @export
#' @return List of results
mv_basic <- function(mvdat, pval_threshold=5e-8)
{
	# This is a matrix of 
	beta.outcome <- mvdat$outcome_beta
	beta.exposure <- mvdat$exposure_beta
	pval.exposure <- mvdat$exposure_pva;

	nexp <- ncol(beta.exposure)
    effs <- array(1:nexp)
    se <- array(1:nexp)
    pval <- array(1:nexp)
    nsnp <- array(1:nexp)
    marginal_outcome <- matrix(0, nrow(beta.exposure), ncol(beta.exposure))
    p <- list()
    nom <- colnames(beta.exposure)
    for (i in 1:nexp) {

    	# For this exposure, only keep SNPs that meet some p-value threshold
    	index <- pval.exposure[,i] < pval_threshold

    	# Get outcome effects adjusted for all effects on all other exposures
    	marginal_outcome[,i] <- lm(beta.outcome ~ beta.exposure[, -c(i)])$res

    	# Get the effect of the exposure on the residuals of the outcome
		mod <- summary(lm(marginal_outcome[index,i] ~ beta.exposure[index, i]))

		effs[i] <- mod$coef[2, 1]
		se[i] <- mod$coef[2, 2]
		pval[i] <- pnorm(abs(effs[i])/se[i], lower.tail = FALSE)
		nsnp[i] <- sum(index)

		# Make scatter plot
		d <- data.frame(outcome=marginal_outcome[,i], exposure=beta.exposure[,i])
		flip <- sign(d$exposure) == -1
		d$outcome[flip] <- d$outcome[flip] * -1
		d$exposure <- abs(d$exposure)
		p[[i]] <- ggplot2::ggplot(d[index,], ggplot2::aes(x=exposure, y=outcome)) +
		ggplot2::geom_point() +
		# geom_abline(intercept=0, slope=effs[i]) +
		ggplot2::stat_smooth(method="lm") +
		ggplot2::labs(x=paste0("SNP effect on ", nom[i]), y="Marginal SNP effect on outcome")
    }

    return(list(result=data.frame(exposure = nom, nsnp = nsnp, b = effs, se = se, pval = pval, stringsAsFactors = FALSE), marginal_outcome=marginal_outcome, plots=p))
}
