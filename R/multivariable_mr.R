#' Extract exposure variables for multivariable MR
#'
#' Requires a list of IDs from available_outcomes. For each ID, it extracts instruments. Then, it gets the full list of all instruments and extracts those SNPs for every exposure. Finally, it keeps only the SNPs that are a) independent and b) present in all exposures, and harmonises them to be all on the same strand. 
#'
#' @param id_exposure Array of IDs (e.g. c(299, 300, 302) for HDL, LDL, trigs)
#' @param clump_r2 The default is `0.01`.
#' @param clump_kb The default is `10000`.
#' @param harmonise_strictness See the `action` option of [`harmonise_data`]. The default is `2`.
#' @param access_token Google OAuth2 access token. Used to authenticate level of access to data.
#' @param find_proxies Look for proxies? This slows everything down but is more accurate. The default is `TRUE`.
#' @param force_server Whether to search through pre-clumped dataset or to re-extract and clump directly from the server. The default is `FALSE`.
#' @param pval_threshold Instrument detection p-value threshold. Default = 5e-8
#'
#' @export
#' @return data frame in `exposure_dat` format
mv_extract_exposures <- function(id_exposure, clump_r2=0.001, clump_kb=10000, harmonise_strictness=2, access_token = ieugwasr::check_access_token(), find_proxies=TRUE, force_server=FALSE, pval_threshold=5e-8)
{
	requireNamespace("reshape2", quietly = TRUE)
	stopifnot(length(id_exposure) > 1)
	id_exposure <- ieugwasr::legacy_ids(id_exposure)

	# Get best instruments for each exposure
	exposure_dat <- extract_instruments(id_exposure, p1 = pval_threshold, r2 = clump_r2, kb=clump_kb, access_token = access_token, force_server=force_server)
	temp <- exposure_dat
	temp$id.exposure <- 1
	temp <- temp[order(temp$pval.exposure, decreasing=FALSE), ]
	temp <- subset(temp, !duplicated(SNP))
	temp <- clump_data(temp, clump_p1=pval_threshold, clump_r2=clump_r2, clump_kb=clump_kb)
	exposure_dat <- subset(exposure_dat, SNP %in% temp$SNP)


	# Get effects of each instrument from each exposure
	d1 <- extract_outcome_data(exposure_dat$SNP, id_exposure, access_token = access_token, proxies=find_proxies)
	stopifnot(length(unique(d1$id)) == length(unique(id_exposure)))
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



#' Attempt to perform MVMR using local data
#'
#' Under construction
#'
#' @param filenames_exposure Filenames for each exposure dataset. Must have header with at least SNP column present. Following arguments are used for determining how to read the filename and clumping etc.
#' @param sep Specify delimeter in file. The default is space, i.e. `sep=" "`.
#' @param phenotype_col Optional column name for the column with phenotype name corresponding the the SNP. If not present then will be created with the value `"Outcome"`. Default is `"Phenotype"`.
#' @param snp_col Required name of column with SNP rs IDs. The default is `"SNP"`.
#' @param beta_col Required for MR. Name of column with effect sizes. THe default is `"beta"`.
#' @param se_col Required for MR. Name of column with standard errors. The default is `"se"`.
#' @param eaf_col Required for MR. Name of column with effect allele frequency. The default is `"eaf"`.
#' @param effect_allele_col Required for MR. Name of column with effect allele. Must be "A", "C", "T" or "G". The default is `"effect_allele"`.
#' @param other_allele_col Required for MR. Name of column with non effect allele. Must be "A", "C", "T" or "G". The default is `"other_allele"`.
#' @param pval_col Required for enrichment tests. Name of column with p-value. The default is `"pval"`.
#' @param units_col Optional column name for units. The default is `"units"`.
#' @param ncase_col Optional column name for number of cases. The default is `"ncase"`.
#' @param ncontrol_col Optional column name for number of controls. The default is `"ncontrol"`.
#' @param samplesize_col Optional column name for sample size. The default is `"samplesize"`.
#' @param gene_col Optional column name for gene name. The default is `"gene"`.
#' @param id_col Optional column name to give the dataset an ID. Will be generated automatically if not provided for every trait / unit combination. The default is `"id"`.
#' @param min_pval Minimum allowed p-value. The default is `1e-200`.
#' @param log_pval The pval is -log10(P). The default is `FALSE`.
#' @param pval_threshold Default=5e-8 for clumping
#' @param clump_r2 Default=0.001 for clumping
#' @param clump_kb Default=10000 for clumping
#' @param harmonise_strictness See action argument in harmonise_data. Default=2
#'
#' @export
#' @return List
mv_extract_exposures_local <- function(filenames_exposure, sep = " ", phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "beta", se_col = "se", eaf_col = "eaf", effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "pval", units_col = "units", ncase_col = "ncase", ncontrol_col = "ncontrol", samplesize_col = "samplesize", gene_col = "gene", id_col = "id", min_pval = 1e-200, log_pval = FALSE, pval_threshold=5e-8, clump_r2=0.001, clump_kb=10000, harmonise_strictness=2)
{
	message("WARNING: Experimental function")
	l_full <- list()
	l_inst <- list()
	for(i in 1:length(filenames_exposure))
	{
		l_full[[i]] <- read_outcome_data(filenames_exposure[i], 
			sep = sep,
			phenotype_col = phenotype_col,
			snp_col = snp_col,
			beta_col = beta_col,
			se_col = se_col,
			eaf_col = eaf_col,
			effect_allele_col = effect_allele_col,
			other_allele_col = other_allele_col,
			pval_col = pval_col,
			units_col = units_col,
			ncase_col = ncase_col,
			ncontrol_col = ncontrol_col,
			samplesize_col = samplesize_col,
			gene_col = gene_col,
			id_col = id_col,
			min_pval = min_pval,
			log_pval = log_pval
		)
		l_inst[[i]] <- subset(l_full[[i]], pval.outcome < pval_threshold)
		l_inst[[i]] <- convert_outcome_to_exposure(l_inst[[i]])
		l_inst[[i]] <- clump_data(l_inst[[i]], clump_p1=pval_threshold, clump_r2=clump_r2, clump_kb=clump_kb)
	}

	exposure_dat <- dplyr::bind_rows(l_inst)
	id_exposure <- unique(exposure_dat$id.exposure)
	temp <- exposure_dat
	temp$id.exposure <- 1
	temp <- temp[order(temp$pval.exposure, decreasing=FALSE), ]
	temp <- subset(temp, !duplicated(SNP))
	temp <- clump_data(temp, clump_p1=pval_threshold, clump_r2=clump_r2, clump_kb=clump_kb)
	exposure_dat <- subset(exposure_dat, SNP %in% temp$SNP)

	d1 <- lapply(l_full, function(x) {
		subset(x, SNP %in% exposure_dat$SNP)
		}) %>% dplyr::bind_rows()

	stopifnot(length(unique(d1$id)) == length(unique(id_exposure)))
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
#' @md
#' @param exposure_dat Output from [`mv_extract_exposures`].
#' @param outcome_dat Output from `extract_outcome_data(exposure_dat$SNP, id_output)`.
#' @param harmonise_strictness See the `action` option of [`harmonise_data`]. The default is `2`.
#'
#' @export
#' @return List of vectors and matrices required for mv analysis. 
#' \describe{
#' \item{exposure_beta}{a matrix of beta coefficients, in which rows correspond to SNPs and columns correspond to exposures.}
#' \item{exposure_se}{is the same as `exposure_beta`, but for standard errors.}
#' \item{exposure_pval}{the same as `exposure_beta`, but for p-values.}
#' \item{expname}{A data frame with two variables, `id.exposure` and `exposure` which are character strings.}
#' \item{outcome_beta}{an array of effects for the outcome, corresponding to the SNPs in exposure_beta.}
#' \item{outcome_se}{an array of standard errors for the outcome.}
#' \item{outcome_pval}{an array of p-values for the outcome.}
#' \item{outname}{A data frame with two variables, `id.outcome` and `outcome` which are character strings.}
#' }
#' 
mv_harmonise_data <- function(exposure_dat, outcome_dat, harmonise_strictness=2)
{

	stopifnot(all(c("SNP", "id.exposure", "exposure", "effect_allele.exposure", "beta.exposure", "se.exposure", "pval.exposure") %in% names(exposure_dat)))
	nexp <- length(unique(exposure_dat$id.exposure))
	stopifnot(nexp > 1)
	tab <- table(exposure_dat$SNP)
	keepsnp <- names(tab)[tab == nexp]
	exposure_dat <- subset(exposure_dat, SNP %in% keepsnp)


	exposure_mat <- reshape2::dcast(exposure_dat, SNP ~ id.exposure, value.var="beta.exposure")


	# Get outcome data
	dat <- harmonise_data(subset(exposure_dat, id.exposure == exposure_dat$id.exposure[1]), outcome_dat, action=harmonise_strictness)
	dat <- subset(dat, mr_keep)
	dat$SNP <- as.character(dat$SNP)

	exposure_beta <- reshape2::dcast(exposure_dat, SNP ~ id.exposure, value.var="beta.exposure")
	exposure_beta <- subset(exposure_beta, SNP %in% dat$SNP)
	exposure_beta$SNP <- as.character(exposure_beta$SNP)

	exposure_pval <- reshape2::dcast(exposure_dat, SNP ~ id.exposure, value.var="pval.exposure")
	exposure_pval <- subset(exposure_pval, SNP %in% dat$SNP)
	exposure_pval$SNP <- as.character(exposure_pval$SNP)

	exposure_se <- reshape2::dcast(exposure_dat, SNP ~ id.exposure, value.var="se.exposure")
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

	expname <- subset(exposure_dat, !duplicated(id.exposure), select=c(id.exposure, exposure))
	outname <- subset(outcome_dat, !duplicated(id.outcome), select=c(id.outcome, outcome))


	return(list(exposure_beta=exposure_beta, exposure_pval=exposure_pval, exposure_se=exposure_se, outcome_beta=outcome_beta, outcome_pval=outcome_pval, outcome_se=outcome_se, expname=expname, outname=outname))
}


#' Perform basic multivariable MR
#'
#' Performs initial multivariable MR analysis from Burgess et al 2015. For each exposure the outcome is residualised for all the other exposures, then unweighted regression is applied.
#'
#' @md
#' @param mvdat Output from [`mv_harmonise_data`].
#' @param intercept Should the intercept by estimated (`TRUE`) or force line through the origin (`FALSE`, default).
#' @param instrument_specific Should the estimate for each exposure be obtained by using all instruments from all exposures (`FALSE`, default) or by using only the instruments specific to each exposure (`TRUE`).
#' @param pval_threshold P-value threshold to include instruments. The default is `5e-8`.
#' @param plots Create plots? The default is `FALSE`.
#'
#' @export
#' @return List of results
#' @importFrom stats lm pnorm
mv_residual <- function(mvdat, intercept=FALSE, instrument_specific=FALSE, pval_threshold=5e-8, plots=FALSE)
{
	# This is a matrix of 
	beta.outcome <- mvdat$outcome_beta
	beta.exposure <- mvdat$exposure_beta
	pval.exposure <- mvdat$exposure_pval

	nexp <- ncol(beta.exposure)
	effs <- array(1:nexp)
	se <- array(1:nexp)
	pval <- array(1:nexp)
	nsnp <- array(1:nexp)
	marginal_outcome <- matrix(0, nrow(beta.exposure), ncol(beta.exposure))
	p <- list()
	nom <- colnames(beta.exposure)
	nom2 <- mvdat$expname$exposure[match(nom, mvdat$expname$id.exposure)]
	for (i in 1:nexp) {

		# For this exposure, only keep SNPs that meet some p-value threshold
		index <- pval.exposure[,i] < pval_threshold

		# Get outcome effects adjusted for all effects on all other exposures
		if(intercept)
		{
			if(instrument_specific)
			{
				marginal_outcome[index,i] <- lm(beta.outcome[index] ~ beta.exposure[index, -c(i), drop=FALSE])$res
				mod <- summary(lm(marginal_outcome[index,i] ~ beta.exposure[index, i]))
			} else {
				marginal_outcome[,i] <- lm(beta.outcome ~ beta.exposure[, -c(i), drop=FALSE])$res
				mod <- summary(lm(marginal_outcome[,i] ~ beta.exposure[,i]))
			}
		} else {
			if(instrument_specific)
			{
				marginal_outcome[index,i] <- lm(beta.outcome[index] ~ 0 + beta.exposure[index, -c(i), drop=FALSE])$res
				mod <- summary(lm(marginal_outcome[index,i] ~ 0 + beta.exposure[index, i]))
			} else {
				marginal_outcome[,i] <- lm(beta.outcome ~ 0 + beta.exposure[, -c(i), drop=FALSE])$res
				mod <- summary(lm(marginal_outcome[,i] ~ 0 + beta.exposure[,i]))
			}			
		}
		if(sum(index) > (nexp + as.numeric(intercept)))
		{
			effs[i] <- mod$coef[as.numeric(intercept) + 1, 1]
			se[i] <- mod$coef[as.numeric(intercept) + 1, 2]
		} else {
			effs[i] <- NA
			se[i] <- NA
		}
		pval[i] <- 2 * pnorm(abs(effs[i])/se[i], lower.tail = FALSE)
		nsnp[i] <- sum(index)

		# Make scatter plot
		d <- data.frame(outcome=marginal_outcome[,i], exposure=beta.exposure[,i])
		flip <- sign(d$exposure) == -1
		d$outcome[flip] <- d$outcome[flip] * -1
		d$exposure <- abs(d$exposure)
		if(plots)
		{
			p[[i]] <- ggplot2::ggplot(d[index,], ggplot2::aes(x=exposure, y=outcome)) +
			ggplot2::geom_point() +
			ggplot2::geom_abline(intercept=0, slope=effs[i]) +
			# ggplot2::stat_smooth(method="lm") +
			ggplot2::labs(x=paste0("SNP effect on ", nom2[i]), y="Marginal SNP effect on outcome")
		}
	}
	result <- data.frame(id.exposure = nom, id.outcome = mvdat$outname$id.outcome, outcome=mvdat$outname$outcome, nsnp = nsnp, b = effs, se = se, pval = pval, stringsAsFactors = FALSE)
	result <- merge(mvdat$expname, result)
	out <- list(
		result=result,
		marginal_outcome=marginal_outcome
	)

	if(plots) out$plots <- p
	return(out)
}



#' Perform IVW multivariable MR
#'
#' Performs modified multivariable MR analysis. For each exposure the instruments are selected then all exposures for those SNPs are regressed against the outcome together, weighting for the inverse variance of the outcome.
#' 
#' @md
#' @param mvdat Output from [`mv_harmonise_data`].
#' @param intercept Should the intercept by estimated (`TRUE`) or force line through the origin (`FALSE`, default).
#' @param instrument_specific Should the estimate for each exposure be obtained by using all instruments from all exposures (`FALSE`, default) or by using only the instruments specific to each exposure (`TRUE`).
#' @param pval_threshold P-value threshold to include instruments. The default is `5e-8`.
#' @param plots Create plots? The default is `FALSE`.
#'
#' @export
#' @return List of results
#' @importFrom stats lm pnorm
mv_multiple <- function(mvdat, intercept=FALSE, instrument_specific=FALSE, pval_threshold=5e-8, plots=FALSE)
{
	# This is a matrix of 
	beta.outcome <- mvdat$outcome_beta
	beta.exposure <- mvdat$exposure_beta
	pval.exposure <- mvdat$exposure_pval
	w <- 1/mvdat$outcome_se^2

	nexp <- ncol(beta.exposure)
	effs <- array(1:nexp)
	se <- array(1:nexp)
	pval <- array(1:nexp)
	nsnp <- array(1:nexp)
	# marginal_outcome <- matrix(0, nrow(beta.exposure), ncol(beta.exposure))
	p <- list()
	nom <- colnames(beta.exposure)
	nom2 <- mvdat$expname$exposure[match(nom, mvdat$expname$id.exposure)]
	for (i in 1:nexp)
	{
		# For this exposure, only keep SNPs that meet some p-value threshold
		index <- pval.exposure[,i] < pval_threshold

		# # Get outcome effects adjusted for all effects on all other exposures
		# marginal_outcome[,i] <- lm(beta.outcome ~ beta.exposure[, -c(i)])$res

		# Get the effect of the exposure on the residuals of the outcome
		if(!intercept)
		{
			if(instrument_specific)
			{
				mod <- summary(lm(beta.outcome[index] ~ 0 + beta.exposure[index, ,drop=FALSE], weights=w[index]))
			} else {
				mod <- summary(lm(beta.outcome ~ 0 + beta.exposure, weights=w))
			}
		} else {
			if(instrument_specific)
			{
				mod <- summary(lm(beta.outcome[index] ~ beta.exposure[index, ,drop=FALSE], weights=w[index]))
			} else {
				mod <- summary(lm(beta.outcome ~ beta.exposure, weights=w))
			}
		}

		if(instrument_specific & sum(index) <= (nexp + as.numeric(intercept)))
		{
			effs[i] <- NA
			se[i] <- NA
		} else {
			effs[i] <- mod$coef[as.numeric(intercept) + i, 1]
			se[i] <- mod$coef[as.numeric(intercept) + i, 2]
		}
		pval[i] <- 2 * pnorm(abs(effs[i])/se[i], lower.tail = FALSE)
		nsnp[i] <- sum(index)

		# Make scatter plot
		d <- data.frame(outcome=beta.outcome, exposure=beta.exposure[,i])
		flip <- sign(d$exposure) == -1
		d$outcome[flip] <- d$outcome[flip] * -1
		d$exposure <- abs(d$exposure)
		if(plots)
		{
			p[[i]] <- ggplot2::ggplot(d[index,], ggplot2::aes(x=exposure, y=outcome)) +
			ggplot2::geom_point() +
			ggplot2::geom_abline(intercept=0, slope=effs[i]) +
			# ggplot2::stat_smooth(method="lm") +
			ggplot2::labs(x=paste0("SNP effect on ", nom2[i]), y="Marginal SNP effect on outcome")
		}
	}
	result <- data.frame(id.exposure = nom, id.outcome = mvdat$outname$id.outcome, outcome=mvdat$outname$outcome, nsnp = nsnp, b = effs, se = se, pval = pval, stringsAsFactors = FALSE)
	result <- merge(mvdat$expname, result)
	out <- list(
		result=result
	)
	if(plots)
		out$plots=p

	return(out)
}

#' Perform basic multivariable MR
#' 
#' Performs initial multivariable MR analysis from Burgess et al 2015. For each exposure the outcome is residualised for all the other exposures, then unweighted regression is applied.
#'
#' @md
#' @param mvdat Output from [`mv_harmonise_data`].
#' @param pval_threshold P-value threshold to include instruments. The default is `5e-8`.
#'
#' @export
#' @return List of results
#' @importFrom stats lm pnorm
mv_basic <- function(mvdat, pval_threshold=5e-8)
{
	# This is a matrix of 
	beta.outcome <- mvdat$outcome_beta
	beta.exposure <- mvdat$exposure_beta
	pval.exposure <- mvdat$exposure_pval

	nexp <- ncol(beta.exposure)
	effs <- array(1:nexp)
	se <- array(1:nexp)
	pval <- array(1:nexp)
	nsnp <- array(1:nexp)
	marginal_outcome <- matrix(0, nrow(beta.exposure), ncol(beta.exposure))
	p <- list()
	nom <- colnames(beta.exposure)
	nom2 <- mvdat$expname$exposure[match(nom, mvdat$expname$id.exposure)]
	for (i in 1:nexp) {

		# For this exposure, only keep SNPs that meet some p-value threshold
		index <- pval.exposure[,i] < pval_threshold

		# Get outcome effects adjusted for all effects on all other exposures
		marginal_outcome[,i] <- lm(beta.outcome ~ beta.exposure[, -c(i)])$res

		# Get the effect of the exposure on the residuals of the outcome
		mod <- summary(lm(marginal_outcome[index,i] ~ beta.exposure[index, i]))

		effs[i] <- mod$coef[2, 1]
		se[i] <- mod$coef[2, 2]
		pval[i] <- 2 * pnorm(abs(effs[i])/se[i], lower.tail = FALSE)
		nsnp[i] <- sum(index)

		# Make scatter plot
		d <- data.frame(outcome=marginal_outcome[,i], exposure=beta.exposure[,i])
		flip <- sign(d$exposure) == -1
		d$outcome[flip] <- d$outcome[flip] * -1
		d$exposure <- abs(d$exposure)
		p[[i]] <- ggplot2::ggplot(d[index,], ggplot2::aes(x=exposure, y=outcome)) +
		ggplot2::geom_point() +
		ggplot2::geom_abline(intercept=0, slope=effs[i]) +
		# ggplot2::stat_smooth(method="lm") +
		ggplot2::labs(x=paste0("SNP effect on ", nom2[i]), y="Marginal SNP effect on outcome")
	}
	result <- data.frame(id.exposure = nom, id.outcome = mvdat$outname$id.outcome, outcome=mvdat$outname$outcome, nsnp = nsnp, b = effs, se = se, pval = pval, stringsAsFactors = FALSE)
	result <- merge(mvdat$expname, result)

	return(list(result=result, marginal_outcome=marginal_outcome, plots=p))
}



#' Perform IVW multivariable MR
#'
#' Performs modified multivariable MR analysis. For each exposure the instruments are selected then all exposures for those SNPs are regressed against the outcome together, weighting for the inverse variance of the outcome.
#'
#' @md
#' @param mvdat Output from [`mv_harmonise_data`].
#' @param pval_threshold P-value threshold to include instruments. The default is `5e-8`.
#'
#' @export
#' @return List of results
#' @importFrom stats pnorm
mv_ivw <- function(mvdat, pval_threshold=5e-8)
{
	# This is a matrix of 
	beta.outcome <- mvdat$outcome_beta
	beta.exposure <- mvdat$exposure_beta
	pval.exposure <- mvdat$exposure_pval
	w <- 1/mvdat$outcome_se^2

	nexp <- ncol(beta.exposure)
	effs <- array(1:nexp)
	se <- array(1:nexp)
	pval <- array(1:nexp)
	nsnp <- array(1:nexp)
	# marginal_outcome <- matrix(0, nrow(beta.exposure), ncol(beta.exposure))
	p <- list()
	nom <- colnames(beta.exposure)
	nom2 <- mvdat$expname$exposure[match(nom, mvdat$expname$id.exposure)]
	for (i in 1:nexp) {

		# For this exposure, only keep SNPs that meet some p-value threshold
		index <- pval.exposure[,i] < pval_threshold

		# # Get outcome effects adjusted for all effects on all other exposures
		# marginal_outcome[,i] <- lm(beta.outcome ~ beta.exposure[, -c(i)])$res

		# Get the effect of the exposure on the residuals of the outcome
		mod <- summary(lm(beta.outcome[index] ~ 0 + beta.exposure[index, ], weights=w[index]))

		effs[i] <- mod$coef[i, 1]
		se[i] <- mod$coef[i, 2]
		pval[i] <- 2 * pnorm(abs(effs[i])/se[i], lower.tail = FALSE)
		nsnp[i] <- sum(index)

		# Make scatter plot
		d <- data.frame(outcome=beta.outcome, exposure=beta.exposure[,i])
		flip <- sign(d$exposure) == -1
		d$outcome[flip] <- d$outcome[flip] * -1
		d$exposure <- abs(d$exposure)
		p[[i]] <- ggplot2::ggplot(d[index,], ggplot2::aes(x=exposure, y=outcome)) +
		ggplot2::geom_point() +
		ggplot2::geom_abline(intercept=0, slope=effs[i]) +
		# ggplot2::stat_smooth(method="lm") +
		ggplot2::labs(x=paste0("SNP effect on ", nom2[i]), y="Marginal SNP effect on outcome")
	}
	result <- data.frame(id.exposure = nom, id.outcome = mvdat$outname$id.outcome, outcome=mvdat$outname$outcome, nsnp = nsnp, b = effs, se = se, pval = pval, stringsAsFactors = FALSE)
	result <- merge(mvdat$expname, result)

	return(list(result=result, plots=p))
}

#' Apply LASSO feature selection to mvdat object
#'
#' @param mvdat Output from [`mv_harmonise_data`].
#'
#' @export
#' @return data frame of retained features
#' @importFrom glmnet cv.glmnet coef.glmnet
mv_lasso_feature_selection <- function(mvdat)
{
	message("Performing feature selection")
	b <- glmnet::cv.glmnet(x=mvdat$exposure_beta, y=mvdat$outcome_beta, weight=1/mvdat$outcome_se^2, intercept=0)
	c <- glmnet::coef.glmnet(b, s = "lambda.min")
  	i <- !c[,1] == 0
  	d <- dplyr::tibble(exposure=rownames(c)[i], b=c[i,])
	return(d)
}

#' Perform multivariable MR on subset of features
#' 
#' The function proceeds as follows:
#' \enumerate{
#' \item Select features (by default this is done using LASSO feature selection).
#' \item Subset the mvdat to only retain relevant features and instruments.
#' \item Perform MVMR on remaining data.
#' }
#' @md
#' @param mvdat Output from [`mv_harmonise_data`].
#' @param features Dataframe of features to retain, must have column with name 'exposure' that has list of exposures tor etain from mvdat. The default is `mvdat_lasso_feature_selection(mvdat)`.
#' @param intercept Should the intercept by estimated (`TRUE`) or force line through the origin (`FALSE`, the default).
#' @param instrument_specific Should the estimate for each exposure be obtained by using all instruments from all exposures (`FALSE`, default) or by using only the instruments specific to each exposure (`TRUE`).
#' @param pval_threshold P-value threshold to include instruments. The default is `5e-8`.
#' @param plots Create plots? The default is `FALSE`.
#'
#' @export
#' @return List of results
#' @importFrom stats lm
mv_subset <- function(mvdat, features=mv_lasso_feature_selection(mvdat), intercept=FALSE, instrument_specific=FALSE, pval_threshold=5e-8, plots=FALSE)
{
	# Update mvdat object
	mvdat$exposure_beta <- mvdat$exposure_beta[, features$exposure, drop=FALSE]
	mvdat$exposure_se <- mvdat$exposure_se[, features$exposure, drop=FALSE]
	mvdat$exposure_pval <- mvdat$exposure_pval[, features$exposure, drop=FALSE]

	# Find relevant instruments
	instruments <- apply(mvdat$exposure_pval, 1, function(x) any(x < pval_threshold))
	stopifnot(sum(instruments) > nrow(features))

	mvdat$exposure_beta <- mvdat$exposure_beta[instruments,,drop=FALSE]
	mvdat$exposure_se <- mvdat$exposure_se[instruments,,drop=FALSE]
	mvdat$exposure_pval <- mvdat$exposure_pval[instruments,,drop=FALSE]	
	mvdat$outcome_beta <- mvdat$outcome_beta[instruments]
	mvdat$outcome_se <- mvdat$outcome_se[instruments]
	mvdat$outcome_pval <- mvdat$outcome_pval[instruments]

	mv_multiple(mvdat, intercept=intercept, instrument_specific=instrument_specific, pval_threshold=pval_threshold, plots=plots)
}

