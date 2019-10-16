#' Convert TwoSampleMR format to MendelianRandomization format
#'
#' The MendelianRandomization package offers MR methods that can 
#' be used with the same data used in the TwoSampleMR package. This
#' function converts from the TwoSampleMR format to the MRInput class.
#'
#' @param dat Output from the \code{harmonise_data} function
#' @param get_correlations Default FALSE. If TRUE then extract the LD matrix for the SNPs from the European 1000 genomes data on the MR-Base server
#'
#' @export
#' @return List of MRInput objects for each exposure/outcome combination
dat_to_MRInput <- function(dat, get_correlations=FALSE)
{
	library(MendelianRandomization)
	out <- plyr::dlply(dat, c("exposure", "outcome"), function(x)
	{
		x <- plyr::mutate(x)
		message("Converting:")
		message(" - exposure: ", x$exposure[1])
		message(" - outcome: ", x$outcome[1])
		if(get_correlations)
		{
			message(" - obtaining LD matrix")
			ld <- ld_matrix(unique(x$SNP))
			out <- harmonise_ld_dat(x, ld)
			if(is.null(out))
			{
				return(NULL)
			}
			x <- out$x
			ld <- out$ld

			MendelianRandomization::mr_input(
				bx = x$beta.exposure,
				bxse = x$se.exposure,
				by = x$beta.outcome,
				byse = x$se.outcome,
				exposure = x$exposure[1],
				outcome = x$outcome[1],
				snps = x$SNP,
				effect_allele=x$effect_allele.exposure,
				other_allele=x$other_allele.exposure,
				eaf = x$eaf.exposure,
				correlation = ld
			)

		} else {
			MendelianRandomization::mr_input(
				bx = x$beta.exposure,
				bxse = x$se.exposure,
				by = x$beta.outcome,
				byse = x$se.outcome,
				exposure = x$exposure[1],
				outcome = x$outcome[1],
				snps = x$SNP,
				effect_allele=x$effect_allele.exposure,
				other_allele=x$other_allele.exposure,
				eaf = x$eaf.exposure
			)
		}
	})
	return(out)
}



harmonise_ld_dat <- function(x, ld)
{
	snpnames <- do.call(rbind, strsplit(rownames(ld), split="_"))
	stopifnot(all(snpnames[,1] == x$SNP))
	x$effect_allele.exposure <- as.character(x$effect_allele.exposure)
	x$other_allele.exposure <- as.character(x$other_allele.exposure)
	# Set1 x and ld alleles match
	snpnames <- data.frame(snpnames, stringsAsFactors=FALSE)
	snpnames <- merge(subset(x, select=c(SNP, effect_allele.exposure, other_allele.exposure)), snpnames, by.x="SNP", by.y="X1")
	snpnames$keep <- (snpnames$X2 == snpnames$effect_allele.exposure & snpnames$X3 == snpnames$other_allele.exposure) |
		(snpnames$X3 == snpnames$effect_allele.exposure & snpnames$X2 == snpnames$other_allele.exposure)

	# What happens if everything is gone?
	if(nrow(x) == 0)
	{
		message(" - none of the SNPs could be aligned to the LD reference panel")
		return(NULL)
	}

	if(any(!snpnames$keep))
	{
		message(" - the following SNPs could not be aligned to the LD reference panel: \n", paste(subset(snpnames, keep)$SNP, collapse="\n - "))
	}


	snpnames$flip1 <- snpnames$X2 != snpnames$effect_allele.exposure
	x <- subset(x, SNP %in% snpnames$SNP)
	temp1 <- x$effect_allele.exposure[snpnames$flip1]
	temp2 <- x$other_allele.exposure[snpnames$flip1]
	x$beta.exposure[snpnames$flip1] <- x$beta.exposure[snpnames$flip1] * -1
	x$beta.outcome[snpnames$flip1] <- x$beta.outcome[snpnames$flip1] * -1
	x$effect_allele.exposure[snpnames$flip1] <- temp2
	x$other_allele.exposure[snpnames$flip1] <- temp1

	rownames(ld) <- snpnames$SNP
	colnames(ld) <- snpnames$SNP

	return(list(x=x, ld=ld))
}



#' Wrapper for MR-PRESSO
#'
#' See https://github.com/rondolab/MR-PRESSO
#'
#' @param dat Output from harmonise_data
#' @param NbDistribution = 1000 Number of bootstraps
#' @param SignifThreshold = 0.05 Outlier significance threshold
#'
#' @export
#' @return List of results for every exposure/outcome combination
run_mr_presso <- function(dat, NbDistribution = 1000,  SignifThreshold = 0.05)
{
	require(MRPRESSO)
	dat <- subset(dat, mr_keep)
	d <- subset(dat, !duplicated(paste(id.exposure, " - ", id.outcome)), select=c(exposure, outcome, id.exposure, id.outcome))
	res <- list()
	attributes(res)$id.exposure <- d$id.exposure
	attributes(res)$id.outcome <- d$id.outcome
	attributes(res)$exposure <- d$exposure
	attributes(res)$outcome <- d$id.exposure
	for(j in 1:nrow(d))
	{
		x <- subset(dat, exposure == d$exposure[j] & outcome == d$outcome[j])
		message(x$exposure[1], " - ", x$outcome[1])
		res[[j]] <- MRPRESSO::mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = x, NbDistribution = NbDistribution,  SignifThreshold = SignifThreshold)
	}
	return(res)
}

#' Convert dat to RadialMR format
#'
#' Creates a list of RadialMR format datasets for each exposure - outcome pair
#'
#' @param dat Output from \code{harmonise_data}
#'
#' @export
#' @return List of RadialMR format datasets
dat_to_RadialMR <- function(dat)
{
	require(RadialMR)
	out <- plyr::dlply(dat, c("exposure", "outcome"), function(x)
	{
		x <- plyr::mutate(x)
		message("Converting:")
		message(" - exposure: ", x$exposure[1])
		message(" - outcome: ", x$outcome[1])
		d <- subset(x, mr_keep=TRUE)
		d <- format_radial(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, RSID=d$SNP)
		return(d)
	})
	return(out)
}

#' Radial IVW analysis
#' 
#' @param b_exp Vector of genetic effects on exposure
#' @param b_out Vector of genetic effects on outcome
#' @param se_exp Standard errors of genetic effects on exposure
#' @param se_out Standard errors of genetic effects on outcome
#'
#' @export
#' @return List with the following elements:
#'         b: causal effect estimate
#'         se: standard error
#'         pval: p-value
mr_ivw_radial <- function(b_exp, b_out, se_exp, se_out, parameters=default_parameters())
{
	require(RadialMR)
	if (sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) &
		!is.na(se_out)) < 2)
		return(list(b = NA, se = NA, pval = NA, nsnp = NA))
	d <- format_radial(BXG=b_exp, BYG=b_out, seBXG=se_exp, seBYG=se_out, RSID=1:length(b_exp))
	out <- ivw_radial(d, alpha=0.05, weights=3)
	b <- out$coef[1,1]
	se <- out$coef[1,2]
	pval <- 2 * pnorm(abs(b/se), low = FALSE)
	Q_df <- out$df
	Q <- out$qstatistic
	Q_pval <- pchisq(Q, Q_df, low=F)
	return(list(b = b, se = se, pval = pval, nsnp = length(b_exp),
        Q = Q, Q_df = Q_df, Q_pval = Q_pval))
}

