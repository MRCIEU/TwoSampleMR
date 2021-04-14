#' Convert TwoSampleMR format to MendelianRandomization format
#'
#' The MendelianRandomization package offers MR methods that can 
#' be used with the same data used in the TwoSampleMR package. This
#' function converts from the TwoSampleMR format to the MRInput class.
#'
#' @param dat Output from the [`harmonise_data`] function.
#' @param get_correlations Default `FALSE`. If `TRUE` then extract the LD matrix for the SNPs from the European 1000 genomes data on the MR-Base server.
#' @param pop If get_correlations is TRUE then use the following 
#'
#' @export
#' @return List of MRInput objects for each exposure/outcome combination
dat_to_MRInput <- function(dat, get_correlations=FALSE, pop="EUR")
{
	out <- plyr::dlply(dat, c("exposure", "outcome"), function(x)
	{
		x <- plyr::mutate(x)
		message("Converting:")
		message(" - exposure: ", x$exposure[1])
		message(" - outcome: ", x$outcome[1])
		if(get_correlations)
		{
			message(" - obtaining LD matrix")
			ld <- ld_matrix(unique(x$SNP), pop=pop)
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



#' Harmonise LD matrix against summary data
#'
#' LD matrix returns with rsid_ea_oa identifiers. Make sure that they are oriented to the same effect allele as the summary dataset. Summary dataset can be exposure dataset or harmonised dartaset
#'
#' @param x Exposure dataset or harmonised dataset
#' @param ld Output from ld_matrix
#'
#' @export
#' @return List of exposure dataset and harmonised LD matrix
harmonise_ld_dat <- function(x, ld)
{
	snpnames <- do.call(rbind, strsplit(rownames(ld), split="_"))
	i1 <- snpnames[,1] %in% x$SNP
	ld <- ld[i1,i1]
	snpnames <- snpnames[i1,]
	i2 <- x$SNP %in% snpnames[,1]
	x <- x[i2,]
	stopifnot(all(snpnames[,1] == x$SNP))
	x$effect_allele.exposure <- as.character(x$effect_allele.exposure)
	x$other_allele.exposure <- as.character(x$other_allele.exposure)
	# Set1 x and ld alleles match
	snpnames <- data.frame(snpnames, stringsAsFactors=FALSE)
	snpnames <- merge(subset(x, select=c(SNP, effect_allele.exposure, other_allele.exposure)), snpnames, by.x="SNP", by.y="X1")
	snpnames <- snpnames[match(x$SNP, snpnames$SNP),]
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
		message(" - the following SNPs could not be aligned to the LD reference panel: \n- ", paste(subset(snpnames, !keep)$SNP, collapse="\n - "))
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

	if(any(!snpnames$keep))
	{
		message("Removing ", sum(!snpnames$keep), " variants due to harmonisation issues")
		ld <- ld[snpnames$keep, snpnames$keep]
		x <- x[snpnames$keep, ]

	}
	return(list(x=x, ld=ld))
}



#' Wrapper for MR-PRESSO
#'
#' See <https://github.com/rondolab/MR-PRESSO> for more details.
#'
#' @param dat Output from [`harmonise_data`].
#' @param NbDistribution Number of bootstrap replications. The default is `1000`.
#' @param SignifThreshold Outlier significance threshold. The default is `0.05`.
#'
#' @export
#' @return List of results for every exposure/outcome combination
run_mr_presso <- function(dat, NbDistribution = 1000,  SignifThreshold = 0.05)
{
	dat <- subset(dat, mr_keep)
	d <- subset(dat, !duplicated(paste(id.exposure, " - ", id.outcome)), select=c(exposure, outcome, id.exposure, id.outcome))
	res <- list()
	attributes(res)$id.exposure <- d$id.exposure
	attributes(res)$id.outcome <- d$id.outcome
	attributes(res)$exposure <- d$exposure
	attributes(res)$outcome <- d$outcome
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
#' Creates a list of RadialMR format datasets for each exposure - outcome pair.
#'
#' @param dat Output from [`harmonise_data`].
#'
#' @export
#' @return List of RadialMR format datasets
dat_to_RadialMR <- function(dat)
{
	out <- plyr::dlply(dat, c("exposure", "outcome"), function(x)
	{
		x <- plyr::mutate(x)
		message("Converting:")
		message(" - exposure: ", x$exposure[1])
		message(" - outcome: ", x$outcome[1])
		d <- subset(x, mr_keep=TRUE)
		d <- RadialMR::format_radial(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, RSID=d$SNP)
		return(d)
	})
	return(out)
}

#' Radial IVW analysis
#' 
#' @md
#' @param b_exp Vector of genetic effects on exposure.
#' @param b_out Vector of genetic effects on outcome.
#' @param se_exp Standard errors of genetic effects on exposure.
#' @param se_out Standard errors of genetic effects on outcome.
#' @param parameters List of parameters.
#'
#' @export
#' @return List with the following elements:
#' \describe{
#' \item{b}{causal effect estimate}
#' \item{se}{standard error}
#' \item{pval}{p-value}
#' }
#' @importFrom stats pchisq pnorm
mr_ivw_radial <- function(b_exp, b_out, se_exp, se_out, parameters=default_parameters())
{
  if (sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) &
		!is.na(se_out)) < 2)
		return(list(b = NA, se = NA, pval = NA, nsnp = NA))
	d <- RadialMR::format_radial(BXG=b_exp, BYG=b_out, seBXG=se_exp, seBYG=se_out, RSID=1:length(b_exp))
	out <- RadialMR::ivw_radial(d, alpha=0.05, weights=3)
	b <- out$coef[1,1]
	se <- out$coef[1,2]
	pval <- 2 * pnorm(abs(b/se), lower.tail = FALSE)
	Q_df <- out$df
	Q <- out$qstatistic
	Q_pval <- pchisq(Q, Q_df, lower.tail=F)
	return(list(b = b, se = se, pval = pval, nsnp = length(b_exp),
        Q = Q, Q_df = Q_df, Q_pval = Q_pval))
}


#' Perform MRMix analysis on harmonised dat object
#'
#' See <https://github.com/gqi/MRMix> for more details.
#'
#' @param dat Output from [`harmonise_data`]. Ensures that no eaf.exposure values are missing.
#'
#' @export
#' @return List of results, with one list item for every exposure/outcome pair in dat object
run_mrmix <- function(dat)
{
	plyr::dlply(dat, c("id.exposure", "id.outcome"), function(x)
	{
		message("Analysing ", x$id.exposure[1], " against ", x$id.outcome[1])
		if(grepl("log odds", x$units.exposure[1]))
		{
			xunit <- "binary"
		} else if(grepl("^SD", x$units.exposure[1]))
		{
			xunit <- "n"
		} else {
			xunit <- "continuous"
		}
		if(grepl("log odds", x$units.outcome[1]))
		{
			yunit <- "binary"
		} else if(grepl("^SD", x$units.outcome[1]))
		{
			yunit <- "n"
		} else {
			yunit <- "continuous"
		}
		index <- is.na(x$eaf.exposure)
		if(any(index))
		{
			warning(paste0(x$id.exposure[1], ".", x$id.outcome[1], " - Some eaf.exposure values are missing, using eaf.outcome in their place"))
			x$eaf.exposure[index] <- x$eaf.outcome[index]
		}
		l <- MRMix::standardize(
			x$beta.exposure, 
			x$beta.outcome, 
			x$se.exposure, 
			x$se.outcome, 
			xunit,
			yunit,
			x$samplesize.exposure,
			x$samplesize.outcome,
			x$eaf.exposure
		)
		x$beta.exposure <- l$betahat_x_std
		x$beta.outcome <- l$betahat_y_std
		x$se.exposure <- l$sx_std
		x$se.outcome <- l$sy_std
		MRMix::MRMix(x$beta.exposure, x$beta.outcome, x$se.exposure, x$se.outcome)
	})
}
