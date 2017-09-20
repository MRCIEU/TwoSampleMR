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
		if(get_correlations)
		{
			ld <- ld_matrix(unique(x$SNP))
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
