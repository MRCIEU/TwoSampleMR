#' Try to standardise continuous traits to be in standard deviation units
#'
#' Uses [`estimate_trait_sd`].
#'
#' @param dat Output from [`harmonise_data`].
#'
#' @export
#' @return Data frame
standardise_units <- function(dat)
{
	out <- plyr::ddply(dat, c("id.exposure", "id.outcome"), function(d)
	{
		if(d$units.exposure[1] != "log odds")
		{
			estsd <- mean(estimate_trait_sd(d$beta.exposure, d$se.exposure, d$samplesize.exposure, d$eaf.exposure), na.rm=TRUE)

			if(!is.na(estsd))
			{
				d$beta.exposure <- d$beta.exposure / estsd
				d$se.exposure <- d$se.exposure / estsd
				d$units.exposure <- "SD"
				d$estimated_sd.exposure <- estsd
			}
		}
		if(d$units.outcome[1] != "log odds")
		{
			estsd <- mean(estimate_trait_sd(d$beta.outcome, d$se.outcome, d$samplesize.outcome, d$eaf.outcome), na.rm=TRUE)
			if(!is.na(estsd))
			{
				d$beta.outcome <- d$beta.outcome / estsd
				d$se.outcome <- d$se.outcome / estsd
				d$units.outcome <- "SD"
				d$estimated_sd.outcome <- estsd
			}
		}
		return(d)
	})
	return(out)
}


#' Estimate trait SD by obtaining beta estimates from z-scores and finding the ratio with original beta values
#'
#' Assumes that sample size and allele frequency is correct for each SNP, and that allele frequency gives a reasonable estimate of the variance of the SNP.
#'
#' @param b vector of effect sizes.
#' @param se vector of standard errors.
#' @param n vector of sample sizes.
#' @param p vector of allele frequencies.
#'
#' @export
#' @return Vector of sd estimates for each association.
estimate_trait_sd <- function(b, se, n, p)
{
	z <- b / se
	standardised_bhat <- sqrt((z^2/(z^2+n-2)) / (2 * p * (1-p))) * sign(z)
	estimated_sd <- b / standardised_bhat
	return(estimated_sd)
}
