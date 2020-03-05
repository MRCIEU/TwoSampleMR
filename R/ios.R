#' Calculate index of suspicion
#'
#' If a SNP influences multiple other traits then it could be 'suspicious', and more likely to be pleiotropic. This function implements two basic approaches to estimate IOS
#'
#' - ios1: A summary of the SNP r2 with the other traits (r2_gu)
#' - ios2: A summary of the ratio of r2_gu / r2_gx, where r2_gx is the variance explained by the SNP on the exposure. Estimates the index of suspicion, whereupon SNPs which have a larger effect on a set of traits given their effect on the exposure are deemed more suspicious
#'
#' Summarising across multiple traits can be dune using mean, sd, iqr, median, 95% value, maximum value
#'
#' @param exposure_dat Instruments for the exposure, obtained using \code{extract_instruments}
#' @param background_dat Effects for the instruments on a set of variables, used to calculate index of suspicion
#'
#' @export
#' @return Data.frame
#' @importFrom stats median quantile sd
ios <- function(exposure_dat, background_dat)
{
	require(dplyr)
	require(reshape2)
	background_dat$vgu <- background_dat$beta.outcome^2 * 2 * background_dat$eaf.outcome * (1 - background_dat$eaf.outcome)
	exposure_dat$vgx <- exposure_dat$beta.exposure^2 * 2 * exposure_dat$eaf.exposure * (1 - exposure_dat$eaf.exposure)
	background_dat <- merge(background_dat, subset(exposure_dat, select=c(SNP, vgx)), by="SNP")
	background_dat$r2_ratio <- background_dat$vgu / background_dat$vgx
	ios <- dplyr::group_by(background_dat, SNP) %>%
		dplyr::summarise(
			ios1_mean = sum(vgu, na.rm=TRUE),
			ios1_sd = sd(vgu, na.rm=TRUE),
			ios1_iqr = quantile(vgu, 0.75, na.rm=TRUE) - quantile(vgu, 0.25, na.rm=TRUE),
			ios1_median = median(vgu, na.rm=TRUE),
			ios1_95 = quantile(vgu, 0.95, na.rm=TRUE),
			ios1_max = max(vgu, na.rm=TRUE),
			ios2_mean = sum(r2_ratio, na.rm=TRUE),
			ios2_sd = sd(r2_ratio, na.rm=TRUE),
			ios2_iqr = quantile(r2_ratio, 0.75, na.rm=TRUE) - quantile(r2_ratio, 0.25, na.rm=TRUE),
			ios2_median = median(r2_ratio, na.rm=TRUE),
			ios2_95 = quantile(r2_ratio, 0.95, na.rm=TRUE),
			ios2_max = max(r2_ratio, na.rm=TRUE)
		) %>% melt
	return(ios)
}

