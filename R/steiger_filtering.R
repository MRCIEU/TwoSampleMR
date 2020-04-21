#' Steiger filtering function
#' 
#' This function takes an object from [`harmonise_data`] and does the following:
#' If there is no rsq.exposure or rsq.outcome column it will try to estimate it. This is done differently for traits that have "log odds" units. 
#' To estimate rsq for quantitative traits there must be either p-values and sample sizes for each SNP, 
#' or effect sizes and standard errors AND the units are in SD units (the column must contain "SD"). 
#' To estimate rsq for binary traits the units must be called "log odds" and there must be beta.exposure, 
#' eaf.exposure, ncase.exposure, ncontrol.exposure, prevalence.exposure. 
#' The same principles apply for calculating the rsq for the outcome trait, except column names are beta.outcome etc. 
#' If prevalence isn't supplied then it uses 0.1 by default.
#' 
#' Once rsq is calculated for the exposure and outcome, it will then perform the Steiger test for each SNP to see if the rsq of the exposure is larger than the rsq of the outcome. 
#' 
#' Note that Steiger filtering, while useful, does have its own pitfalls. 
#' Try to use replication effect estimates for the exposure (which are not biased by winner's curse), 
#' and note that if there is strong antagonistic confounding or differential measurement error between the exposure and outcome then the causal directions could be inferred incorrectly.
#'
#' @param dat Output from [`harmonise_data`].
#' 
#' @export
#' @return [`harmonise_data`] style data frame with additional columns rsq.exposure, rsq.outcome, steiger_dir (which is `TRUE` if the rsq.exposure is larger than rsq.outcome) and steiger_pval which is a test to see if rsq.exposure is significantly larger than rsq.outcome.
steiger_filtering <- function(dat)
{
	plyr::ddply(dat, c("id.exposure", "id.outcome"), steiger_filtering_internal)
}



#' @importFrom stats pnorm
steiger_filtering_internal <- function(dat)
{
	if(! "units.outcome" %in% names(dat))
	{
		dat$units.outcome <- NA
	}
	if(! "units.exposure" %in% names(dat))
	{
		dat$units.exposure <- NA
	}
	stopifnot(length(unique(dat$exposure)) == 1)
	stopifnot(length(unique(dat$outcome)) == 1)
	stopifnot(length(unique(dat$units.exposure)) == 1)
	stopifnot(length(unique(dat$units.outcome)) == 1)

	dat <- add_rsq(dat)

	st <- psych::r.test(
		n = dat$samplesize.exposure, 
		n2 = dat$samplesize.outcome, 
		r12 = sqrt(dat$rsq.exposure), 
		r34 = sqrt(dat$rsq.outcome)
	)
	dat$steiger_dir <- dat$rsq.exposure > dat$rsq.outcome
	dat$steiger_pval <- pnorm(-abs(st$z)) * 2

	return(dat)
}


compareNA <- function(v1,v2) {
    same <- (v1 == v2) | (is.na(v1) & is.na(v2))
    same[is.na(same)] <- FALSE
    return(same)
}


#' Estimate r-square of each association
#'
#' Can be applied to exposure_dat, outcome_dat or harmonised_data. Note that it will be beneficial in some circumstances to add the meta data to the data object using \code{add_metadata()} before running this function.
#'
#' @param dat exposure_dat, outcome_dat or harmonised_data
#'
#' @export
#' @return data frame
add_rsq <- function(dat)
{
	if("id.exposure" %in% names(dat))
	{
		dat <- add_rsq_one(dat, "exposure")
	}
	if("id.outcome" %in% names(dat))
	{
		dat <- add_rsq_one(dat, "outcome")
	}
	return(dat)
}

add_rsq_one <- function(dat, what="exposure")
{
	if(! paste0("units.", what) %in% names(dat))
	{
		dat[[paste0("units.", what)]] <- NA
	}
	stopifnot(length(unique(dat[[what]])) == 1)
	stopifnot(length(unique(dat[[paste0("units.", what)]])) == 1)

	if(! paste0("rsq.", what) %in% names(dat))
	{
		dat[[paste0("pval.", what)]][dat[[paste0("pval.", what)]] < 1e-300] <- 1e-300
		if(compareNA(dat[[paste0("units.", what)]][1], "log odds"))
		{
			# message("Estimating rsq.exposure for binary trait")
			# message("Ensure that beta.exposure, eaf.exposure, ncase.exposure, ncontrol.exposure are all specified with no missing values")
			if(! paste0("prevalence.", what) %in% names(dat))
			{
				dat[[paste0("prevalence.", what)]] <- 0.1
				warning(paste0("Assuming ", what, " prevalence of 0.1. Alternatively, add prevalence.", what, " column and re-run."))
			}
			ind1 <- !is.na(dat[[paste0("beta.", what)]]) &
				!is.na(dat[[paste0("eaf.", what)]]) &
				!is.na(dat[[paste0("ncase.", what)]]) &
				!is.na(dat[[paste0("ncontrol.", what)]]) &
				!is.na(dat[[paste0("prevalence.", what)]])
			dat[[paste0("rsq.", what)]] <- NA
			if(sum(ind1) > 0)
			{
				dat[[paste0("rsq.", what)]][ind1] <- get_r_from_lor(
					dat[[paste0("beta.", what)]][ind1],
					dat[[paste0("eaf.", what)]][ind1],
					dat[[paste0("ncase.", what)]][ind1],
					dat[[paste0("ncontrol.", what)]][ind1],
					dat[[paste0("prevalence.", what)]]
				)^2
			}
		} else if(all(grepl("SD", dat[[paste0("units.", what)]])) & all(!is.na(dat[[paste0("eaf.", what)]]))) {
			dat[[paste0("rsq.", what)]] <- NA
			dat[[paste0("rsq.", what)]] <- 2 * dat[[paste0("beta.", what)]]^2 * dat[[paste0("eaf.", what)]] * (1-dat[[paste0("eaf.", what)]])
		} else {
			ind1 <- !is.na(dat[[paste0("pval.", what)]]) & !is.na(dat[[paste0("samplesize.", what)]])
			dat[[paste0("rsq.", what)]] <- NA
			if(sum(ind1) > 0)
			{		
				dat[[paste0("rsq.", what)]][ind1] <- get_r_from_pn(
					dat[[paste0("pval.", what)]][ind1],
					dat[[paste0("samplesize.", what)]][ind1]
				)^2
			}
		}
	}
	return(dat)
}
