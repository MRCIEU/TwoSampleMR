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


	if(!"rsq.exposure" %in% names(dat))
	{
		dat$pval.exposure[dat$pval.exposure < 1e-300] <- 1e-300
		if(compareNA(dat$units.exposure[1], "log odds"))
		{
			# message("Estimating rsq.exposure for binary trait")
			# message("Ensure that beta.exposure, eaf.exposure, ncase.exposure, ncontrol.exposure are all specified with no missing values")
			if(! "prevalence.exposure" %in% names(dat))
			{
				dat$prevalence.exposure <- 0.1
				warning("Assuming exposure prevalence of 0.1. Alternatively, add prevalence.exposure column and re-run.")
			}
			ind1 <- !is.na(dat$beta.exposure) &
				!is.na(dat$eaf.exposure) &
				!is.na(dat$ncase.exposure) &
				!is.na(dat$ncontrol.exposure) &
				!is.na(dat$prevalence.exposure)
			dat$rsq.exposure <- NA
			if(sum(ind1) > 0)
			{
				dat$rsq.exposure[ind1] <- get_r_from_lor(
					dat$beta.exposure[ind1],
					dat$eaf.exposure[ind1],
					dat$ncase.exposure[ind1],
					dat$ncontrol.exposure[ind1],
					dat$prevalence.exposure
				)^2
			}
		} else if(all(grepl("SD", dat$units.exposure)) & all(!is.na(dat$eaf.exposure))) {
			dat$rsq.exposure <- NA
			dat$rsq.exposure <- 2 * dat$beta.exposure^2 * dat$eaf.exposure * (1-dat$eaf.exposure)
		} else {
			ind1 <- !is.na(dat$pval.exposure) & !is.na(dat$samplesize.exposure)
			dat$rsq.exposure <- NA
			if(sum(ind1) > 0)
			{		
				dat$rsq.exposure[ind1] <- get_r_from_pn(
					dat$pval.exposure[ind1],
					dat$samplesize.exposure[ind1]
				)^2
			}
		}
	}

	if(!"rsq.outcome" %in% names(dat))
	{
		dat$pval.outcome[dat$pval.outcome < 1e-300] <- 1e-300
		if(compareNA(dat$units.outcome[1], "log odds"))
		{
			if(! "prevalence.outcome" %in% names(dat))
			{
				dat$prevalence.outcome <- 0.1
				warning("Assuming outcome prevalence of 0.1. Alternatively, add prevalence.outcome column and re-run.")
			}
			ind1 <- !is.na(dat$beta.outcome) &
				!is.na(dat$eaf.outcome) &
				!is.na(dat$ncase.outcome) &
				!is.na(dat$ncontrol.outcome) &
				!is.na(dat$prevalence.outcome)
			dat$rsq.outcome <- NA
			if(sum(ind1) > 0)
			{
				dat$rsq.outcome[ind1] <- get_r_from_lor(
					dat$beta.outcome[ind1],
					dat$eaf.outcome[ind1],
					dat$ncase.outcome[ind1],
					dat$ncontrol.outcome[ind1],
					dat$prevalence.outcome
				)^2
			}
		} else if(all(grepl("SD", dat$units.outcome)) & all(!is.na(dat$eaf.outcome))) {
			dat$rsq.outcome <- NA
			dat$rsq.outcome <- 2 * dat$beta.outcome^2 * dat$eaf.outcome * (1-dat$eaf.outcome)
		} else {
			ind1 <- !is.na(dat$pval.outcome) & !is.na(dat$samplesize.outcome)
			dat$rsq.outcome <- NA
			if(sum(ind1) > 0)
			{		
				dat$rsq.outcome[ind1] <- get_r_from_pn(
					dat$pval.outcome[ind1],
					dat$samplesize.outcome[ind1]
				)^2
			}
		}
	}
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