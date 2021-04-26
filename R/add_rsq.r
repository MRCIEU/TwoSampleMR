#' Estimate r-square of each association
#'
#' Can be applied to exposure_dat, outcome_dat or harmonised_data. Note that it will be beneficial in some circumstances to add the meta data to the data object using \code{add_metadata()} before running this function. Also adds effective sample size for case control data
#'
#' @param dat exposure_dat, outcome_dat or harmonised_data
#'
#' @export
#' @return data frame
add_rsq <- function(dat)
{
	if("id.exposure" %in% names(dat))
	{
		dat <- 	plyr::ddply(dat, c("id.exposure"), function(x)
			{
				add_rsq_one(x, "exposure")
			})
	}
	if("id.outcome" %in% names(dat))
	{
		dat <- 	plyr::ddply(dat, c("id.outcome"), function(x)
			{
				add_rsq_one(x, "outcome")
			})
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
				dat[[paste0("effective_n.", what)]][ind1] <- effective_n(dat[[paste0("ncase.", what)]][ind1], dat[[paste0("ncontrol.", what)]][ind1])
			}
		} else if(all(grepl("SD", dat[[paste0("units.", what)]])) & all(!is.na(dat[[paste0("eaf.", what)]]))) {
			dat[[paste0("rsq.", what)]] <- NA
			dat[[paste0("rsq.", what)]] <- 2 * dat[[paste0("beta.", what)]]^2 * dat[[paste0("eaf.", what)]] * (1-dat[[paste0("eaf.", what)]])
			dat[[paste0("effective_n.", what)]] <- dat[[paste0("samplesize.", what)]]
		} else {
			ind1 <- !is.na(dat[[paste0("pval.", what)]]) & !is.na(dat[[paste0("samplesize.", what)]])
			dat[[paste0("rsq.", what)]] <- NA
			if(sum(ind1) > 0)
			{		
				dat[[paste0("rsq.", what)]][ind1] <- get_r_from_pn(
					dat[[paste0("pval.", what)]][ind1],
					dat[[paste0("samplesize.", what)]][ind1]
				)^2
				dat[[paste0("effective_n.", what)]] <- dat[[paste0("samplesize.", what)]]
			}
		}
	}
	return(dat)
}


#' @importFrom stats qchisq
get_r_from_pn_less_accurate <- function(p, n)
{
	# qval <- qf(p, 1, n-2, low=FALSE)
	p[p == 1] <- 0.999
	p[p == 0] <- 1e-200
	qval <- qchisq(p, 1, lower.tail = FALSE) / (qchisq(p, n-2, lower.tail = FALSE)/(n-2))
	r <- sqrt(sum(qval / (n - qval)))

	if(r >= 1)
	{
		warning("Correlation greater than 1, make sure SNPs are pruned for LD.")
	}
	return(r)
}

#' @importFrom stats coefficients cor lm rnorm
test_r_from_pn <- function()
{
	requireNamespace("ggplot2", quietly = TRUE)
	# requireNamespace("tidyr", quietly = TRUE) # not used

	param <- expand.grid(
		n = c(10, 100, 1000, 10000, 100000),
		rsq = 10^seq(-4,-0.5, length.out=30)
	)

	for(i in 1:nrow(param))
	{
		message(i)
		x <- scale(rnorm(param$n[i]))
		y <- x * sqrt(param$rsq[i]) + scale(rnorm(param$n[i])) * sqrt(1 - param$rsq[i])
		param$rsq_emp[i] <- cor(x, y)^2
		param$pval[i] <- max(coefficients(summary(lm(y ~ x)))[2,4], 1e-300)
		param$rsq1[i] <- get_r_from_pn_less_accurate(param$pval[i], param$n[i])^2
		param$rsq2[i] <- get_r_from_pn(param$pval[i], param$n[i])^2
	}

	param <- gather(param, key=out, value=value, rsq1, rsq2)

	p <- ggplot2::ggplot(param, ggplot2::aes(x=rsq_emp, value)) +
	  ggplot2::geom_abline(slope=1, linetype="dotted") +
	  ggplot2::geom_line(aes(colour=out)) +
	  ggplot2::facet_grid(. ~ n) +
	  ggplot2::scale_x_log10() +
	  ggplot2::scale_y_log10()
	return(list(dat=param, p=p))
}


#' Calculate p-value from rsq and sample size
#'
#' @param r2 Rsq
#' @param n Sample size
#'
#' @export
#' @return P-value
#' @importFrom stats pf df
get_p_from_r2n <- function(r2, n)
{
	fval <- r2 * (n-2) / (1 - r2)
	pval <- pf(fval, 1, n-1, lower.tail=FALSE)
	return(pval)
}


#' Calculate variance explained from p vals and sample size
#'
#' @param p Array of pvals
#' @param n Array of sample sizes
#'
#' @export
#' @return r value
#' @importFrom stats optim qf
get_r_from_pn <- function(p, n)
{
	message("Estimating correlation for quantitative trait.")
	message("This method is an approximation, and may be numerically unstable.")
	message("Ideally you should estimate r directly from independent replication samples.")
	message("Use get_r_from_lor for binary traits.")
	optim.get_p_from_rn <- function(x, sample_size, pvalue)
	{
		abs(-log10(suppressWarnings(get_p_from_r2n(x, sample_size))) - -log10(pvalue))
	}

	if(length(p) > 1 & length(n) == 1)
	{
		message("Assuming n the same for all p values")
		n <- rep(n, length(p))
	}

	Fval <- suppressWarnings(qf(p, 1, n-1, lower.tail=FALSE))
	R2 <- Fval / (n - 2 + Fval)
	index <- !is.finite(Fval)
	if(any(index))
	{
		index <- which(index)
		for(i in 1:length(index))
		{
			if(p[index[i]] == 0)
			{
				R2[index[i]] <- NA
				warning("P-value of 0 cannot be converted to R value")
			} else {
				R2[index[i]] <- suppressWarnings(optim(0.001, optim.get_p_from_rn, sample_size=n[index[i]], pvalue=p[index[i]])$par)
			}
		}
	}
	return(sqrt(R2))
}

compareNA <- function(v1,v2) {
    same <- (v1 == v2) | (is.na(v1) & is.na(v2))
    same[is.na(same)] <- FALSE
    return(same)
}



#' Estimate proportion of variance of liability explained by SNP in general population
#'
#' This uses equation 10 in Lee et al. A Better Coefficient of Determination for Genetic Profile Analysis. 
#' Genetic Epidemiology 36: 214–224 (2012) <https://doi.org/10.1002/gepi.21614>.
#'
#' @param lor Vector of Log odds ratio.
#' @param af Vector of allele frequencies.
#' @param ncase Vector of Number of cases.
#' @param ncontrol Vector of Number of controls.
#' @param prevalence Vector of Disease prevalence in the population.
#' @param model Is the effect size estimated from the `"logit"` (default) or `"probit"` model.
#' @param correction Scale the estimated r by correction value. The default is `FALSE`.
#'
#' @export
#' @return Vector of r values
get_r_from_lor <- function(lor, af, ncase, ncontrol, prevalence, model="logit", correction=FALSE)
{
	stopifnot(length(lor) == length(af))
	stopifnot(length(ncase) == 1 | length(ncase) == length(lor))
	stopifnot(length(ncontrol) == 1 | length(ncontrol) == length(lor))
	stopifnot(length(prevalence) == 1 | length(prevalence) == length(lor))
	if(length(prevalence) == 1 & length(lor) != 1)
	{
		prevalence <- rep(prevalence, length(lor))
	}
	if(length(ncase) == 1 & length(lor) != 1)
	{
		ncase <- rep(ncase, length(lor))
	}
	if(length(ncontrol) == 1 & length(lor) != 1)
	{
		ncontrol <- rep(ncontrol, length(lor))
	}

	nsnp <- length(lor)
	r <- array(NA, nsnp)
	for(i in 1:nsnp)
	{
		if(model == "logit")
		{
			ve <- pi^2/3
		} else if(model == "probit") {
			ve <- 1
		} else {
			stop("Model must be probit or logit")
		}
		popaf <- get_population_allele_frequency(af[i], ncase[i] / (ncase[i] + ncontrol[i]), exp(lor[i]), prevalence[i])
		vg <- lor[i]^2 * popaf * (1-popaf)
		r[i] <- vg / (vg + ve)
		if(correction)
		{
			r[i] <- r[i] / 0.58
		}
		r[i] <- sqrt(r[i])
	}
	return(r)
}



#' Obtain 2x2 contingency table from marginal parameters and odds ratio
#'
#' Columns are the case and control frequencies.
#' Rows are the frequencies for allele 1 and allele 2.
#'
#' @param af Allele frequency of effect allele.
#' @param prop Proportion of cases.
#' @param odds_ratio Odds ratio.
#' @param eps tolerance, default is `1e-15`.
#'
#' @export
#' @return 2x2 contingency table as matrix
contingency <- function(af, prop, odds_ratio, eps=1e-15)
{
	a <- odds_ratio-1
	b <- (af+prop)*(1-odds_ratio)-1
	c_ <- odds_ratio*af*prop

	if (abs(a) < eps)
	{
		z <- -c_ / b
	} else {
		d <- b^2 - 4*a*c_
		if (d < eps*eps) 
		{
			s <- 0
		} else {
			s <- c(-1,1)
		}
		z <- (-b + s*sqrt(max(0, d))) / (2*a)
	}
	y <- vapply(z, function(a) zapsmall(matrix(c(a, prop-a, af-a, 1+a-af-prop), 2, 2)), matrix(0.0, 2, 2))
	i <- apply(y, 3, function(u) all(u >= 0))
	return(y[,,i])
}

#' Estimate allele frequency from SNP
#'
#' @param g Vector of 0/1/2
#'
#' @export
#' @return Allele frequency 
allele_frequency <- function(g)
{
	(sum(g == 1) + 2 * sum(g == 2)) / (2 * sum(!is.na(g)))
}


#' Estimate the allele frequency in population from case/control summary data
#'
#' @param af Effect allele frequency (or MAF)
#' @param prop Proportion of samples that are cases
#' @param odds_ratio Odds ratio
#' @param prevalence Population disease prevalence
#'
#' @export
#' @return Population allele frequency
get_population_allele_frequency <- function(af, prop, odds_ratio, prevalence)
{
	stopifnot(length(af) == length(odds_ratio))
	stopifnot(length(prop) == length(odds_ratio))
	for(i in 1:length(odds_ratio))
	{
		co <- contingency(af[i], prop[i], odds_ratio[i])
		af_controls <- co[1,2] / (co[1,2] + co[2,2])
		af_cases <- co[1,1] / (co[1,1] + co[2,1])
		af[i] <- af_controls * (1 - prevalence) + af_cases * prevalence
	}
	return(af)
}


#' Estimate the effective sample size in a case control study
#'
#' Taken from https://www.nature.com/articles/nprot.2014.071
#'
#' @param ncase Vector of number of cases
#' @param ncontrol Vector of number of controls
#'
#' @export
#' @return Vector of effective sample size
effective_n <- function(ncase, ncontrol)
{
	return(2 / (1/ncase + 1/ncontrol))
}
