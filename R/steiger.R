

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


#' Evaluate the Steiger test's sensitivity to measurement error
#'
#' @param rgx_o Observed variance of exposure explained by SNPs
#' @param rgy_o Observed variance of outcome explained by SNPs
#' @param ... Further arguments to be passed to wireframe
#'
#' @export
#' @return List with the following elements:
#' \describe{
#' \item{vz}{Total volume of the error parameter space}
#' \item{vz0}{Volume of the parameter space that gives the incorrect answer}
#' \item{vz1}{Volume of the paramtere space that gives the correct answer}
#' \item{sensitivity_ratio}{Ratio of vz1/vz0. Higher means inferred direction is less susceptible to measurement error}
#' \item{pl}{plot of parameter space}
#' }
#' @importFrom grDevices rgb
#' @importFrom lattice wireframe
steiger_sensitivity <- function(rgx_o, rgy_o, ...)
{
	requireNamespace("lattice", quietly=TRUE)
	if(rgy_o > rgx_o)
	{
		a <- rgy_o
		b <- rgx_o
	} else {
		a <- rgx_o
		b <- rgy_o
	}

	d <- expand.grid(rxx_o=seq(rgx_o,1,length.out=50), ryy_o=seq(rgy_o,1,length.out=50), type=c("A","B"))
	d$rgy <- rgy_o / d$ryy_o
	d$rgx <- rgx_o / d$rxx_o
	d$z <- d$rgy - d$rgx
	d$z[d$type=="A"] <- 0
	mycolors.trans = rgb(c(255,0), c(0,0), 
               c(0,255),alpha = c(70,255), maxColorValue = 255) 

	temp <- lattice::wireframe(
		z ~ rxx_o * ryy_o, 
		groups=type, 
		data=d, 
		scales=list(arrows=FALSE), 
		col.groups = mycolors.trans, 
		drape=FALSE, 
		ylab=expression(rho[xx[o]]), 
		xlab=expression(rho[yy[o]]),
		zlab=expression(rho[gy]-rho[gx]),
		par.settings = list(axis.line=list(col="transparent")),
		...
	)

	vz <- a * log(a) - b * log(b) + a*b*(log(b)-log(a))
	vz0 <- -2*b - b * log(a) - a*b*log(a) + 2*a*b

	vz1 <- abs(vz - vz0)

	sensitivity <- vz0 / (2 * vz0 + abs(vz))
	sensitivity_ratio <- vz1 / vz0

	return(list(
		vz = vz,
		vz0 = vz0,
		vz1 = vz1,
		# sensitivity = sensitivity,
		sensitivity_ratio = sensitivity_ratio,
		pl = temp
	))
}


#' MR Steiger test of directionality
#'
#' A statistical test for whether the assumption that exposure causes outcome is valid
#'
#' @param p_exp Vector of p-values of SNP-exposure
#' @param p_out Vector of p-values of SNP-outcome
#' @param n_exp Sample sizes for p_exp
#' @param n_out Sample sizes for p_out
#' @param r_exp Vector of absolute correlations for SNP-exposure
#' @param r_out Vector of absolute correlations for SNP-outcome
#' @param r_xxo Measurememt precision of exposure
#' @param r_yyo Measurement precision of outcome
#' @param ... Further arguments to be passed to wireframe
#'
#' @export
#' @return List with the following elements:
#' \describe{
#' \item{r2_exp}{Estimated variance explained in x}
#' \item{r2_out}{Estimated variance explained in y}
#' \item{r2_exp_adj}{Predicted variance explained in x accounting for estimated measurement error}
#' \item{r2_out_adj}{Predicted variance explained in y accounting for estimated measurement error}
#' \item{correct_causal_direction}{TRUE/FALSE}
#' \item{steiger_test}{p-value for inference of direction}
#' \item{correct_causal_direction_adj}{TRUE/FALSE, direction of causality for given measurement error parameters}
#' \item{steiger_test_adj}{p-value for inference of direction of causality for given measurement error parameters}
#' \item{vz}{Total volume of the error parameter space}
#' \item{vz0}{Volume of the parameter space that gives the incorrect answer}
#' \item{vz1}{Volume of the paramtere space that gives the correct answer}
#' \item{sensitivity_ratio}{Ratio of vz1/vz0. Higher means inferred direction is less susceptible to measurement error}
#' \item{sensitivity_plot}{Plot of parameter space of causal directions and measurement error}
#' }
#' @importFrom stats pnorm
mr_steiger <- function(p_exp, p_out, n_exp, n_out, r_exp, r_out, r_xxo = 1, r_yyo=1, ...)
{
	requireNamespace("psych", quietly=TRUE)

	r_exp <- abs(r_exp)
	r_out <- abs(r_out)

	ir_exp <- is.na(r_exp)
	ir_out <- is.na(r_out)

	ip_exp <- is.na(p_exp) | is.na(n_exp)
	ip_out <- is.na(p_out) | is.na(n_out)

	if(any(ir_exp))
	{
		r_exp[ir_exp] <- get_r_from_pn(p_exp[ir_exp & !ip_exp], n_exp[ir_exp & !ip_exp])
	}
	if(any(ir_out))
	{
		r_out[ir_out] <- get_r_from_pn(p_out[ir_out & !ip_out], n_out[ir_out & !ip_out])
	}

	r_exp <- sqrt(sum(r_exp[!is.na(r_exp) | is.na(r_out)]^2))
	r_out <- sqrt(sum(r_out[!is.na(r_exp) | is.na(r_out)]^2))

	stopifnot(r_xxo <= 1 & r_xxo >= 0)
	stopifnot(r_yyo <= 1 & r_yyo >= 0)

	r_exp_adj <- sqrt(r_exp^2 / r_xxo^2)
	r_out_adj <- sqrt(r_out^2 / r_yyo^2)

	sensitivity <- steiger_sensitivity(r_exp, r_out, ...)

	rtest <- psych::r.test(n = mean(n_exp), n2 = mean(n_out), r12 = r_exp, r34 = r_out)
	rtest_adj <- psych::r.test(n = mean(n_exp), n2 = mean(n_out), r12 = r_exp_adj, r34 = r_out_adj)
	l <- list(
		r2_exp = r_exp^2, 
		r2_out = r_out^2, 
		r2_exp_adj = r_exp_adj^2, 
		r2_out_adj = r_out_adj^2, 
		correct_causal_direction = r_exp > r_out, 
		steiger_test = pnorm(-abs(rtest[["z"]])) * 2,
		correct_causal_direction_adj = r_exp_adj > r_out_adj, 
		steiger_test_adj = pnorm(-abs(rtest_adj[["z"]])) * 2,
		vz = sensitivity$vz,
		vz0 = sensitivity$vz0,
		vz1 = sensitivity$vz1,
		sensitivity_ratio = sensitivity$sensitivity_ratio,
		sensitivity_plot = sensitivity$pl
	)
	return(l)
}


#' Perform MR Steiger test of directionality
#'
#' A statistical test for whether the assumption that exposure causes outcome is valid.
#'
#' @param dat Harmonised exposure and outcome data. Output from \code{\link{harmonise_data}}.
#'
#' @export
#' @return List
#' 
directionality_test <- function(dat)
{
	if(! all(c("r.exposure", "r.outcome") %in% names(dat)))
	{
		message("r.exposure and/or r.outcome not present.")
		if(! all(c("pval.exposure", "pval.outcome", "samplesize.exposure", "samplesize.outcome") %in% names(dat)))
		{
			message("Can't calculate approximate SNP-exposure and SNP-outcome correlations without pval.exposure, pval.outcome, samplesize.exposure, samplesize.outcome")
			message("Either supply these values, or supply the r.exposure and r.outcome values")
			message("Note, automated correlations assume quantitative traits. For binary traits please pre-calculate in r.exposure and r.outcome e.g. using get_r_from_lor()")
			return(NULL)
		} else {
			message("Calculating approximate SNP-exposure and/or SNP-outcome correlations, assuming all are quantitative traits. Please pre-calculate r.exposure and/or r.outcome using get_r_from_lor() for any binary traits")
		}
	}
	dtest <- plyr::ddply(dat, c("id.exposure", "id.outcome"), function(x)
	{
		if(!"r.exposure" %in% names(x)) 
		{
			x$r.exposure <- NA
		}
		if(!"r.outcome" %in% names(x)) 
		{
			x$r.outcome <- NA
		}
		b <- mr_steiger(x$pval.exposure, x$pval.outcome, x$samplesize.exposure, x$samplesize.outcome, x$r.exposure, x$r.outcome)
		a <- data.frame(
			exposure = x$exposure[1],
			outcome = x$outcome[1],
			snp_r2.exposure = b$r2_exp,
			snp_r2.outcome = b$r2_out,
			correct_causal_direction = b$correct_causal_direction,
			steiger_pval = b$steiger_test
		)
		return(a)
	})

	return(dtest)
}

#' @importFrom stats qchisq
get_r_from_pn_less_accurate <- function(p, n)
{
	# qval <- qf(p, 1, n-2, low=FALSE)
	p[p == 1] <- 0.999
	p[p == 0] <- 1e-200
	qval <- qchisq(p, 1, lower.tail=F) / (qchisq(p, n-2, lower.tail=F)/(n-2))
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




#' MR Steiger test of directionality
#'
#' A statistical test for whether the assumption that exposure causes outcome is valid
#'
#' @param r_exp Vector of correlations of SNP-exposure
#' @param r_out Vector of correlations of SNP-outcome
#' @param n_exp Sample sizes for p_exp
#' @param n_out Sample sizes for p_out
#' @param r_xxo Measurememt precision of exposure
#' @param r_yyo Measurement precision of outcome
#' @param ... Further arguments to be passed to wireframe
#'
#' @export
#' @return List with the following elements:
#' \describe{
#' \item{r2_exp}{Estimated variance explained in x}
#' \item{r2_out}{Estimated variance explained in y}
#' \item{r2_exp_adj}{Predicted variance explained in x accounting for estimated measurement error}
#' \item{r2_out_adj}{Predicted variance explained in y accounting for estimated measurement error}
#' \item{correct_causal_direction}{TRUE/FALSE}
#' \item{steiger_test}{p-value for inference of direction}
#' \item{correct_causal_direction_adj}{TRUE/FALSE, direction of causality for given measurement error parameters}
#' \item{steiger_test_adj}{p-value for inference of direction of causality for given measurement error parameters}
#' \item{vz}{Total volume of the error parameter space}
#' \item{vz0}{Volume of the parameter space that gives the incorrect answer}
#' \item{vz1}{Volume of the paramtere space that gives the correct answer}
#' \item{sensitivity_ratio}{Ratio of vz1/vz0. Higher means inferred direction is less susceptible to measurement error}
#' \item{sensitivity_plot}{Plot of parameter space of causal directions and measurement error}
#' }
#' @importFrom stats pnorm rnorm
mr_steiger2 <- function(r_exp, r_out, n_exp, n_out, r_xxo = 1, r_yyo=1, ...)
{
	requireNamespace("psych", quietly=TRUE)
	index <- any(is.na(r_exp)) | any(is.na(r_out)) | any(is.na(n_exp)) | any(is.na(n_out))
	n_exp <- n_exp[!index]
	n_out <- n_out[!index]

	r_exp <- sqrt(sum(r_exp^2))
	r_out <- sqrt(sum(r_out^2))

	stopifnot(r_xxo <= 1 & r_xxo >= 0)
	stopifnot(r_yyo <= 1 & r_yyo >= 0)

	r_exp_adj <- sqrt(r_exp^2 / r_xxo^2)
	r_out_adj <- sqrt(r_out^2 / r_yyo^2)

	sensitivity <- steiger_sensitivity(r_exp, r_out, ...)

	rtest <- psych::r.test(n = mean(n_exp), n2 = mean(n_out), r12 = r_exp, r34 = r_out)
	rtest_adj <- psych::r.test(n = mean(n_exp), n2 = mean(n_out), r12 = r_exp_adj, r34 = r_out_adj)
	l <- list(
		r2_exp = r_exp^2, 
		r2_out = r_out^2, 
		r2_exp_adj = r_exp_adj^2, 
		r2_out_adj = r_out_adj^2, 
		correct_causal_direction = r_exp > r_out, 
		steiger_test = pnorm(-abs(rtest[["z"]]))*2,
		correct_causal_direction_adj = r_exp_adj > r_out_adj, 
		steiger_test_adj = pnorm(-abs(rtest_adj[["z"]]))*2,
		vz = sensitivity$vz,
		vz0 = sensitivity$vz0,
		vz1 = sensitivity$vz1,
		sensitivity_ratio = sensitivity$sensitivity_ratio,
		sensitivity_plot = sensitivity$pl
	)
	return(l)
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
	co <- contingency(af, prop, odds_ratio)
	af_controls <- co[1,2] / (co[1,2] + co[2,2])
	af_cases <- co[1,1] / (co[1,1] + co[2,1])
	af <- af_controls * (1 - prevalence) + af_cases * prevalence
	return(af)
}


#' Estimate proportion of variance of liability explained by SNP in general population
#'
#' This uses equation 10 in Genetic Epidemiology 36 : 214–224 (2012)
#'
#' @param lor Vector of Log odds ratio
#' @param af Vector of allele frequencies 
#' @param ncase Vector of Number of cases
#' @param ncontrol Vector of Number of controls
#' @param prevalence Vector of Disease prevalence in the population
#' @param model Is the effect size estiamted in "logit" (default) or "probit" model
#' @param correction Scale the estimated r by correction value. Default = FALSE
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
