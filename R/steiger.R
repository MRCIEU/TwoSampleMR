

#' Calculate p-value from rsq and sample size
#'
#' @param r2 Rsq
#' @param n Sample size
#'
#' @export
#' @return P-value
get_p_from_r2n <- function(r2, n)
{
	fval <- r2 * (n-2) / (1 - r2)
	pval <- pf(fval, 1, n-1, low=FALSE)
	return(pval)
}


#' Calculate variance explained from p vals and sample size
#'
#' @param p Array of pvals
#' @param n Array of sample sizes
#'
#' @export
#' @return r value
get_r_from_pn <- function(p, n)
{
	optim.get_p_from_rn <- function(x, sample_size, pvalue)
	{
		abs(-log10(get_p_from_r2n(x, sample_size)) - -log10(pvalue))
	}

	if(length(p) > 1 & length(n) == 1)
	{
		message("Assuming n the same for all p values")
		n <- rep(n, length(p))
	}

	Fval <- suppressWarnings(qf(p, 1, n-1, low=FALSE))
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
#' @return List
#' - vz: Total volume of the error parameter space
#' - vz0: Volume of the parameter space that gives the incorrect answer
#' - vz1: Volume of the paramtere space that gives the correct answer
#' - sensitivity_ratio: Ratio of vz1/vz0. Higher means inferred direction is less susceptible to measurement error
#' - pl: plot of parameter sapce
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
	temp <- lattice::wireframe(
		z ~ rxx_o * ryy_o, 
		groups=type, 
		data=d, 
		scales=list(arrows=FALSE), 
		col.groups = colorRampPalette(c("red", "blue"))(2), 
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
#' @param r_xxo Measurememt precision of exposure
#' @param r_yyo Measurement precision of outcome
#' @param ... Further arguments to be passed to wireframe
#'
#' @export
#' @return List
#' - r2_exp: Estimated variance explained in x
#' - r2_out: Estimated variance explained in y
#' - r2_exp_adj: Predicted variance explained in x accounting for estimated measurement error
#' - r2_out_adj: Predicted variance explained in y accounting for estimated measurement error
#' - correct_causal_direction: TRUE/FALSE
#' - steiger_test: p-value for inference of direction
#' - correct_causal_direction_adj: TRUE/FALSE, direction of causality for given measurement error parameters
#' - steiger_test_adj: p-value for inference of direction of causality for given measurement error parameters
#' - vz: Total volume of the error parameter space
#' - vz0: Volume of the parameter space that gives the incorrect answer
#' - vz1: Volume of the paramtere space that gives the correct answer
#' - sensitivity_ratio: Ratio of vz1/vz0. Higher means inferred direction is less susceptible to measurement error
#' - sensitivity_plot: Plot of parameter space of causal directions and measurement error
mr_steiger <- function(p_exp, p_out, n_exp, n_out, r_xxo = 1, r_yyo=1, ...)
{
	requireNamespace("psych", quietly=TRUE)
	index <- any(is.na(p_exp)) | any(is.na(p_out)) | any(is.na(n_exp)) | any(is.na(n_out))
	p_exp <- p_exp[!index]
	p_out <- p_out[!index]
	n_exp <- n_exp[!index]
	n_out <- n_out[!index]

	r_exp <- sqrt(sum(get_r_from_pn(p_exp, n_exp)^2))
	r_out <- sqrt(sum(get_r_from_pn(p_out, n_out)^2))

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
		steiger_test = rtest$p,
		correct_causal_direction_adj = r_exp_adj > r_out_adj, 
		steiger_test_adj = rtest_adj$p,
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
#' A statistical test for whether the assumption that exposure causes outcome is valid
#'
#' @param dat Harmonised exposure and outcome data. Output from \code{harmonise_exposure_outcome}
#'
#' @export
#' @return List
#' - 
directionality_test <- function(dat)
{
	if(! all(c("pval.exposure", "pval.outcome", "samplesize.exposure", "samplesize.outcome") %in% names(dat)))
	{
		message("Data requires p-values and sample sizes for outcomes and exposures")
		return(NULL)
	}
	dtest <- plyr::ddply(dat, c("id.exposure", "id.outcome"), function(x)
	{
		b <- mr_steiger(x$pval.exposure, x$pval.outcome, x$samplesize.exposure, x$samplesize.outcome)
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


get_r_from_pn_less_accurate <- function(p, n)
{
	# qval <- qf(p, 1, n-2, low=FALSE)
	p[p == 1] <- 0.999
	p[p == 0] <- 1e-200
	qval <- qchisq(p, 1, low=F) / (qchisq(p, n-2, low=F)/(n-2))
	r <- sqrt(sum(qval / (n - qval)))

	if(r >= 1)
	{
		warning("Correlation greater than 1, make sure SNPs are pruned for LD.")
	}
	return(r)
}


test_r_from_pn <- function()
{
	require(ggplot2)
	require(tidyr)

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

	p <- ggplot(param, aes(x=rsq_emp, value)) +
	geom_abline(slope=1, linetype="dotted") +
	geom_line(aes(colour=out)) +
	facet_grid(. ~ n) +
	scale_x_log10() +
	scale_y_log10()
	return(list(dat=param, p=p))
}
