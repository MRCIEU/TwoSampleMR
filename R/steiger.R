

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
