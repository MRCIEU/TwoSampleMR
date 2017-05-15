
#' Perform all Mendelian randomization tests
#'
#' @param dat Harmonised exposure and outcome data. Output from \code{harmonise_exposure_outcome}
#' @param parameters Parameters to be used for various MR methods. Default is output from \code{dafault_param}.
#' @param method_list List of methods to use in analysis. See \code{mr_method_list()} for details.
#'
#' @export
#' @return List with the following elements:
#'         mr: Table of MR results
#'         extra: Table of extra results
mr <- function(dat, parameters=default_parameters(), method_list=subset(mr_method_list(), use_by_default)$obj)
{
	mr_tab <- plyr::ddply(dat, c("id.exposure", "id.outcome"), function(x1)
	{
		# message("Performing MR analysis of '", x1$id.exposure[1], "' on '", x18WII58$id.outcome[1], "'")
		x <- subset(x1, mr_keep)
		if(nrow(x) == 0)
		{
			message("No SNPs available for MR analysis of '", x1$id.exposure[1], "' on '", x1$id.outcome[1], "'")
			return(NULL)
		} else {
			message("Analysing '", x1$id.exposure[1], "' on '", x1$id.outcome[1], "'")
		}
		res <- lapply(method_list, function(meth)
		{
			get(meth)(x$beta.exposure, x$beta.outcome, x$se.exposure, x$se.outcome, parameters)
		})
		methl <- mr_method_list()
		mr_tab <- data.frame(
			outcome = x$outcome[1],
			exposure = x$exposure[1],
			method = methl$name[match(method_list, methl$obj)],
			nsnp = sapply(res, function(x) x$nsnp),
			b = sapply(res, function(x) x$b),
			se = sapply(res, function(x) x$se),
			pval = sapply(res, function(x) x$pval)
		)
		mr_tab <- subset(mr_tab, !(is.na(b) & is.na(se) & is.na(pval)))
		return(mr_tab)
	})

	return(mr_tab)
}


#' Get list of available MR methods
#'
#' @export
#' @return character vector of method names
mr_method_list <- function()
{
	a <- list(
		list(
			obj="mr_wald_ratio",
			name="Wald ratio",
			PubmedID="",
			Description="",
			use_by_default=TRUE,
			heterogeneity_test=FALSE
		),
		list(
			obj="mr_meta_fixed_simple",
			name="Fixed effects meta analysis (simple SE)",
			PubmedID="",
			Description="",
			use_by_default=FALSE,
			heterogeneity_test=FALSE
		),
		list(
			obj="mr_meta_fixed",
			name="Fixed effects meta analysis (delta method)",
			PubmedID="",
			Description="",
			use_by_default=FALSE,
			heterogeneity_test=TRUE
		),
		list(
			obj="mr_meta_random",
			name="Random effects meta analysis (delta method)",
			PubmedID="",
			Description="",
			use_by_default=FALSE,
			heterogeneity_test=TRUE
		),
		list(
			obj="mr_two_sample_ml",
			name="Maximum likelihood",
			PubmedID="",
			Description="",
			use_by_default=FALSE,
			heterogeneity_test=TRUE
		),
		list(
			obj="mr_egger_regression",
			name="MR Egger",
			PubmedID="26050253",
			Description="",
			use_by_default=TRUE,
			heterogeneity_test=TRUE
		),
		list(
			obj="mr_egger_regression_bootstrap",
			name="MR Egger (bootstrap)",
			PubmedID="26050253",
			Description="",
			use_by_default=FALSE,
			heterogeneity_test=FALSE
		),
		list(
			obj="mr_simple_median",
			name="Simple median",
			PubmedID="",
			Description="",
			use_by_default=FALSE,
			heterogeneity_test=FALSE
		),
		list(
			obj="mr_weighted_median",
			name="Weighted median",
			PubmedID="",
			Description="",
			use_by_default=TRUE,
			heterogeneity_test=FALSE
		),
		list(
			obj="mr_penalised_weighted_median",
			name="Penalised weighted median",
			PubmedID="",
			Description="",
			use_by_default=FALSE,
			heterogeneity_test=FALSE
		),
		list(
			obj="mr_ivw",
			name="Inverse variance weighted",
			PubmedID="",
			Description="",
			use_by_default=TRUE,
			heterogeneity_test=TRUE
		),
		list(
			obj="mr_weighted_mode",
			name="Weighted mode",
			PubmedID="",
			Description="",
			use_by_default=TRUE,
			heterogeneity_test=FALSE
		)
	)
	a <- lapply(a, as.data.frame)
	a <- plyr::rbind.fill(a)
	a <- as.data.frame(lapply(a, as.character), stringsAsFactors=FALSE)
	a$heterogeneity_test <- as.logical(a$heterogeneity_test)
	a$use_by_default <- as.logical(a$use_by_default)
	return(a)
}


#' List of parameters for use with MR functions
#'
#' @export
default_parameters <- function()
{
	list(
		nboot = 1000,
		Cov = 0,
		penk = 20,
		phi = 1,
		alpha = 0.05,
		Qthresh = 0.05
	)
}


#' Perform 2 sample IV using Wald ratio
#'
#' @param b_exp Vector of genetic effects on exposure
#' @param b_out Vector of genetic effects on outcome
#' @param se_exp Standard errors of genetic effects on exposure
#' @param se_out Standard errors of genetic effects on outcome
#'
#' @export
#' @return List with the following elements:
#'         b: causal effect estimate
#'         se: standard error
#'         pval: p-value
mr_wald_ratio <- function(b_exp, b_out, se_exp, se_out, parameters)
{
	if(length(b_exp) > 1)
	{
		return(list(b=NA, se=NA, pval=NA, nsnp=NA))
	}
	b <- b_out / b_exp
	se <- se_out / abs(b_exp)
	# sqrt((segd^2/gp^2) + (gd^2/gp^4)*segp^2 - 2*(gd/gp^3)) #full delta method with cov set to 0
	pval <- pnorm(abs(b) / se, lower.tail=FALSE) * 2
	return(list(b=b, se=se, pval=pval, nsnp=1))
}


#' Perform 2 sample IV using simple standard error
#'
#' @param b_exp Vector of genetic effects on exposure
#' @param b_out Vector of genetic effects on outcome
#' @param se_exp Standard errors of genetic effects on exposure
#' @param se_out Standard errors of genetic effects on outcome
#'
#' @export
#' @return List with the following elements:
#'         b: causal effect estimate
#'         se: standard error
#'         pval: p-value
mr_meta_fixed_simple <- function(b_exp, b_out, se_exp, se_out, parameters)
{
	if(sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 2)
	{
		return(list(b=NA, se=NA, pval=NA, nsnp=NA))
	}
	b <- sum(b_exp*b_out / se_out^2) / sum(b_exp^2/se_out^2)
	se <- sqrt(1 / sum(b_exp^2/se_out^2))
	pval <- 2 * pnorm(abs(b) / se, lower.tail=FALSE)
	return(list(b=b, se=se, pval=pval, nsnp=length(b_exp)))
}



#' Perform 2 sample IV using fixed effects meta analysis and delta method for standard errors
#'
#' @param b_exp Vector of genetic effects on exposure
#' @param b_out Vector of genetic effects on outcome
#' @param se_exp Standard errors of genetic effects on exposure
#' @param se_out Standard errors of genetic effects on outcome
#'
#' @export
#' @return List with the following elements:
#'         b: causal effect estimate
#'         se: standard error
#'         pval: p-value
#'         Q, Q_df, Q_pval: Heterogeneity stats
mr_meta_fixed <- function(b_exp, b_out, se_exp, se_out, parameters)
{
	if(sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 1)
	{
		return(list(b=NA, se=NA, pval=NA, nsnp=NA, Q =NA, Q_df =NA, Q_pval =NA))
	}
	ratio <- b_out / b_exp
	ratio.se <- sqrt((se_out^2/b_exp^2) + (b_out^2/b_exp^4)*se_exp^2 - 2*(b_out/b_exp^3)*parameters$Cov)
	res <- meta::metagen(ratio, ratio.se)
	b <- res$TE.fixed
	se <- res$seTE.fixed
	pval <- res$pval.fixed
	Q_pval <- pchisq(res$Q, res$df.Q, low=FALSE)
	return(list(b=b, se=se, pval=pval, nsnp=length(b_exp), Q = res$Q, Q_df = res$df.Q, Q_pval = Q_pval))
}



#' Perform 2 sample IV using random effects meta analysis and delta method for standard errors
#'
#' @param b_exp Vector of genetic effects on exposure
#' @param b_out Vector of genetic effects on outcome
#' @param se_exp Standard errors of genetic effects on exposure
#' @param se_out Standard errors of genetic effects on outcome
#'
#' @export
#' @return List with the following elements:
#'         b: causal effect estimate
#'         se: standard error
#'         pval: p-value
#'         Q, Q_df, Q_pval: Heterogeneity stats
mr_meta_random <- function(b_exp, b_out, se_exp, se_out, parameters)
{
	if(sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 2)
	{
		return(list(b=NA, se=NA, pval=NA, nsnp=NA, Q =NA, Q_df =NA, Q_pval =NA))
	}
	ratio <- b_out / b_exp
	ratio.se <- sqrt((se_out^2/b_exp^2) + (b_out^2/b_exp^4)*se_exp^2 - 2*(b_out/b_exp^3)*parameters$Cov)
	res <- meta::metagen(ratio, ratio.se)
	b <- res$TE.random
	se <- res$seTE.random
	pval <- res$pval.random
	Q_pval <- pchisq(res$Q, res$df.Q, low=FALSE)
	return(list(b=b, se=se, pval=pval, nsnp=length(b_exp), Q = res$Q, Q_df = res$df.Q, Q_pval = Q_pval))
}



#' Maximum likelihood MR method
#'
#' @param b_exp Vector of genetic effects on exposure
#' @param b_out Vector of genetic effects on outcome
#' @param se_exp Standard errors of genetic effects on exposure
#' @param se_out Standard errors of genetic effects on outcome
#'
#' @export
#' @return List with the following elements:
#'         b: causal effect estimate
#'         se: standard error
#'         pval: p-value
#'         Q, Q_df, Q_pval: Heterogeneity stats
mr_two_sample_ml <- function(b_exp, b_out, se_exp, se_out, parameters)
{
	if(sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 2)
	{
		return(list(b=NA, se=NA, pval=NA, nsnp=NA, Q=NA, Q_df=NA, Q_pval=NA))
	}
	loglikelihood <- function(param) {
		return(1/2*sum((b_exp-param[1:length(b_exp)])^2/se_exp^2)+1/2*sum((b_out-param[length(b_exp)+1]*param[1:length(b_exp)])^2/se_out^2))
	}
	opt <- try(optim(
		c(b_exp, sum(b_exp*b_out/se_out^2)/sum(b_exp^2/se_out^2)),
		loglikelihood, 
		hessian=TRUE, 
		control = list(maxit=25000)), silent=TRUE)
	if(class(opt)=="try-error")
	{
		message("mr_two_sample_ml failed to converge")
		return(list(b=NA, se=NA, pval=NA, nsnp=NA, Q=NA, Q_df=NA, Q_pval=NA))
	}

	b <- opt$par[length(b_exp)+1]
	se <- try(sqrt(solve(opt$hessian)[length(b_exp)+1,length(b_exp)+1]))
	if(class(se)=="try-error")
	{
		message("mr_two_sample_ml failed to converge")
		return(list(b=NA, se=NA, pval=NA, nsnp=NA, Q=NA, Q_df=NA, Q_pval=NA))
	}

	pval <- 2 * pnorm(abs(b) / se, low=FALSE)

	Q <- 2 * opt$value
	Q_df <- length(b_exp) - 1
	Q_pval <- pchisq(Q, Q_df, low=FALSE)

	return(list(b=b, se=se, pval=pval, nsnp=length(b_exp), Q = Q, Q_df = Q_df, Q_pval = Q_pval))
}



#' Egger's regression for Mendelian randomization
#'
#' @param b_exp Vector of genetic effects on exposure
#' @param b_out Vector of genetic effects on outcome
#' @param se_exp Standard errors of genetic effects on exposure
#' @param se_out Standard errors of genetic effects on outcome
#' @param bootstrap Number of bootstraps to estimate standard error. If NULL then don't use bootstrap
#'
#' @export
#' @return List of with the following elements:
#'         b: MR estimate
#'         se: Standard error of MR estimate
#'         pval: p-value of MR estimate
#'         b_i: Estimate of horizontal pleiotropy (intercept)
#'         se_i: Standard error of intercept
#'         pval_i: p-value of intercept
#'         Q, Q_df, Q_pval: Heterogeneity stats
#'         mod: Summary of regression
#'         dat: Original data used for MR Egger regression
mr_egger_regression <- function(b_exp, b_out, se_exp, se_out, parameters)
{
	stopifnot(length(b_exp) == length(b_out))
	stopifnot(length(se_exp) == length(se_out))
	stopifnot(length(b_exp) == length(se_out))

	# print(b_exp)

	nulllist <- list(
			b = NA,
			se = NA,
			pval = NA,
			nsnp = NA,
			b_i = NA,
			se_i = NA,
			pval_i = NA,
			Q = NA,
			Q_df = NA,
			Q_pval = NA,
			mod = NA,
			smod = NA,
			dat = NA
		)
	if(sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 3)
	{
		return(nulllist)
	}

	sign0 <- function(x)
	{
		x[x==0] <- 1
		return(sign(x))
	}

	to_flip <- sign0(b_exp) == -1
	b_out = b_out*sign0(b_exp)
	b_exp = abs(b_exp) 
	dat <- data.frame(b_out=b_out, b_exp=b_exp, se_exp=se_exp, se_out=se_out, flipped=to_flip)
	mod <- lm(b_out ~ b_exp, weights=1/se_out^2)
	smod <- summary(mod)
	if(nrow(coefficients(smod)) > 1)
	{
		b <- coefficients(smod)[2,1]
		se <- coefficients(smod)[2,2] / min(1,smod$sigma)
		pval <- 2 * pt(abs(b / se), length(b_exp) - 2, lower.tail = FALSE)
		b_i <- coefficients(smod)[1,1]
		se_i <- coefficients(smod)[1,2] / min(1,smod$sigma)
		pval_i <- 2 * pt(abs(b_i / se_i), length(b_exp) - 2, lower.tail = FALSE)

		Q <- smod$sigma^2 * (length(b_exp) - 2)
		Q_df <- length(b_exp) - 2
		Q_pval <- pchisq(Q, Q_df, low=FALSE)
	} else {
		warning("Collinearities in MR Egger, try LD pruning the exposure variables.")
		return(nulllist)
	}
 	return(list(b = b, se = se, pval = pval, nsnp = length(b_exp), b_i = b_i, se_i = se_i, pval_i = pval_i, Q = Q, Q_df = Q_df, Q_pval = Q_pval, mod = smod, dat = dat))
}




linreg <- function(x, y, w=rep(x,1))
{
	xp <- w*x
	yp <- w*y
	t(xp) %*% yp / (t(xp)%*%xp)

	bhat <- cov(x*w,y*w, use="pair") / var(x*w, na.rm=T)
	ahat <- mean(y, na.rm=T) - mean(x, na.rm=T) * bhat
	yhat <- ahat + bhat * x
	se <- sqrt(sum((yp - yhat)^2) / (sum(!is.na(yhat)) - 2) / t(x)%*%x )

	sum(w * (y-yhat)^2)
	se <- sqrt(sum(w*(y-yhat)^2) /  (sum(!is.na(yhat)) - 2) / (sum(w*x^2)))
	pval <- 2 * pnorm(abs(bhat / se), low=FALSE)
	return(list(ahat=ahat,bhat=bhat,se=se, pval=pval))
}


#' Run bootstrap to generate standard errors for MR
#'
#' @param b_exp Vector of genetic effects on exposure
#' @param b_out Vector of genetic effects on outcome
#' @param se_exp Standard errors of genetic effects on exposure
#' @param se_out Standard errors of genetic effects on outcome
#' @param nboot Number of bootstraps. Default 1000
#'
#' @export
#' @return List of with the following elements:
#'         b: MR estimate
#'         se: Standard error of MR estimate
#'         pval: p-value of MR estimate
#'         b_i: Estimate of horizontal pleiotropy (intercept)
#'         se_i: Standard error of intercept
#'         pval_i: p-value of intercept
#'         mod: Summary of regression
#'         dat: Original data used for MR Egger regression
mr_egger_regression_bootstrap <- function(b_exp, b_out, se_exp, se_out, parameters)
{
	if(sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 3)
	{
		return(list(
			b = NA,
			se = NA,
			pval = NA,
			nsnp = NA,
			b_i = NA,
			se_i = NA,
			pval_i = NA,
			mod = NA,
			smod = NA,
			dat = NA
		))
	}
	nboot <- parameters$nboot
	# Do bootstraps
	res <- array(0, c(nboot+1, 2))
	# pb <- txtProgressBar(min = 0, max = nboot, initial = 0, style=3) 
	for (i in 1:nboot)
	{
		# setTxtProgressBar(pb, i)
		#sample from distributions of SNP betas
		xs <- rnorm(length(b_exp),b_exp,se_exp)
		ys <- rnorm(length(b_out),b_out,se_out)

		# Use absolute values for Egger reg
		ys <- ys*sign(xs)
		xs <- abs(xs)

		#weighted regression with given formula
		# r <- summary(lm(ys ~ xs, weights=1/se_out^2))
		r <- linreg(xs, ys, 1/se_out^2)

		#collect coefficient from given line.
		res[i, 1] <- r$ahat
		res[i, 2] <- r$bhat
		# res[i, 1] <- r$coefficients[1,1]
		# res[i, 2] <- r$coefficients[2,1]

	}
	cat("\n")

	return(list(b = mean(res[,2], na.rm=T), se = sd(res[,2], na.rm=T), pval = sum(sign(mean(res[,2],na.rm=T)) * res[,2] < 0)/nboot, nsnp = length(b_exp), b_i = mean(res[,1], na.rm=T), se_i = sd(res[,1], na.rm=T), pval_i = sum(sign(mean(res[,1],na.rm=T)) * res[,1] < 0)/nboot))
}


#' Weighted median method
#'
#' Perform MR using summary statistics. Bootstraps used to calculate standard error.
#'
#' @param b_exp Vector of genetic effects on exposure
#' @param b_out Vector of genetic effects on outcome
#' @param se_exp Standard errors of genetic effects on exposure
#' @param se_out Standard errors of genetic effects on outcome
#' @param nboot Number of bootstraps to calculate se. Default 1000
#'
#' @export
#' @return List with the following elements:
#'         b: MR estimate
#'         se: Standard error
#'         pval: p-value
mr_weighted_median <- function(b_exp, b_out, se_exp, se_out, parameters=default_parameters())
{
	if(sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 3)
	return(list(b=NA, se=NA, pval=NA, nsnp=NA))

	b_iv <- b_out / b_exp
	VBj <- ((se_out)^2)/(b_exp)^2 + (b_out^2)*((se_exp^2))/(b_exp)^4
	b <- weighted_median(b_iv, 1 / VBj)
	se <- weighted_median_bootstrap(b_exp, b_out, se_exp, se_out, 1 / VBj, parameters$nboot)
	pval <- 2 * pnorm(abs(b/se), low=FALSE)
	return(list(b=b, se=se, pval=pval, Q=NA, Q_df=NA, Q_pval=NA, nsnp=length(b_exp)))
}


#' Simple median method
#'
#' Perform MR using summary statistics. Bootstraps used to calculate standard error.
#'
#' @param b_exp Vector of genetic effects on exposure
#' @param b_out Vector of genetic effects on outcome
#' @param se_exp Standard errors of genetic effects on exposure
#' @param se_out Standard errors of genetic effects on outcome
#' @param nboot Number of bootstraps to calculate se. Default 1000
#'
#' @export
#' @return List with the following elements:
#'         b: MR estimate
#'         se: Standard error
#'         pval: p-value
mr_simple_median <- function(b_exp, b_out, se_exp, se_out, parameters=default_parameters())
{
	if(sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 3)
	return(list(b=NA, se=NA, pval=NA, nsnp=NA))

	b_iv <- b_out / b_exp
	b <- weighted_median(b_iv, rep(1/length(b_exp), length(b_exp)))
	se <- weighted_median_bootstrap(b_exp, b_out, se_exp, se_out, rep(1/length(b_exp), length(b_exp)), parameters$nboot)
	pval <- 2 * pnorm(abs(b/se), low=FALSE)
	return(list(b=b, se=se, pval=pval, nsnp=length(b_exp)))
}



#' Weighted median method
#'
#' New method from Jack
#'
#' @param b_iv Wald ratios
#' @param weights Weights of each SNP
#'
#' @export
#' @return MR estimate
weighted_median <- function(b_iv, weights)
{
	betaIV.order <- b_iv[order(b_iv)]
	weights.order <- weights[order(b_iv)]
	weights.sum <- cumsum(weights.order)-0.5*weights.order
	weights.sum <- weights.sum/sum(weights.order)
	below <- max(which(weights.sum<0.5))
	b = betaIV.order[below] + (betaIV.order[below+1]-betaIV.order[below])*
	(0.5-weights.sum[below])/(weights.sum[below+1]-weights.sum[below])
	return(b)
}

weighted_median <- function(b_iv, weights)
{
	betaIV.order <- b_iv[order(b_iv)]
	weights.order <- weights[order(b_iv)]
	weights.sum <- cumsum(weights.order)-0.5*weights.order
	weights.sum <- weights.sum/sum(weights.order)
	below <- max(which(weights.sum<0.5))
	b = betaIV.order[below] + (betaIV.order[below+1]-betaIV.order[below])*
	(0.5-weights.sum[below])/(weights.sum[below+1]-weights.sum[below])
	return(b)
}





#' Calculate standard errors for weighted median method using bootstrap
#'
#' Based on new script for weighted median confidence interval, update 31 July 2015
#'
#' @param b_exp Vector of genetic effects on exposure
#' @param b_out Vector of genetic effects on outcome
#' @param se_exp Standard errors of genetic effects on exposure
#' @param se_out Standard errors of genetic effects on outcome
#' @param weights Weights to apply to each SNP
#' @param nboot Number of bootstraps. Default 1000
#'
#' @export
#' @return Empirical standard error
weighted_median_bootstrap <- function(b_exp, b_out, se_exp, se_out, weights, nboot)
{
	med <- rep(0, nboot)
	for(i in 1:nboot){
		b_exp.boot = rnorm(length(b_exp), mean=b_exp, sd=se_exp)
		b_out.boot = rnorm(length(b_out), mean=b_out, sd=se_out)
		betaIV.boot = b_out.boot/b_exp.boot
		med[i] = weighted_median(betaIV.boot, weights)
	}
	return(sd(med)) 
}



#' Penalised weighted median MR
#'
#' Modification to standard weighted median MR
#' Updated based on Burgess 2016 "Robust instrumental variable methods using multiple candidate instruments with application to Mendelian randomization"
#'
#' @param b_exp Vector of genetic effects on exposure
#' @param b_out Vector of genetic effects on outcome
#' @param se_exp Standard errors of genetic effects on exposure
#' @param se_out Standard errors of genetic effects on outcome
#' @param penk Constant term in penalisation. Default=20 
#' @param nboot Number of bootstraps to calculate SE. Default 1000
#'
#' @export
#' @return List with the following elements:
#'         b: MR estimate
#'         se: Standard error
#'         pval: p-value
mr_penalised_weighted_median <- function(b_exp, b_out, se_exp, se_out, parameters=default_parameters())
{
	if(sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 3)
	return(list(b=NA, se=NA, pval=NA, nsnp=NA))
	betaIV <- b_out/b_exp # ratio estimates
	betaIVW <- sum(b_out*b_exp*se_out^-2)/sum(b_exp^2*se_out^-2) # IVW estimate
	VBj <- ((se_out)^2)/(b_exp)^2 + (b_out^2)*((se_exp^2))/(b_exp)^4
	weights <- 1/VBj
	bwm <- mr_weighted_median(b_exp, b_out, se_exp, se_out, parameters)
	penalty <- pchisq(weights*(betaIV-bwm$b)^2, df=1, lower.tail=FALSE)
	pen.weights <- weights*pmin(1, penalty*parameters$penk) # penalized weights
	b <- weighted_median(betaIV, pen.weights) # penalized weighted median estimate
	se <- weighted_median_bootstrap(b_exp, b_out, se_exp, se_out, pen.weights, parameters$nboot)
	pval <- 2 * pnorm(abs(b/se), low=FALSE)
	return(list(b = b, se = se, pval=pval, nsnp=length(b_exp)))
}


#' MR weighted mode estimators
#'
#' Weighted mode estimator
#'
#' @param dat Output from harmonise_data()
#' @param parameters=default_parameters() <what param does>
#'
#' @export
#' @return data frame
mr_weighted_mode <- function(b_exp, b_out, se_exp, se_out, parameters=default_parameters()) 
{
	index <- !is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)
	if(sum(index) < 3)
	return(list(b=NA, se=NA, pval=NA, nsnp=NA))

	b_exp <- b_exp[index]
	b_out <- b_out[index]
	se_exp <- se_exp[index]
	se_out <- se_out[index]

	#--------------------------------------#
	#Function to compute the point estimate#
	#--------------------------------------#
	#BetaIV.in: ratio estimates
	#seBetaIV.in: standard errors of ratio estimates
	beta <- function(BetaIV.in, seBetaIV.in, phi)
	{
		#Bandwidth rule - modified Silverman's rule proposed by Bickel (2002)
		s <- 0.9*(min(sd(BetaIV.in), mad(BetaIV.in)))/length(BetaIV.in)^(1/5)

		#Standardised weights
		weights <- seBetaIV.in^-2/sum(seBetaIV.in^-2)

		beta <- NULL

		for(cur_phi in phi)
		{
			#Define the actual bandwidth
			h <- s*cur_phi
			#Compute the smoothed empirical density function
			densityIV <- density(BetaIV.in, weights=weights, bw=h)
			#Extract the point with the highest density as the point estimate 
			beta[length(beta)+1] <- densityIV$x[densityIV$y==max(densityIV$y)]
		}
		return(beta)
	}

	#------------------------------------------#
	#Function to estimate SEs through bootstrap#
	#------------------------------------------#
	#BetaIV.in: ratio estimates
	#seBetaIV.in: standard errors of ratio estimates
	#beta_Mode.in: point causal effect estimates
	boot <- function(BetaIV.in, seBetaIV.in, beta_Mode.in, nboot)
	{
		#Set up a matrix to store the results from each bootstrap iteration
		beta.boot <- matrix(nrow=nboot, ncol=length(beta_Mode.in))

		for(i in 1:nboot) 
		{
			#Re-sample each ratio estimate using SEs derived not assuming NOME
			BetaIV.boot      <- rnorm(length(BetaIV.in), mean=BetaIV.in, sd=seBetaIV.in[,1])
			#Re-sample each ratio estimate using SEs derived under NOME
			BetaIV.boot_NOME <- rnorm(length(BetaIV.in), mean=BetaIV.in, sd=seBetaIV.in[,2])

			#Simple mode, not assuming NOME
			beta.boot[i,1:length(phi)] <- beta(BetaIV.in=BetaIV.boot, seBetaIV.in=rep(1, length(BetaIV)), phi=phi)
			#Weighted mode, not assuming NOME
			beta.boot[i,(length(phi)+1):(2*length(phi))] <- beta(BetaIV.in=BetaIV.boot, seBetaIV.in=seBetaIV.in[,1], phi=phi)
			#Simple mode, assuming NOME
			beta.boot[i,(2*length(phi)+1):(3*length(phi))] <- beta(BetaIV.in=BetaIV.boot_NOME, seBetaIV.in=rep(1, length(BetaIV)), phi=phi)
			#Weighted mode, assuming NOME
			beta.boot[i,(3*length(phi)+1):(4*length(phi))] <- beta(BetaIV.in=BetaIV.boot_NOME, seBetaIV.in=seBetaIV.in[,2], phi=phi)
		}
		return(beta.boot)
	}

	# Parameters
	phi <- parameters$phi
	nboot <- parameters$nboot
	alpha <- parameters$alpha

	#Ratio estimates
	BetaIV   <- b_out/b_exp

	#SEs of ratio estimates
	seBetaIV <- cbind(sqrt((se_out^2)/(b_exp^2) + ((b_out^2)*(se_exp^2))/(b_exp^4)), #SEs NOT assuming NOME
	se_out/abs(b_exp)) #SEs ASSUMING NOME

	#Point causal effect estimate using the simple mode
	beta_SimpleMode <- beta(BetaIV.in=BetaIV, seBetaIV.in=rep(1, length(BetaIV)), phi=phi)

	#Point causal effect estimate using the weighted mode (not asusming NOME)
	beta_WeightedMode <- beta(BetaIV.in=BetaIV, seBetaIV.in=seBetaIV[,1], phi=phi)

	#Point causal effect estimate using the weighted mode (asusming NOME)
	beta_WeightedMode_NOME <- beta(BetaIV.in=BetaIV, seBetaIV.in=seBetaIV[,2], phi=phi)

	#Combine all point effect estimates in a single vector
	beta_Mode <- rep(c(beta_SimpleMode, beta_WeightedMode,
	beta_SimpleMode, beta_WeightedMode_NOME))

	#Compute SEs, confidence intervals and P-value
	beta_Mode.boot <- boot(BetaIV.in=BetaIV, seBetaIV.in=seBetaIV, beta_Mode.in=beta_Mode, nboot=nboot)
	se_Mode <- apply(beta_Mode.boot, 2, mad)

	CIlow_Mode <- beta_Mode-qnorm(1-alpha/2)*se_Mode
	CIupp_Mode <- beta_Mode+qnorm(1-alpha/2)*se_Mode

	P_Mode <- pt(abs(beta_Mode/se_Mode), df=length(b_exp)-1, lower.tail=F)*2

	#Vector to indicate the method referring to each row
	Method <- rep(c('Simple mode', 'Weighted mode', 'Simple mode (NOME)', 'Weighted mode (NOME)'), each=length(phi))

	#Return a data frame containing the results
	Results <- data.frame(Method, length(b_exp), beta_Mode, se_Mode, CIlow_Mode, CIupp_Mode, P_Mode)
	colnames(Results) <- c('Method', 'nsnp', 'Estimate', 'SE', 'CI_low', 'CI_upp', 'P')

	i <- which(Method == "Weighted mode")

	return(list(b = Results$Estimate[i], se = Results$SE[i], pval=Results$P[i], nsnp=length(b_exp)))
}



#' MR median estimators
#'
#' @param dat Output from harmonise_data()
#' @param parameters=default_parameters() <what param does>
#'
#' @export
#' @return data frame
mr_median <- function(dat, parameters=default_parameters())
{
	if("mr_keep" %in% names(dat)) dat <- subset(dat, mr_keep)

	if(nrow(dat) < 3) 
	{
		warning("Need at least 3 SNPs")
		return(NULL)
	}

	b_exp <- dat$beta.exposure
	b_out <- dat$beta.outcome
	se_exp <- dat$se.exposure
	se_out <- dat$se.outcome

	sm <- mr_simple_median(b_exp, b_out, se_exp, se_out, parameters)
	wm <- mr_weighted_median(b_exp, b_out, se_exp, se_out, parameters)
	pm <- mr_penalised_weighted_median(b_exp, b_out, se_exp, se_out, parameters)

	res <- data.frame(
		Method = c("Simple median", "Weighted median", "Penalised median"),
		nsnp = length(b_exp),
		Estimate = c(sm$b, wm$b, pm$b),
		SE = c(sm$se, wm$se, pm$se)
	)
	res$CI_low <- res$Estimate - qnorm(1-parameters$alpha/2) * res$SE
	res$CI_upp <- res$Estimate + qnorm(1-parameters$alpha/2) * res$SE
	res$P <- c(sm$pval, wm$pval, pm$pval)
	return(res)
}


#' Inverse variance weighted regression
#'
#' @param b_exp Vector of genetic effects on exposure
#' @param b_out Vector of genetic effects on outcome
#' @param se_exp Standard errors of genetic effects on exposure
#' @param se_out Standard errors of genetic effects on outcome
#'
#' @export
#' @return List with the following elements:
#'         b: MR estimate
#'         se: Standard error
#'         pval: p-value
#'         Q, Q_df, Q_pval: Heterogeneity stats
mr_ivw <- function(b_exp, b_out, se_exp, se_out, parameters=default_parameters())
{
	if(sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 2)
	return(list(b=NA, se=NA, pval=NA, nsnp=NA))

	ivw.res <- summary(lm(b_out ~ -1 + b_exp, weights = 1/se_out^2))
	b <- ivw.res$coef["b_exp","Estimate"]
	se <- ivw.res$coef["b_exp","Std. Error"]/min(1,ivw.res$sigma) #sigma is the residual standard error 
	pval <- 2 * pnorm(abs(b/se), low=FALSE)
	Q <- ivw.res$sigma^2*(length(b_exp)-2)
	Q_df <- length(b_exp) - 1
	Q_pval <- pchisq(Q, Q_df, low=FALSE)
	# from formula phi =  Q/DF rearranged to to Q = phi*DF, where phi is sigma^2
	# Q.ivw<-sum((1/(se_out/b_exp)^2)*(b_out/b_exp-ivw.reg.beta)^2)
	return(list(b = b, se = se, pval = pval, nsnp=length(b_exp), Q = Q, Q_df = Q_df, Q_pval = Q_pval))
}



#' Leave one out sensitivity analysis
#'
#' @param dat Output from \code{harmonise_exposure_outcome}
#' @param method=mr_ivw Choose which method to use
#'
#' @export
#' @return List of data frames
mr_leaveoneout <- function(dat, parameters=default_parameters(), method=mr_ivw)
{
	res <- plyr::ddply(dat, c("id.exposure", "id.outcome"), function(X)
	{
		x <- subset(X, mr_keep)
		nsnp <- nrow(x)
		if(nsnp == 0)
		{
			x <- X[1,]
			d <- data.frame(
				SNP = "All",
				b = NA,
				se = NA,
				p = NA,
				samplesize = NA,
				outcome = x$outcome[1],
				exposure = x$exposure[1]
			)
			return(d)
		}
		if(nsnp > 2)
		{
			l <- lapply(1:nsnp, function(i)
			{
				with(x, method(beta.exposure[-i], beta.outcome[-i], se.exposure[-i], se.outcome[-i], parameters))
			})
			l[[nsnp+1]] <- with(x, method(beta.exposure, beta.outcome, se.exposure, se.outcome, parameters))
			d <- data.frame(
				SNP = c(as.character(x$SNP), "All"),
				b = sapply(l, function(y) y$b),
				se = sapply(l, function(y) y$se),
				p = sapply(l, function(y) y$pval),
				samplesize = x$samplesize.outcome[1]
			)
			d$outcome <- x$outcome[1]
			d$exposure <- x$exposure[1]

		} else {
			a <- with(x, method(beta.exposure, beta.outcome, se.exposure, se.outcome, parameters))
			d <- data.frame(
				SNP = "All",
				b = a$b,
				se = a$se,
				p = a$pval,
				samplesize = x$samplesize.outcome[1]
			)
			d$outcome <- x$outcome[1]
			d$exposure <- x$exposure[1]
		}
		return(d)
	})
	res <- subset(res, select=c(exposure, outcome, id.exposure, id.outcome, samplesize, SNP, b, se, p))
	return(res)
}


#' Perform 2 sample MR on each SNP individually
#'
#' @param dat Output from \code{harmonise_exposure_outcome}
#' @param method=mr_two_sample_ml Function to use for MR analysis
#'
#' @export
#' @return List of data frames
mr_singlesnp <- function(dat, parameters=default_parameters(), single_method="mr_wald_ratio", all_method=c("mr_ivw", "mr_egger_regression"))
{
	res <- plyr::ddply(dat, c("id.exposure", "id.outcome"), function(X)
	{
		x <- subset(X, mr_keep)
		nsnp <- nrow(x)
		if(nsnp == 0)
		{
			x <- X[1,]
			d <- data.frame(
				SNP = "No available data",
				b = NA,
				se = NA,
				p = NA,
				samplesize = NA,
				outcome = x$outcome[1],
				exposure = x$exposure[1]
			)
			return(d)
		}
		l <- lapply(1:nsnp, function(i)
		{
			with(x, get(single_method)(beta.exposure[i], beta.outcome[i], se.exposure[i], se.outcome[i], parameters))
		})
		nom <- c()
		for(i in 1:length(all_method))
		{
			l[[nsnp+i]] <- with(x, get(all_method[i])(beta.exposure, beta.outcome, se.exposure, se.outcome, parameters))

			nom <- c(nom, paste0("All - ", subset(mr_method_list(), obj==all_method[i])$name))
		}

		d <- data.frame(
			SNP = c(as.character(x$SNP), nom),
			b = sapply(l, function(y) y$b),
			se = sapply(l, function(y) y$se),
			p = sapply(l, function(y) y$pval),
			samplesize = x$samplesize.outcome[1]
		)
		d$outcome <- x$outcome[1]
		d$exposure <- x$exposure[1]
		return(d)
	})
	res <- subset(res, select=c(exposure, outcome, id.exposure, id.outcome, samplesize, SNP, b, se, p))
	return(res)
}


#' Get heterogeneity stats
#'
#' @param dat Harmonised exposure and outcome data. Output from \code{harmonise_exposure_outcome}
#' @param parameters Parameters to be used for various MR methods. Default is output from \code{dafault_param}.
#' @param method_list List of methods to use in analysis. See \code{mr_method_list()} for details.
#'
#' @export
#' @return Data frame
mr_heterogeneity <- function(dat, parameters=default_parameters(), method_list = subset(mr_method_list(), heterogeneity_test)$obj)
{
	het_tab <- plyr::ddply(dat, c("id.exposure", "id.outcome"), function(x1)
	{
		# message("Performing MR analysis of '", x$id.exposure[1], "' on '", x$id.outcome[1], "'")
		x <- subset(x1, mr_keep)
		if(nrow(x) < 2)
		{
			message("Not enough SNPs available for Heterogeneity analysis of '", x1$id.exposure[1], "' on '", x1$id.outcome[1], "'")
			return(NULL)
		}
		res <- lapply(method_list, function(meth)
		{
			get(meth)(x$beta.exposure, x$beta.outcome, x$se.exposure, x$se.outcome, parameters)	
		})
		methl <- mr_method_list()
		het_tab <- data.frame(
			outcome = x$outcome[1],
			exposure = x$exposure[1],
			method = methl$name[methl$obj %in% method_list],
			Q = sapply(res, function(x) x$Q),
			Q_df = sapply(res, function(x) x$Q_df),
			Q_pval = sapply(res, function(x) x$Q_pval)
		)
		het_tab <- subset(het_tab, !(is.na(Q) & is.na(Q_df) & is.na(Q_pval)))
		return(het_tab)
	})

	return(het_tab)
}


Isq <- function(y,s)
{
	k = length(y)
	w = 1/s^2; 
	sum.w = sum(w)
	mu.hat = sum(y*w)/sum.w
	Q = sum(w*(y-mu.hat)^2)
	Isq = (Q - (k-1))/Q
	Isq = max(0,Isq)
	return(Isq)
}



#' Test for horizontal pleiotropy in MR analysis
#'
#' Performs MR Egger and returns intercept values
#'
#' @param dat Harmonised exposure and outcome data. Output from \code{harmonise_exposure_outcome}
#'
#' @export
#' @return data frame
mr_pleiotropy_test <- function(dat)
{
	ptab <- plyr::ddply(dat, c("id.exposure", "id.outcome"), function(x1)
	{
		x <- subset(x1, mr_keep)
		if(nrow(x) < 2)
		{
			message("Not enough SNPs available for pleiotropy analysis of '", x1$id.exposure[1], "' on '", x1$id.outcome[1], "'")
			return(NULL)
		}
		res <- mr_egger_regression(x$beta.exposure, x$beta.outcome, x$se.exposure, x$se.outcome, default_parameters())
		out <- data.frame(
			outcome = x$outcome[1],
			exposure = x$exposure[1],
			egger_intercept = res$b_i,
			se = res$se_i,
			pval = res$pval_i
		)
		return(out)
	})
	return(ptab)
}


#' MR mode estimators
#'
#' <full description>
#'
#' @param dat Output from harmonise_data()
#' @param parameters=default_parameters() <what param does>
#'
#' @export
#' @return data frame
mr_mode <- function(dat, parameters=default_parameters()) 
{
	if("mr_keep" %in% names(dat)) dat <- subset(dat, mr_keep)

	if(nrow(dat) < 3) 
	{
		warning("Need at least 3 SNPs")
		return(NULL)
	}

	b_exp <- dat$beta.exposure
	b_out <- dat$beta.outcome
	se_exp <- dat$se.exposure
	se_out <- dat$se.outcome

	#--------------------------------------#
	#Function to compute the point estimate#
	#--------------------------------------#
	#BetaIV.in: ratio estimates
	#seBetaIV.in: standard errors of ratio estimates
	beta <- function(BetaIV.in, seBetaIV.in, phi)
	{
		#Bandwidth rule - modified Silverman's rule proposed by Bickel (2002)
		s <- 0.9*(min(sd(BetaIV.in), mad(BetaIV.in)))/length(BetaIV.in)^(1/5)

		#Standardised weights
		weights <- seBetaIV.in^-2/sum(seBetaIV.in^-2)

		beta <- NULL

		for(cur_phi in phi)
		{
			#Define the actual bandwidth
			h <- s*cur_phi
			#Compute the smoothed empirical density function
			densityIV <- density(BetaIV.in, weights=weights, bw=h)
			#Extract the point with the highest density as the point estimate 
			beta[length(beta)+1] <- densityIV$x[densityIV$y==max(densityIV$y)]
		}
		return(beta)
	}

	#------------------------------------------#
	#Function to estimate SEs through bootstrap#
	#------------------------------------------#
	#BetaIV.in: ratio estimates
	#seBetaIV.in: standard errors of ratio estimates
	#beta_Mode.in: point causal effect estimates
	boot <- function(BetaIV.in, seBetaIV.in, beta_Mode.in, nboot)
	{
		#Set up a matrix to store the results from each bootstrap iteration
		beta.boot <- matrix(nrow=nboot, ncol=length(beta_Mode.in))

		for(i in 1:nboot) 
		{
			#Re-sample each ratio estimate using SEs derived not assuming NOME
			BetaIV.boot      <- rnorm(length(BetaIV.in), mean=BetaIV.in, sd=seBetaIV.in[,1])
			#Re-sample each ratio estimate using SEs derived under NOME
			BetaIV.boot_NOME <- rnorm(length(BetaIV.in), mean=BetaIV.in, sd=seBetaIV.in[,2])

			#Simple mode, not assuming NOME
			beta.boot[i,1:length(phi)] <- beta(BetaIV.in=BetaIV.boot, seBetaIV.in=rep(1, length(BetaIV)), phi=phi)
			#Weighted mode, not assuming NOME
			beta.boot[i,(length(phi)+1):(2*length(phi))] <- beta(BetaIV.in=BetaIV.boot, seBetaIV.in=seBetaIV.in[,1], phi=phi)
			#Simple mode, assuming NOME
			beta.boot[i,(2*length(phi)+1):(3*length(phi))] <- beta(BetaIV.in=BetaIV.boot_NOME, seBetaIV.in=rep(1, length(BetaIV)), phi=phi)
			#Weighted mode, assuming NOME
			beta.boot[i,(3*length(phi)+1):(4*length(phi))] <- beta(BetaIV.in=BetaIV.boot_NOME, seBetaIV.in=seBetaIV.in[,2], phi=phi)
		}
		return(beta.boot)
	}

	# Parameters
	phi <- parameters$phi
	nboot <- parameters$nboot
	alpha <- parameters$alpha

	#Ratio estimates
	BetaIV   <- b_out/b_exp

	#SEs of ratio estimates
	seBetaIV <- cbind(sqrt((se_out^2)/(b_exp^2) + ((b_out^2)*(se_exp^2))/(b_exp^4)), #SEs NOT assuming NOME
	se_out/abs(b_exp)) #SEs ASSUMING NOME

	#Point causal effect estimate using the simple mode
	beta_SimpleMode <- beta(BetaIV.in=BetaIV, seBetaIV.in=rep(1, length(BetaIV)), phi=phi)

	#Point causal effect estimate using the weighted mode (not asusming NOME)
	beta_WeightedMode <- beta(BetaIV.in=BetaIV, seBetaIV.in=seBetaIV[,1], phi=phi)

	#Point causal effect estimate using the weighted mode (asusming NOME)
	beta_WeightedMode_NOME <- beta(BetaIV.in=BetaIV, seBetaIV.in=seBetaIV[,2], phi=phi)

	#Combine all point effect estimates in a single vector
	beta_Mode <- rep(c(beta_SimpleMode, beta_WeightedMode,
	beta_SimpleMode, beta_WeightedMode_NOME))

	#Compute SEs, confidence intervals and P-value
	beta_Mode.boot <- boot(BetaIV.in=BetaIV, seBetaIV.in=seBetaIV, beta_Mode.in=beta_Mode, nboot=nboot)
	se_Mode <- apply(beta_Mode.boot, 2, mad)

	CIlow_Mode <- beta_Mode-qnorm(1-alpha/2)*se_Mode
	CIupp_Mode <- beta_Mode+qnorm(1-alpha/2)*se_Mode

	P_Mode <- pt(abs(beta_Mode/se_Mode), df=length(b_exp)-1, lower.tail=F)*2

	#Vector to indicate the method referring to each row
	Method <- rep(c('Simple mode', 'Weighted mode', 'Simple mode (NOME)', 'Weighted mode (NOME)'), each=length(phi))

	#Return a data frame containing the results
	Results <- data.frame(Method, length(b_exp), beta_Mode, se_Mode, CIlow_Mode, CIupp_Mode, P_Mode)
	colnames(Results) <- c('Method', 'nsnp', 'Estimate', 'SE', 'CI_low', 'CI_upp', 'P')

	return(Results)
}






#' Perform Rucker, Median and Mode
#'
#' @param dat <what param does>
#' @param parameters=default_parameters() <what param does>
#'
#' @export
#' @return list
mr_all <- function(dat, parameters=default_parameters())
{
	dat <- subset(dat, mr_keep)
	res <- dplyr::group_by(dat, exposure, outcome) %>%
		do({
			x <- .
			if(nrow(x) == 1)
			{
				a <- mr_wald_ratio(x$beta.exposure, x$beta.outcome, x$se.exposure, x$se.outcome)
				out <- tibble(
					exposure = x$exposure[1],
					outcome = x$outcome[1],
					id.exposure = x$id.exposure[1],
					id.outcome = x$id.outcome[1],
					Method = "Wald ratio",
					Estimate = a$b,
					SE = a$se,
					CI_low = a$b - 1.96 * a$se,
					CI_upp = a$b + 1.96 * a$se,
					P = a$pval,
					nsnp = 1
				)
				return(out)
			} else if(nrow(x) <= 3) {
				a <- mr_ivw(x$beta.exposure, x$beta.outcome, x$se.exposure, x$se.outcome)
				out <- tibble(
					exposure = x$exposure[1],
					outcome = x$outcome[1],
					id.exposure = x$id.exposure[1],
					id.outcome = x$id.outcome[1],
					Method = "IVW fixed effects",
					Estimate = a$b,
					SE = a$se,
					CI_low = a$b - 1.96 * a$se,
					CI_upp = a$b + 1.96 * a$se,
					P = a$pval,
					nsnp = nrow(x)
				)
				return(out)
			} else {
				mrrucker <- mr_rucker(x, parameters)
				mrruckercd <- mr_rucker_cooksdistance(x, parameters)
				mrmode <- mr_mode(x, parameters)
				mrmedian <- mr_median(x, parameters)
				out <- suppressWarnings(dplyr::bind_rows(
					mrrucker$rucker,
					mrrucker$selected, 
					mrruckercd$rucker,
					mrruckercd$selected,
					mrmode, 
					mrmedian
				))
				out$exposure <- x$exposure[1]
				out$outcome <- x$outcome[1]
				out$id.exposure <- x$id.exposure[1]
				out$id.outcome <- x$id.outcome[1]
				return(out)
			}
		})

	return(res)
}
