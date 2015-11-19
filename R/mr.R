
#' Perform all Mendelian randomization tests
#'
#' @param dat Harmonised exposure and outcome data. Output from \code{harmonise_exposure_outcome}
#' @param param Parameters to be used for various MR methods. Default is output from \code{dafault_param}.
#' @param method_list List of methods to use in analysis. See \code{mr_method_list()} for details.
#'
#' @export
#' @return List with the following elements:
#'         mr: Table of MR results
#'         extra: Table of extra results
mr <- function(dat, parameters=default_parameters(), method_list=mr_method_list())
{
	require(plyr)
	res <- dlply(subset(dat, mr_keep), .(id.outcome, exposure), function(x)
	{
		res <- lapply(method_list, function(meth)
		{
			print(meth)
			get(meth)(x$beta.exposure, x$beta.outcome, x$se.exposure, x$se.outcome, parameters)	
		})

		mr_tab <- data.frame(
			Study.ID = x$id.outcome[1],
			Exposure = x$exposure[1],
			Test = sapply(res, function(x) x$testname),
			n.SNPs = sapply(res, function(x) x$nsnp),
			b = sapply(res, function(x) x$b),
			se = sapply(res, function(x) x$se),
			pval = sapply(res, function(x) x$pval)
		)
		mr_tab <- subset(mr_tab, !(is.na(b) & is.na(se) & is.na(pval)))
		return(mr_tab)
	})
	mr_tab <- rbind.fill(lapply(res, function(x) x))

	ao <- available_outcomes()
	ao <- subset(ao, select=c(id, trait, trait_strict, consortium, ethnic, gender, ncase, ncontrol, sample_size, pmid, unit, sd, year))

	mr_tab$ord <- 1:nrow(mr_tab)
	mr_tab <- merge(mr_tab, ao, by.x="Study.ID", by.y="id")
	mr_tab <- mr_tab[order(mr_tab$ord), ]
	mr_tab <- subset(mr_tab, select=-c(ord))

	return(mr_tab)
}


#' Get list of available MR methods
#'
#' @export
#' @return character vector of method names
mr_method_list <- function()
{
	c(
		"mr_wald_ratio",
		"mr_meta_fixed_simple",
		"mr_meta_fixed",
		"mr_meta_random",
		"mr_two_sample_ml",
		"mr_eggers_regression",
		"mr_eggers_regression_bootstrap",
		"mr_weighted_median",
		"mr_penalised_weighted_median",
		"mr_ivw"
	)
}


#' List of parameters for use with MR functions
#'
#' @export
default_parameters <- function()
{
	list(
		nboot = 1000,
		Cov = 0,
		penk = 20
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
#'         testname: Name of test
mr_wald_ratio <- function(b_exp, b_out, se_exp, se_out, parameters)
{
	if(length(b_exp) > 1)
	{
		return(list(b=NA, se=NA, pval=NA, nsnp=NA, testname="Wald ratio"))
	}
	b <- b_out / b_exp
	se <- se_out / b_exp
	# sqrt((segd^2/gp^2) + (gd^2/gp^4)*segp^2 - 2*(gd/gp^3)) #full delta method with cov set to 0
	pval <- pnorm(abs(b)/se,lower.tail=F)*2
	return(list(b=b, se=se, pval=pval, nsnp=1, testname="Wald ratio"))
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
#'         testname: Name of test
mr_meta_fixed_simple <- function(b_exp, b_out, se_exp, se_out, parameters)
{
	b <- sum(b_exp*b_out / se_out^2) / sum(b_exp^2/se_out^2)
	se <- sqrt(1 / sum(b_exp^2/se_out^2))
	pval <- pt(abs(b) / se, df = length(b_exp)-1, low=FALSE)
	return(list(b=b, se=se, pval=pval, nsnp=length(b_exp), testname="Fixed effects meta analysis (simple SE)"))
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
#'         testname: Name of test
mr_meta_fixed <- function(b_exp, b_out, se_exp, se_out, parameters)
{
	require(meta)
	ratio <- b_out / b_exp
	ratio.se <- sqrt((se_out^2/b_exp^2) + (b_out^2/b_exp^4)*se_exp^2 - 2*(b_out/b_exp^3)*parameters$Cov)
	res <- metagen(ratio, ratio.se)
	b <- res$TE.fixed
	se <- res$seTE.fixed
	pval <- res$pval.fixed
	return(list(b=b, se=se, pval=pval, nsnp=length(b_exp), testname="Fixed effects meta analysis (delta method)"))
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
#'         testname: Name of test
mr_meta_random <- function(b_exp, b_out, se_exp, se_out, parameters)
{
	require(meta)
	ratio <- b_out / b_exp
	ratio.se <- sqrt((se_out^2/b_exp^2) + (b_out^2/b_exp^4)*se_exp^2 - 2*(b_out/b_exp^3)*parameters$Cov)
	res <- metagen(ratio, ratio.se)
	b <- res$TE.random
	se <- res$seTE.random
	pval <- res$pval.random
	return(list(b=b, se=se, pval=pval, nsnp=length(b_exp), testname="Random effects meta analysis (delta method)"))
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
#'         testname: Name of test
mr_two_sample_ml <- function(b_exp, b_out, se_exp, se_out, parameters)
{
	if(length(b_exp)<2)
	{
		return(list(b=NA, se=NA, pval=NA, nsnp=NA, testname="Maximum likelihood"))
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
		return(list(b=NA, se=NA, pval=NA, nsnp=NA, testname="Maximum likelihood"))
	}

	b <- opt$par[length(b_exp)+1]
	se <- sqrt(solve(opt$hessian)[length(b_exp)+1,length(b_exp)+1])
	pval <- pt(abs(b) / se, df = length(b_exp)-1, low=FALSE)
	return(list(b=b, se=se, pval=pval, nsnp=length(b_exp), testname="Maximum likelihood"))
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
#'         mod: Summary of regression
#'         dat: Original data used for MR Egger regression
mr_eggers_regression <- function(b_exp, b_out, se_exp, se_out, parameters)
{
	stopifnot(length(b_exp) == length(b_out))
	stopifnot(length(se_exp) == length(se_out))
	stopifnot(length(b_exp) == length(se_out))

	# print(b_exp)

	if(length(b_exp) < 2)
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
			dat = NA,
			testname = "Egger regression"
		))
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

	b <- coefficients(smod)[2,1]
	se <- coefficients(smod)[2,2]
	pval <- coefficients(smod)[2,4]
	b_i <- coefficients(smod)[1,1]
	se_i <- coefficients(smod)[1,2]
	pval_i <- coefficients(smod)[1,4]

	testname <- "Egger regression"
 	return(list(b = b, se = se, pval = pval, nsnp = length(b_exp), b_i = b_i, se_i = se_i, pval_i = pval_i, mod = smod, dat = dat, testname=testname))
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
	pval <- pt(abs(bhat / se), df=sum(!is.na(yhat)), low=FALSE)
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
mr_eggers_regression_bootstrap <- function(b_exp, b_out, se_exp, se_out, parameters)
{
	if(length(b_exp) < 2)
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
			dat = NA,
			testname = "Egger regression (bootstrap)"
		))
	}
	nboot <- parameters$nboot
	require(reshape2)
	require(plyr)
	# Do bootstraps
	res <- array(0, c(nboot+1, 2))
	pb <- txtProgressBar(min = 0, max = nboot, initial = 0, style=3) 
	for (i in 1:nboot)
	{
		setTxtProgressBar(pb, i)
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

	return(list(b = mean(res[,2], na.rm=T), se = sd(res[,2], na.rm=T), pval = sum(sign(mean(res[,2],na.rm=T)) * res[,2] < 0)/nboot, nsnp = length(b_exp), b_i = mean(res[,1], na.rm=T), se_i = sd(res[,1], na.rm=T), pval_i = sum(sign(mean(res[,1],na.rm=T)) * res[,1] < 0)/nboot, testname="Egger regression (bootstrap)"))
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
#'         testname: Name of the test
mr_weighted_median <- function(b_exp, b_out, se_exp, se_out, parameters)
{
	if(sum(!is.na(b_exp)) < 1)
	return(list(b=NA, se=NA, pval=NA, nsnp=NA, testname="Weighted median"))

	b_iv <- b_out / b_exp
	VBj <- ((se_out)^2)/(b_exp)^2 + (b_out^2)*((se_exp^2))/(b_exp)^4
	b <- weighted_median(b_iv, 1 / VBj)
	se <- weighted_median_bootstrap(b_exp, b_out, se_exp, se_out, 1 / VBj, parameters$nboot)
	pval <- pt(abs(b/se), df = length(b_exp)-1, low=FALSE)
	return(list(b=b, se=se, pval=pval, nsnp=length(b_exp), testname="Weighted median"))
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
weighted_median_bootstrap = function(b_exp, b_out, se_exp, se_out, weights, nboot)
{
	med <- rep(0, nboot)
	for(i in 1:nboot){
		betaXG.boot = rnorm(length(b_exp), mean=b_exp, sd=se_exp)
		betaYG.boot = rnorm(length(b_out), mean=b_out, sd=se_out)
		betaIV.boot = betaYG.boot/betaXG.boot
		med[i] = weighted_median(betaIV.boot, weights)
	}
	return(sd(med)) 
}




#' Penalised weighted median MR
#'
#' Modification to standard weighted median MR
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
#'         testname: Name of the test
mr_penalised_weighted_median <- function(b_exp, b_out, se_exp, se_out, parameters)
{
	if(sum(!is.na(b_exp)) < 1)
	return(list(b=NA, se=NA, pval=NA, nsnp=NA, testname="Penalized weighted median"))
	betaIV <- b_out/b_exp # ratio estimates
	betaIVW <- sum(b_out*b_exp*se_out^-2)/sum(b_exp^2*se_out^-2) # IVW estimate
	VBj <- ((se_out)^2)/(b_exp)^2 + (b_out^2)*((se_exp^2))/(b_exp)^4
	weights <- 1/VBj
	penalty <- pchisq(weights*(betaIV-betaIVW)^2, df=1, lower.tail=FALSE)
	pen.weights <- weights*pmin(1, penalty*parameters$penk) # penalized weights
	b <- weighted_median(betaIV, pen.weights) # penalized weighted median estimate
	se <- weighted_median_bootstrap(b_exp, b_out, se_out, se_exp, pen.weights, parameters$nboot)
	pval <- pt(abs(b/se), df=length(b_exp)-1, low=FALSE)
	return(list(b = b, se = se, pval=pval, nsnp=length(b_exp), testname="Penalized weighted median"))
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
#'         testname: Name of the test
mr_ivw <- function(b_exp, b_out, se_exp, se_out, parameters)
{
	if(sum(!is.na(b_exp)) < 1)
	return(list(b=NA, se=NA, pval=NA, nsnp=NA, testname="Inverse variance weighted regression"))

	ivw.res <- summary(lm(b_out ~ -1 + b_exp, weights = 1/se_out^2))
	b <- ivw.res$coef["b_exp","Estimate"]
	se <- ivw.res$coef["b_exp","Std. Error"]/min(1,ivw.res$sigma) #sigma is the residual standard error 
	pval <- pt(abs(b/se), df = length(b_exp)-1, low=FALSE)
	Q.ivw <- ivw.res$sigma^2*(length(b_exp)-2) 
	# from formula phi =  Q/DF rearranged to to Q = phi*DF, where phi is sigma^2
	# Q.ivw<-sum((1/(se_out/b_exp)^2)*(b_out/b_exp-ivw.reg.beta)^2)
	return(list(b = b, se = se, pval = pval, nsnp=length(b_exp), testname="Inverse variance weighted regression"))
}



#' Leave one out sensitivity analysis
#'
#' @param dat Output from \code{harmonise_exposure_outcome}
#' @param method=mr_meta_fixed_simple Choose which method to use
#'
#' @export
#' @return List of data frames
mr_leaveoneout <- function(dat, method=mr_meta_fixed_simple)
{
	res <- ddply(dat, .(exposure, id.outcome), function(x)
	{
		nsnp <- nrow(x)
		if(nsnp > 1)
		{
			l <- lapply(1:nsnp, function(i)
			{
				with(x, method(beta.exposure[-i], beta.outcome[-i], se.exposure[-i], se.outcome[-i]))
			})
			l[[nsnp+1]] <- with(x, method(beta.exposure, beta.outcome, se.exposure, se.outcome))

			d <- data.frame(
				SNP = c(x$SNP, "All"),
				b = sapply(l, function(y) y$b),
				se = sapply(l, function(y) y$se),
				p = sapply(l, function(y) y$pval),
				n = x$samplesize.outcome[1]
			)

		} else {
			a <- with(x, method(beta.exposure, beta.outcome, se.exposure, se.outcome))
			d <- data.frame(
				SNP = "All",
				b = a$b,
				se = a$se,
				p = a$pval,
				n = x$samplesize.outcome[1]
			)
		}
		return(d)
	})
	return(res)
}


#' Perform 2 sample MR on each SNP individually
#'
#' @param dat Output from \code{harmonise_exposure_outcome}
#' @param method=mr_two_sample_ml Function to use for MR analysis
#'
#' @export
#' @return List of data frames
mr_singlesnp <- function(dat, parameters=default_parameters(), single_method=mr_wald_ratio, all_method=mr_two_sample_ml)
{
	res <- ddply(subset(dat, mr_keep), .(exposure, id.outcome), function(x)
	{
		nsnp <- nrow(x)
		l <- lapply(1:nsnp, function(i)
		{
			with(x, single_method(beta.exposure[i], beta.outcome[i], se.exposure[i], se.outcome[i], parameters))
		})
		l[[nsnp+1]] <- with(x, all_method(beta.exposure, beta.outcome, se.exposure, se.outcome, parameters))

		d <- data.frame(
			SNP = c(as.character(x$SNP), "All"),
			b = sapply(l, function(y) y$b),
			se = sapply(l, function(y) y$se),
			p = sapply(l, function(y) y$pval),
			Outcome.sample.size = x$samplesize.outcome[1],
			Outcome.n.case = x$ncase.outcome[1],
			Outcome.n.control = x$ncontrol.outcome[1]
		)
		d$id.outcome <- x$displayname.outcome[1]
		return(d)
	})
		res <- subset(res, select=c(exposure, id.outcome, Outcome.n.case, Outcome.n.control, Outcome.sample.size, SNP, b, se, p))
		names(res)[2] <- "Outcome"
	return(res)
}



#' Fisher's combined test
#'
#' @param pval Vector of outcome p-values
#' @export
#' @return List with the following elements:
#'         b: MR estimate
#'         se: Standard error
#'         pval: p-value
#'         testname: Name of the test
fishers_combined_test <- function(pval)
{
	p <- pchisq(-2 * sum(log(pval)), df=2*length(pval), low=FALSE)
	return(list(b=NA, se=NA, pval=p, nsnp=length(pval), testname="Fisher's combined test"))
}