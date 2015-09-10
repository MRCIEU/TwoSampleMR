
#' Perform all Mendelian randomization tests
#'
#' @param dat Harmonised exposure and outcome data. Output from \code{harmonise_exposure_outcome}
#' @param bootstrap Number of bootstraps to estimate standard error. If NULL then don't use bootstrap
#'
#' @export
#' @return List with the following elements:
#'         mr: Table of MR results
#'         extra: Table of extra results
mr <- function(dat, bootstrap=1000)
{
	require(plyr)
	res <- dlply(dat, .(outcome, exposure), function(x)
	{
		b_exp <- x$beta.exposure
		b_out <- x$beta.outcome
		se_exp <- x$se.exposure
		se_out <- x$se.outcome

		res <- list()

		res[[1]] <- mr_meta_fixed_simple(b_exp, b_out, se_exp, se_out)
		res[[2]] <- mr_meta_fixed(b_exp, b_out, se_exp, se_out)
		res[[3]] <- mr_meta_random(b_exp, b_out, se_exp, se_out)
		res[[4]] <- mr_two_sample_ml(b_exp, b_out, se_exp, se_out)
		res[[5]] <- mr_eggers_regression(b_exp, b_out, se_exp, se_out)
		res[[6]] <- mr_weighted_median(b_exp, b_out, se_exp, se_out, bootstrap)
		res[[7]] <- mr_penalised_weighted_median(b_exp, b_out, se_exp, se_out, bootstrap)
		res[[8]] <- mr_ivw(b_exp, b_out, se_exp, se_out)

		mr_tab <- data.frame(
			Exposure = x$exposure[1],
			Outcome = x$outcome[1],
			Test = sapply(res, function(x) x$testname),
			b = sapply(res, function(x) x$b),
			se = sapply(res, function(x) x$se),
			pval = sapply(res, function(x) x$pval)
		)

		extra_tab <- data.frame(
			Exposure = x$exposure[1],
			Outcome = x$outcome[1],
			Test = c("Egger regression intercept"),
			b = res[[5]]$b_i,
			se = res[[5]]$se_i,
			pval = res[[5]]$pval_i
		)

		l <- list(mr_tab=mr_tab, extra_tab=extra_tab)
		return(l)
	})
	mr_tab <- rbind.fill(lapply(res, function(x) x$mr_tab))
	extra_tab <- rbind.fill(lapply(res, function(x) x$extra_tab))
	return(list(mr=mr_tab, extra=extra_tab))
}

mr_sensitivity_analysis <- function(dat)
{
	res <- dlply(dat, .(outcome, exposure), function(x)
	{
		b_exp <- x$beta.exposure
		b_out <- x$beta.outcome
		se_exp <- x$se.exposure
		se_out <- x$se.outcome

		res <- list()

		res[[1]] <- mr_leaveoneout(b_exp, b_out, se_exp, se_out)

	})
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
mr_meta_fixed_simple <- function(b_exp, b_out, se_exp, se_out)
{
	b <- sum(b_exp*b_out / se_out^2) / sum(b_exp^2/se_out^2)
	se <- sqrt(1 / sum(b_exp^2/se_out^2))
	pval <- pt(abs(b) / se, df = length(b_exp)-1, low=FALSE)
	return(list(b=b, se=se, pval=pval, testname="Fixed effects meta analysis (simple SE)"))
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
mr_meta_fixed <- function(b_exp, b_out, se_exp, se_out, Cov=0)
{
	require(meta)
	ratio <- b_out / b_exp
	ratio.se <- sqrt((se_out^2/b_exp^2) + (b_out^2/b_exp^4)*se_exp^2 - 2*(b_out/b_exp^3)*Cov)
	res <- metagen(ratio, ratio.se)
	b <- res$TE.fixed
	se <- res$seTE.fixed
	pval <- res$pval.fixed
	return(list(b=b, se=se, pval=pval, testname="Fixed effects meta analysis (delta method)"))
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
mr_meta_random <- function(b_exp, b_out, se_exp, se_out, Cov=0)
{
	require(meta)
	ratio <- b_out / b_exp
	ratio.se <- sqrt((se_out^2/b_exp^2) + (b_out^2/b_exp^4)*se_exp^2 - 2*(b_out/b_exp^3)*Cov)
	res <- metagen(ratio, ratio.se)
	b <- res$TE.random
	se <- res$seTE.random
	pval <- res$pval.random
	return(list(b=b, se=se, pval=pval, testname="Random effects meta analysis (delta method)"))
}



mr_two_sample_ml <- function(b_exp, b_out, se_exp, se_out)
{
	loglikelihood <- function(param) {
		return(1/2*sum((b_exp-param[1:length(b_exp)])^2/se_exp^2)+1/2*sum((b_out-param[length(b_exp)+1]*param[1:length(b_exp)])^2/se_out^2))
	}
	opt <- optim(
		c(b_exp, sum(b_exp*b_out/se_out^2)/sum(b_exp^2/se_out^2)),
		loglikelihood, 
		hessian=TRUE, 
		control = list(maxit=25000))

	b <- opt$par[length(b_exp)+1]
	se <- sqrt(solve(opt$hessian)[length(b_exp)+1,length(b_exp)+1])
	pval <- pt(abs(b) / se, df = length(b_exp)-1, low=FALSE)
	return(list(b=b, se=se, pval=pval, testname="Wald ratio (random)"))
}



#' Egger's regression for Mendelian randomization
#'
#' @param b_exp Vector of genetic effects on exposure
#' @param b_out Vector of genetic effects on outcome
#' @param se_exp Standard errors of genetic effects on exposure
#' @param se_out Standard errors of genetic effects on outcome
#' @param bootstrap Number of bootstraps to estimate standard error. If NULL then don't use bootstrap
#' @param alpha Quantiles to use for calculating confidence intervals using bootstraps. Default 0.05
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
mr_eggers_regression <- function(b_exp, b_out, se_exp, se_out, bootstrap=NULL, alpha=0.05)
{
	stopifnot(length(b_exp) == length(b_out))
	stopifnot(length(se_exp) == length(se_out))
	stopifnot(length(b_exp) == length(se_out))

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

	if(!is.null(bootstrap))
	{
		boots <- eggers_regression_bootstrap(b_exp, b_out, se_exp, se_out, bootstrap)
		se <- boots$boots$se[boots$boots$param == "slope" & boots$boots$stat == "b" & boots$boots$what == "bootstrap"]
		se_i <- boots$boots$se[boots$boots$param == "intercept" & boots$boots$stat == "b" & boots$boots$what == "bootstrap"]
		testname <- "Egger regression (bootstrapped SE)"
	}
	return(list(b = b, se = se, pval = pval, b_i = b_i, se_i = se_i, pval_i = pval_i, mod = smod, dat = dat, testname=testname))
}



#' Run bootstrap to generate standard errors for MR
#'
#' @param b_exp Vector of genetic effects on exposure
#' @param b_out Vector of genetic effects on outcome
#' @param se_exp Standard errors of genetic effects on exposure
#' @param se_out Standard errors of genetic effects on outcome
#' @param nboot Number of bootstraps. Default 1000
#' @param alpha Quantile to use for confidence intervals. Default 0.05
#'
#' @export
#' @return A list with the following elements:
#'         boots: data frame of bootstrap results
#'         res: data frame of summary results from each bootstrap
#'         data: data frame of input summary stats
eggers_regression_bootstrap <- function(b_exp, b_out, se_exp, se_out, nboot = 1000, alpha = 0.05)
{
	require(reshape2)
	require(plyr)
	# Do bootstraps
	res <- array(0, c(n+1, 4))
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
		r <- summary(lm(ys ~ xs, weights=1/se_out^2))

		#collect coefficient from given line.
		res[i, 1] <- r$coefficients[1,1]
		res[i, 2] <- r$coefficients[1,2]
		res[i, 3] <- r$coefficients[2,1]
		res[i, 4] <- r$coefficients[2,2]
	}
	cat("\n")

	# Run original analysis
	b_out <- b_out*sign(b_exp)
	b_exp <- abs(b_exp)
	dat <- data.frame(b_out, b_exp, se_out, se_exp)
	orig <- coefficients(summary(lm(b_out ~ b_exp, weights=1/se_out^2)))
	res[n+1, ] <- c(orig[1,1], orig[1,2], orig[2,1], orig[2,2])
	res <- as.data.frame(res)
	res$what <- "bootstrap"
	res$what[n+1] <- "original"

	datl <- melt(res, measure.vars=c("V1", "V2", "V3", "V4"))
	datl$param <- "slope"
	datl$param[datl$variable %in% c("V1", "V2")] <- "intercept"
	datl$stat <- "b"
	datl$stat[datl$variable %in% c("V2", "V4")] <- "se"

	qu <- ddply(datl, .(param, stat, what), summarise, 
		m=mean(value),
		se=sd(value),
		qupper=quantile(value, alpha),
		qlower=quantile(value, 1-alpha),
		pval=sum(value < 0)/length(value))

	res <- as.data.frame(res)
	names(res) <- c("b_i", "se_i", "b", "se")

	return(list(boots=qu, res=res, data=dat))
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
mr_weighted_median <- function(b_exp, b_out, se_exp, se_out, nboot=1000)
{
	b_iv <- b_out / b_exp
	VBj <- ((se_out)^2)/(b_exp)^2 + (b_out^2)*((se_exp^2))/(b_exp)^4
	b <- weighted_median(b_iv, 1 / VBj)
	se <- weighted_median_bootstrap(b_exp, b_out, se_exp, se_out, 1 / VBj)
	pval <- pt(abs(b/se), df = length(b_exp)-1, low=FALSE)
	return(list(b=b, se=se, pval=pval, testname="Weighted median"))
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
weighted_median_bootstrap = function(b_exp, b_out, se_exp, se_out, weights, nboot=1000)
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
mr_penalised_weighted_median <- function(b_exp, b_out, se_exp, se_out, penk=20, nboot=1000)
{
	betaIV <- b_out/b_exp # ratio estimates
	betaIVW <- sum(b_out*b_exp*se_out^-2)/sum(b_exp^2*se_out^-2) # IVW estimate
	VBj <- ((se_out)^2)/(b_exp)^2 + (b_out^2)*((se_exp^2))/(b_exp)^4
	weights <- 1/VBj
	penalty <- pchisq(weights*(betaIV-betaIVW)^2, df=1, lower.tail=FALSE)
	pen.weights <- weights*pmin(1, penalty*penk) # penalized weights
	b <- weighted_median(betaIV, pen.weights) # penalized weighted median estimate
	se <- weighted_median_bootstrap(b_exp, b_out, se_out, se_exp, pen.weights)
	pval <- pt(abs(b/se), df=length(b_exp)-1, low=FALSE)
	return(list(b = b, se = se, pval=pval, testname="Penalized weighted median"))
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
mr_ivw <- function(b_exp, b_out, se_exp, se_out)
{
	ivw.res <- summary(lm(b_out ~ -1 + b_exp, weights = 1/se_out^2))
	b <- ivw.res$coef["b_exp","Estimate"]
	se <- ivw.res$coef["b_exp","Std. Error"]/min(1,ivw.res$sigma) #sigma is the residual standard error 
	pval <- pt(abs(b/se), df = length(b_exp)-1, low=FALSE)
	Q.ivw <- ivw.res$sigma^2*(length(b_exp)-2) 
	# from formula phi =  Q/DF rearranged to to Q = phi*DF, where phi is sigma^2
	# Q.ivw<-sum((1/(se_out/b_exp)^2)*(b_out/b_exp-ivw.reg.beta)^2)
	return(list(b = b, se = se, pval = pval, testname="Inverse variance weighted regression"))
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
	res <- dlply(dat, .(exposure, outcome), function(x)
	{
		nsnp <- nrow(x)
		l <- lapply(1:nsnp, function(i)
		{
			with(x, method(beta.exposure[-i], beta.outcome[-i], se.exposure[-i], se.outcome[-i]))
		})

		d <- data.frame(
			SNP = x$SNP,
			b = sapply(l, function(x) x$b),
			se = sapply(l, function(x) x$se)
		)
		return(d)
	})
	return(res)
}

