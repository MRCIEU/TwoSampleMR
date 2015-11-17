
#' Perform all Mendelian randomization tests
#'
#' @param dat Harmonised exposure and outcome data. Output from \code{harmonise_exposure_outcome}
#' @param nboot Number of bootstraps to estimate standard error. Default is 1000.
#' @param method_list List of methods to use in analysis. See \code{mr_method_list()} for details.
#'
#' @export
#' @return List with the following elements:
#'         mr: Table of MR results
#'         extra: Table of extra results
mr <- function(dat, nboot=1000, method_list=mr_method_list())
{
	require(plyr)
	res <-dlply(dat, .(consortium,outcome, exposure), function(x)
	# x<-dlply(dat, .(outcome, exposure))
	{
		x <- mutate(x)

		keep_mr <- rep(TRUE, nrow(x))
		keep_pval <- rep(FALSE, nrow(x))
		keep_mr[!is.finite(x$beta.exposure)] <- FALSE #
		keep_mr[!is.finite(x$beta.outcome)] <- FALSE #
		keep_mr[!is.finite(x$se.exposure)] <- FALSE #
		keep_mr[!is.finite(x$se.outcome)] <- FALSE #
		keep_pval[is.finite(x$pval.outcome) & x$pval.outcome > 0 & x$pval.outcome <= 1] <- TRUE

		print(sum(keep_mr))
		print(sum(keep_pval))

		if(sum(keep_mr) == 0) return(NULL)

		b_exp <- x$beta.exposure[keep_mr]
		b_out <- x$beta.outcome[keep_mr]
		se_exp <- x$se.exposure[keep_mr]
		se_out <- x$se.outcome[keep_mr]
		p_out <- x$pval.outcome[keep_pval]

		res <- lapply(method_list, function(meth)
		{
			if(meth %in% c("mr_weighted_median", "mr_penalised_weighted_median", "mr_eggers_regression_bootstrap"))
			{
				get(meth)(b_exp, b_out, se_exp, se_out, nboot)
			}
			else if(meth == "fishers_combined_test")
			{
				get(meth)(p_out)
			} else {
				get(meth)(b_exp, b_out, se_exp, se_out)	
			}
		})
		mr_tab <- data.frame(
			consortium = x$consortium[1],
			Exposure = x$exposure[1],
			Outcome = x$outcome[1],
			Test = sapply(res, function(x) x$testname),
			"n SNPs" = sapply(res, function(x) x$nsnp),
			b = sapply(res, function(x) x$b),
			se = sapply(res, function(x) x$se),
			pval = sapply(res, function(x) x$pval)
		)

		# mregger <- res[[which(method_list == "mr_eggers_regression")[1]]]
		# mreggerb <- res[[which(method_list == "mr_eggers_regression_bootstrap")[1]]]
		# if(is.null(mregger)) mregger <- mr_eggers_regression(b_exp, b_out, se_exp, se_out)
		# if(is.null(mreggerb)) mreggerb <- mr_eggers_regression_bootstrap(b_exp, b_out, se_exp, se_out, nboot)
		# extra_tab <- data.frame(
		# 	Exposure = x$exposure[1],
		# 	Outcome = x$outcome[1],
		# 	Test = c("Egger regression intercept", "... Using bootstrap"),
		# 	"n SNPs" = c(mregger$nsnp, mreggerb$nsnp),
		# 	b = c(mregger$b_i, mreggerb$b_i),
		# 	se = c(mregger$se_i, mreggerb$se_i),
		# 	pval = c(mregger$pval_i, mreggerb$pval_i)
		# )

		# l <- list(mr_tab=mr_tab, extra_tab=extra_tab)
		# return(l)
		return(mr_tab)
	})
	mr_tab <- rbind.fill(lapply(res, function(x) x))
	# mr_tab <- rbind.fill(lapply(res, function(x) x$mr_tab))
# 	# extra_tab <- rbind.fill(lapply(res, function(x) x$extra_tab))
# 	# return(list(mr=mr_tab, extra=extra_tab))
# 	return(mr_tab)
}


#' Get list of available MR methods
#'
#' @export
#' @return character vector of method names
mr_method_list <- function()
{
	c(
		"mr_meta_fixed_simple",
		"mr_meta_fixed",
		"mr_meta_random",
		"mr_two_sample_ml",
		"mr_eggers_regression",
		"mr_eggers_regression_bootstrap",
		"mr_weighted_median",
		"mr_penalised_weighted_median",
		"mr_ivw",
		"fishers_combined_test"
	)
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
mr_meta_fixed_simple <- function(b_exp, b_out, se_exp, se_out, ...)
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
mr_meta_fixed <- function(b_exp, b_out, se_exp, se_out, Cov=0, ...)
{
	require(meta)
	ratio <- b_out / b_exp
	ratio.se <- sqrt((se_out^2/b_exp^2) + (b_out^2/b_exp^4)*se_exp^2 - 2*(b_out/b_exp^3)*Cov)
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
mr_meta_random <- function(b_exp, b_out, se_exp, se_out, Cov=0, ...)
{
	require(meta)
	ratio <- b_out / b_exp
	ratio.se <- sqrt((se_out^2/b_exp^2) + (b_out^2/b_exp^4)*se_exp^2 - 2*(b_out/b_exp^3)*Cov)
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
mr_two_sample_ml <- function(b_exp, b_out, se_exp, se_out, ...)
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
mr_eggers_regression <- function(b_exp, b_out, se_exp, se_out, bootstrap=NULL, ...)
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
mr_eggers_regression_bootstrap <- function(b_exp, b_out, se_exp, se_out, nboot = 1000, ...)
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
mr_weighted_median <- function(b_exp, b_out, se_exp, se_out, nboot=1000, ...)
{
	if(sum(!is.na(b_exp)) < 1)
	return(list(b=NA, se=NA, pval=NA, nsnp=NA, testname="Weighted median"))

	b_iv <- b_out / b_exp
	VBj <- ((se_out)^2)/(b_exp)^2 + (b_out^2)*((se_exp^2))/(b_exp)^4
	b <- weighted_median(b_iv, 1 / VBj)
	se <- weighted_median_bootstrap(b_exp, b_out, se_exp, se_out, 1 / VBj)
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
weighted_median_bootstrap = function(b_exp, b_out, se_exp, se_out, weights, nboot=1000, ...)
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
	if(sum(!is.na(b_exp)) < 1)
	return(list(b=NA, se=NA, pval=NA, nsnp=NA, testname="Penalized weighted median"))

	betaIV <- b_out/b_exp # ratio estimates
	betaIVW <- sum(b_out*b_exp*se_out^-2)/sum(b_exp^2*se_out^-2) # IVW estimate
	VBj <- ((se_out)^2)/(b_exp)^2 + (b_out^2)*((se_exp^2))/(b_exp)^4
	weights <- 1/VBj
	penalty <- pchisq(weights*(betaIV-betaIVW)^2, df=1, lower.tail=FALSE)
	pen.weights <- weights*pmin(1, penalty*penk) # penalized weights
	b <- weighted_median(betaIV, pen.weights) # penalized weighted median estimate
	se <- weighted_median_bootstrap(b_exp, b_out, se_out, se_exp, pen.weights)
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
mr_ivw <- function(b_exp, b_out, se_exp, se_out, ...)
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
	res <- ddply(dat, .(exposure, outcome), function(x)
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
mr_singlesnp <- function(dat, method=mr_meta_fixed)
{
	res <- ddply(dat, .(exposure, outcome), function(x)
	{
		nsnp <- nrow(x)
		l <- lapply(1:nsnp, function(i)
		{
			with(x, method(beta.exposure[i], beta.outcome[i], se.exposure[i], se.outcome[i]))
		})
		l[[nsnp+1]] <- with(x, method(beta.exposure, beta.outcome, se.exposure, se.outcome))

		d <- data.frame(
			SNP = c(x$SNP, "All"),
			b = sapply(l, function(y) y$b),
			se = sapply(l, function(y) y$se),
			p = sapply(l, function(y) y$pval),
			n = x$samplesize.outcome[1]
		)
		return(d)
	})
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
fishers_combined_test <- function(pval, ...)
{
	p <- pchisq(-2 * sum(log(pval)), df=2*length(pval), low=FALSE)
	return(list(b=NA, se=NA, pval=p, nsnp=length(pval), testname="Fisher's combined test"))
}