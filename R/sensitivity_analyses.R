#' Leave one out sensitivity analysis
#'
#' @param dat Output from \code{harmonise_exposure_outcome}
#' @param method=mr_ivw Choose which method to use
#'
#' @export
#' @return List of data frames
mr_leaveoneout <- function(dat, parameters=default_parameters(), method=mr_ivw)
{
	if(!"samplesize.outcome" %in% names(dat))
	{
		dat$samplesize.outcome <- NA
	}

	stopifnot("outcome" %in% names(dat))
	stopifnot("exposure" %in% names(dat))
	stopifnot("beta.exposure" %in% names(dat))
	stopifnot("beta.outcome" %in% names(dat))
	stopifnot("se.exposure" %in% names(dat))
	stopifnot("se.outcome" %in% names(dat))


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

	if(!"samplesize.outcome" %in% names(dat))
	{
		dat$samplesize.outcome <- NA
	}

	stopifnot("outcome" %in% names(dat))
	stopifnot("exposure" %in% names(dat))
	stopifnot("beta.exposure" %in% names(dat))
	stopifnot("beta.outcome" %in% names(dat))
	stopifnot("se.exposure" %in% names(dat))
	stopifnot("se.outcome" %in% names(dat))

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
	warning("Prior to version 0.4.9 there was a bug in the IVW Q statistic estimate, leading to a slight underestimation in heterogeneity. This has now been resolved.")
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


