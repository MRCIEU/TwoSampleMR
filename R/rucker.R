#' I-square calculation
#'
#' This function calculates the I-squared_GX statistic.
#' 
#' @param y Vector of effects.
#' @param s Vector of standard errors.
#'
#' @export
#' @return Isq value
Isq <- function(y,s)
{
	k <- length(y)
	w <- 1/s^2
	sum.w <- sum(w)
	mu.hat <- sum(y*w)/sum.w
	Q <- sum(w*(y-mu.hat)^2)
	Isq <- (Q - (k-1))/Q
	Isq <- max(0,Isq)
	return(Isq)
}

#' @importFrom stats qchisq qnorm
PM <- function(y = y, s = s, Alpha = 0.1)
{
	k = length(y)
	df = k - 1
	sig = qnorm(1-Alpha/2)
	low = qchisq((Alpha/2), df)
	up = qchisq(1-(Alpha/2), df)
	med = qchisq(0.5, df)
	mn = df
	mode = df-1
	Quant = c(low, mode, mn, med, up)
	L = length(Quant)
	Tausq = NULL
	Isq = NULL
	CI = matrix(nrow = L, ncol = 2)
	MU = NULL
	v = 1/s^2
	sum.v = sum(v)
	typS = sum(v*(k-1))/(sum.v^2 - sum(v^2))
	for(j in 1:L)
	{
		tausq = 0 ; F = 1 ;TAUsq = NULL
		while(F>0)
		{
			TAUsq = c(TAUsq, tausq)
			w = 1/(s^2+tausq)
			sum.w = sum(w)
			w2 = w^2
			yW = sum(y*w)/sum.w
			Q1 = sum(w*(y-yW)^2)
			Q2 = sum(w2*(y-yW)^2)
			F = Q1-Quant[j]
			Ftau = max(F,0)
			delta = F/Q2
			tausq = tausq + delta
		}
		MU[j] = yW
		V = 1/sum(w)
		Tausq[j] = max(tausq,0)
		Isq[j] = Tausq[j]/(Tausq[j]+typS)
		CI[j,] = yW + sig*c(-1,1) *sqrt(V)
	}
	return(list(tausq = Tausq, muhat = MU, Isq = Isq, CI = CI, quant = Quant))
}


#' MR Rucker framework
#'
#' MR Rucker framework.
#'
#' @md
#' @param dat Output from [`harmonise_data()`].
#' @param parameters List of Qthresh for determing transition between models, and alpha values for calculating confidence intervals. Defaults to 0.05 for both in `default_parameters()`.
#'
#' @export
#' @return list
mr_rucker <- function(dat, parameters=default_parameters())
{
	dat <- subset(dat, mr_keep)
	d <- subset(dat, !duplicated(paste(id.exposure, " - ", id.outcome)), select=c(exposure, outcome, id.exposure, id.outcome))
	res <- list()
	attributes(res)$id.exposure <- d$id.exposure
	attributes(res)$id.outcome <- d$id.outcome
	attributes(res)$exposure <- d$exposure
	attributes(res)$outcome <- d$outcome
	for(j in 1:nrow(d))
	{
		x <- subset(dat, exposure == d$exposure[j] & outcome == d$outcome[j])
		message(x$exposure[1], " - ", x$outcome[1])
		res[[j]] <- mr_rucker_internal(x, parameters)
	}
	return(res)
}

#' @importFrom stats coefficients cooks.distance lm mad pchisq pnorm pt qnorm
mr_rucker_internal <- function(dat, parameters=default_parameters())
{
	if("mr_keep" %in% names(dat)) dat <- subset(dat, mr_keep)

	if(nrow(dat) < 3) 
	{
		warning("Need at least 3 SNPs")
		return(NULL)
	}


    sign0 <- function(x) {
		x[x == 0] <- 1
		return(sign(x))
    }
	dat$beta.outcome <- dat$beta.outcome * sign0(dat$beta.exposure)
	dat$beta.exposure <- abs(dat$beta.exposure)

	Qthresh <- parameters$Qthresh
	alpha <- parameters$alpha

	nsnp <- nrow(dat)
	b_exp <- dat$beta.exposure
	b_out <- dat$beta.outcome
	se_exp <- dat$se.exposure
	se_out <- dat$se.outcome
	w <- b_exp^2 / se_out^2
	y <- b_out / se_out
	x <- b_exp / se_out
	i <- 1 / se_out


	# IVW FE
	lmod_ivw <- lm(y ~ 0 + x)
	mod_ivw <- summary(lmod_ivw)
	b_ivw_fe <- coefficients(mod_ivw)[1,1]

	# Q_ivw <- sum((y - x*b_ivw_fe)^2)
	Q_ivw <- mod_ivw$sigma^2 * (nsnp - 1)
	Q_df_ivw <- length(b_exp) - 1
	Q_pval_ivw <- pchisq(Q_ivw, Q_df_ivw, lower.tail = FALSE)
	phi_ivw <- Q_ivw / (nsnp - 1)

	se_ivw_fe <- coefficients(mod_ivw)[1,2] / max(mod_ivw$sigma, 1)
	if(parameters$test_dist == "z")
	{
		pval_ivw_fe <- pnorm(abs(b_ivw_fe/se_ivw_fe), lower.tail=FALSE) * 2
	} else {
		pval_ivw_fe <- pt(abs(b_ivw_fe/se_ivw_fe), nsnp-1, lower.tail=FALSE) * 2
	}


	# IVW MRE
	b_ivw_re <- b_ivw_fe
	# se_ivw_re <- sqrt(phi_ivw / sum(w))
	se_ivw_re <- coefficients(mod_ivw)[1,2]
	# pval_ivw_re <- pt(abs(b_ivw_re/se_ivw_re), nsnp-1, lower.tail=FALSE) * 2
	if(parameters$test_dist == "z")
	{
		pval_ivw_re <- pnorm(abs(coefficients(mod_ivw)[1,1]/coefficients(mod_ivw)[1,2]), lower.tail=FALSE) * 2
	} else {
		pval_ivw_re <- coefficients(mod_ivw)[1,4]
	}


	# Egger FE
	lmod_egger <- lm(y ~ 0 + i + x)
	mod_egger <- summary(lmod_egger)

	b1_egger_fe <- coefficients(mod_egger)[2,1]
	b0_egger_fe <- coefficients(mod_egger)[1,1]

	# This is equivalent to mod$sigma^2
	# Q_egger <- sum(
	# 	1 / se_out^2 * (b_out - (b0_egger_fe + b1_egger_fe * b_exp))^2
	# )
	Q_egger <- mod_egger$sigma^2 * (nsnp - 2)
	Q_df_egger <- nsnp - 2
	Q_pval_egger <- pchisq(Q_egger, Q_df_egger, lower.tail=FALSE)
	phi_egger <- Q_egger / (nsnp - 2)

	se1_egger_fe <- coefficients(mod_egger)[2,2] / max(mod_egger$sigma, 1)
	pval1_egger_fe <- pt(abs(b1_egger_fe/se1_egger_fe), nsnp-2, lower.tail=FALSE) * 2
	se0_egger_fe <- coefficients(mod_egger)[1,2] / max(mod_egger$sigma, 1)
	if(parameters$test_dist == "z")
	{
		pval0_egger_fe <- pnorm(abs(b0_egger_fe/se0_egger_fe), lower.tail=FALSE) * 2	
	} else {
		pval0_egger_fe <- pt(abs(b0_egger_fe/se0_egger_fe), nsnp-2, lower.tail=FALSE) * 2
	}

	# Egger RE
	b1_egger_re <- coefficients(mod_egger)[2,1]
	se1_egger_re <- coefficients(mod_egger)[2,2]
	pval1_egger_re <- coefficients(mod_egger)[2,4]
	b0_egger_re <- coefficients(mod_egger)[1,1]
	se0_egger_re <- coefficients(mod_egger)[1,2]
	if(parameters$test_dist == "z")
	{
		pval0_egger_re <- pnorm(coefficients(mod_egger)[1,1]/coefficients(mod_egger)[1,2], lower.tail=FALSE)
	} else {
		pval0_egger_re <- coefficients(mod_egger)[1,4]
	}

	results <- data.frame(
		Method = c("IVW fixed effects", "IVW random effects", "Egger fixed effects", "Egger random effects"),
		nsnp = nsnp,
		Estimate = c(b_ivw_fe, b_ivw_re, b1_egger_fe, b1_egger_re),
		SE = c(se_ivw_fe, se_ivw_re, se1_egger_fe, se1_egger_re)
	)
	results$CI_low <- results$Estimate - qnorm(1-alpha/2) * results$SE
	results$CI_upp <- results$Estimate + qnorm(1-alpha/2) * results$SE
	results$P <- c(pval_ivw_fe, pval_ivw_re, pval1_egger_fe, pval1_egger_re)

	Qdiff <- max(0, Q_ivw - Q_egger)
	Qdiff_p <- pchisq(Qdiff, 1, lower.tail=FALSE)


	Q <- data.frame(
		Method=c("Q_ivw", "Q_egger", "Q_diff"),
		Q=c(Q_ivw, Q_egger, Qdiff),
		df=c(Q_df_ivw, Q_df_egger, 1),
		P=c(Q_pval_ivw, Q_pval_egger, Qdiff_p)
	)

	intercept <- data.frame(
		Method=c("Egger fixed effects", "Egger random effects"),
		Estimate = c(b0_egger_fe, b0_egger_fe),
		SE = c(se0_egger_fe, se0_egger_re)
	)
	intercept$CI_low <- intercept$Estimate - qnorm(1-alpha/2) * intercept$SE
	intercept$CI_upp <- intercept$Estimate + qnorm(1-alpha/2) * intercept$SE
	intercept$P <- c(pval0_egger_fe, pval0_egger_re)

	if(Q_pval_ivw <= Qthresh)
	{
		if(Qdiff_p <= Qthresh)
		{
			if(Q_pval_egger <= Qthresh)
			{
				res <- "D"
			} else {
				res <- "C"
			}
		} else {
			res <- "B"
		}
	} else {
		res <- "A"
	}

	selected <- results[c("A", "B", "C", "D") %in% res, ]
	selected$Method <- "Rucker"

	if(res %in% c("A", "B"))
	{
		cd <- cooks.distance(lmod_ivw)
	} else {
		cd <- cooks.distance(lmod_egger)
	}

	return(list(rucker=results, intercept=intercept, Q=Q, res=res, selected=selected, cooksdistance=cd, lmod_ivw=lmod_ivw, lmod_egger=lmod_egger))
}



#' Run rucker with bootstrap estimates
#'
#' Run Rucker with bootstrap estimates.
#'
#' @md
#' @param dat Output from [`harmonise_data`].
#' @param parameters List of parameters. The default is `default_parameters()`.
#'
#' @return List
#' @export
#' @importFrom stats median pt qchisq qnorm quantile rnorm sd
mr_rucker_bootstrap <- function(dat, parameters=default_parameters())
{
	requireNamespace("ggplot2", quietly=TRUE)
	requireNamespace("plyr", quietly=TRUE)

	if("mr_keep" %in% names(dat)) dat <- subset(dat, mr_keep)

	nboot <- parameters$nboot
	nsnp <- nrow(dat)
	Qthresh <- parameters$Qthresh

	# Main result
	rucker <- mr_rucker(dat, parameters)
	dat2 <- dat
	l <- list()
	for(i in 1:nboot)
	{
		dat2$beta.exposure <- rnorm(nsnp, mean=dat$beta.exposure, sd=dat$se.exposure)
		dat2$beta.outcome <- rnorm(nsnp, mean=dat$beta.outcome, sd=dat$se.outcome)
		l[[i]] <- mr_rucker(dat2, parameters)
	}

	modsel <- plyr::rbind.fill(lapply(l, function(x) x$selected))
	modsel$model <- sapply(l, function(x) x$res)

	bootstrap <- data.frame(
		Q = c(rucker$Q$Q[1], sapply(l, function(x) x$Q$Q[1])),
		Qdash = c(rucker$Q$Q[2], sapply(l, function(x) x$Q$Q[2])),
		model = c(rucker$res, sapply(l, function(x) x$res)),
		i = c("Full", rep("Bootstrap", nboot))
	)

	# Get the median estimate
	rucker_point <- rucker$selected
	rucker_point$Method <- "Rucker point estimate"

	rucker_median <- data.frame(
		Method = "Rucker median",
		nsnp = nsnp,
		Estimate = median(modsel$Estimate),
		SE = mad(modsel$Estimate),
		CI_low = quantile(modsel$Estimate, 0.025),
		CI_upp = quantile(modsel$Estimate, 0.975)
	)
	rucker_median$P <- 2 * pt(abs(rucker_median$Estimate/rucker_median$SE), nsnp-1, lower.tail=FALSE)

	rucker_mean <- data.frame(
		Method = "Rucker mean",
		nsnp = nsnp,
		Estimate = mean(modsel$Estimate),
		SE = sd(modsel$Estimate)
	)
	rucker_mean$CI_low <- rucker_mean$Estimate - qnorm(Qthresh/2, lower.tail=TRUE) * rucker_mean$SE
	rucker_mean$CI_upp <- rucker_mean$Estimate + qnorm(Qthresh/2, lower.tail=TRUE) * rucker_mean$SE
	rucker_mean$P <- 2 * pt(abs(rucker_mean$Estimate/rucker_mean$SE), nsnp-1, lower.tail=FALSE)


	res <- rbind(rucker$rucker, rucker_point, rucker_mean, rucker_median)
	rownames(res) <- NULL

	p1 <- ggplot2::ggplot(bootstrap, ggplot2::aes_string(x="Q", y="Qdash")) +
		ggplot2::geom_point(ggplot2::aes_string(colour="model")) +
		ggplot2::geom_point(data=subset(bootstrap, i=="Full")) +
		ggplot2::scale_colour_brewer(type="qual") +
		ggplot2::xlim(0, max(bootstrap$Q, bootstrap$Qdash)) +
		ggplot2::ylim(0, max(bootstrap$Q, bootstrap$Qdash)) +
		ggplot2::geom_abline(slope=1, colour="grey") +
		ggplot2::geom_abline(slope=1, intercept=-qchisq(Qthresh, 1, lower.tail=FALSE), linetype="dotted") +
		ggplot2::geom_hline(yintercept = qchisq(Qthresh, nsnp - 2, lower.tail=FALSE), linetype="dotted") +
		ggplot2::geom_vline(xintercept = qchisq(Qthresh, nsnp - 1, lower.tail=FALSE), linetype="dotted") +
		ggplot2::labs(x="Q", y="Q'")

	modsel$model_name <- "IVW"
	modsel$model_name[modsel$model %in% c("C", "D")] <- "Egger"

	p2 <- ggplot2::ggplot(modsel, ggplot2::aes_string(x="Estimate")) +
		ggplot2::geom_density(ggplot2::aes_string(fill="model_name"), alpha=0.4) +
		ggplot2::geom_vline(data=res, ggplot2::aes_string(xintercept="Estimate", colour="Method")) +
		ggplot2::scale_colour_brewer(type="qual") +
		ggplot2::scale_fill_brewer(type="qual") + 
		ggplot2::labs(fill="Bootstrap estimates", colour="")

	return(list(rucker=rucker, res=res, bootstrap_estimates=modsel, boostrap_q=bootstrap, q_plot=p1, e_plot=p2))
}


#' Run rucker with jackknife estimates
#'
#' Run rucker with jackknife estimates.
#'
#' @md
#' @param dat Output from harmonise_data.
#' @param parameters List of parameters. The default is `default_parameters()`.
#'
#' @export
#' @return List
mr_rucker_jackknife <- function(dat, parameters=default_parameters())
{
	dat <- subset(dat, mr_keep)
	d <- subset(dat, !duplicated(paste(id.exposure, " - ", id.outcome)), select=c(exposure, outcome, id.exposure, id.outcome))
	res <- list()
	attributes(res)$id.exposure <- d$id.exposure
	attributes(res)$id.outcome <- d$id.outcome
	attributes(res)$exposure <- d$exposure
	attributes(res)$outcome <- d$outcome
	for(j in 1:nrow(d))
	{
		x <- subset(dat, exposure == d$exposure[j] & outcome == d$outcome[j])
		message(x$exposure[1], " - ", x$outcome[1])
		res[[j]] <- mr_rucker_jackknife_internal(x, parameters)
	}
	return(res)
}

#' @importFrom stats mad median pt qchisq qnorm quantile sd
mr_rucker_jackknife_internal <- function(dat, parameters=default_parameters())
{
	requireNamespace("ggplot2", quietly=TRUE)

	if("mr_keep" %in% names(dat)) dat <- subset(dat, mr_keep)

	nboot <- parameters$nboot
	nsnp <- nrow(dat)
	Qthresh <- parameters$Qthresh


	# Main result
	rucker <- mr_rucker_internal(dat, parameters)
	rucker_point <- rucker$selected
	rucker_point$Method <- "Rucker point estimate"


	if(nrow(dat) < 15)
	{
		message("Too few SNPs for jackknife")
		res <- rbind(rucker$rucker, rucker_point)
		return(list(rucker=rucker, res=res, bootstrap_estimates=NULL, boostrap_q=NULL, q_plot=NULL, e_plot=NULL))

	} else {

		l <- list()
		for(i in 1:nboot)
		{
			# dat2$beta.exposure <- rnorm(nsnp, mean=dat$beta.exposure, sd=dat$se.exposure)
			# dat2$beta.outcome <- rnorm(nsnp, mean=dat$beta.outcome, sd=dat$se.outcome)
			dat2 <- dat[sample(1:nrow(dat), nrow(dat), replace=TRUE), ]
			l[[i]] <- mr_rucker_internal(dat2, parameters)
		}

		modsel <- plyr::rbind.fill(lapply(l, function(x) x$selected))
		modsel$model <- sapply(l, function(x) x$res)

		bootstrap <- data.frame(
			Q = c(rucker$Q$Q[1], sapply(l, function(x) x$Q$Q[1])),
			Qdash = c(rucker$Q$Q[2], sapply(l, function(x) x$Q$Q[2])),
			model = c(rucker$res, sapply(l, function(x) x$res)),
			i = c("Full", rep("Jackknife", nboot))
		)

		# Get the median estimate

		rucker_median <- data.frame(
			Method = "Rucker median (JK)",
			nsnp = nsnp,
			Estimate = median(modsel$Estimate),
			SE = mad(modsel$Estimate),
			CI_low = quantile(modsel$Estimate, 0.025),
			CI_upp = quantile(modsel$Estimate, 0.975)
		)
		rucker_median$P <- 2 * pt(abs(rucker_median$Estimate/rucker_median$SE), nsnp-1, lower.tail=FALSE)

		rucker_mean <- data.frame(
			Method = "Rucker mean (JK)",
			nsnp = nsnp,
			Estimate = mean(modsel$Estimate),
			SE = sd(modsel$Estimate)
		)
		rucker_mean$CI_low <- rucker_mean$Estimate - qnorm(Qthresh/2, lower.tail=TRUE) * rucker_mean$SE
		rucker_mean$CI_upp <- rucker_mean$Estimate + qnorm(Qthresh/2, lower.tail=TRUE) * rucker_mean$SE
		rucker_mean$P <- 2 * pt(abs(rucker_mean$Estimate/rucker_mean$SE), nsnp-1, lower.tail=FALSE)

		res <- rbind(rucker$rucker, rucker_point, rucker_mean, rucker_median)
		rownames(res) <- NULL

		p1 <- ggplot2::ggplot(bootstrap, ggplot2::aes_string(x="Q", y="Qdash")) +
			ggplot2::geom_point(ggplot2::aes_string(colour="model")) +
			ggplot2::geom_point(data=subset(bootstrap, i=="Full")) +
			ggplot2::scale_colour_brewer(type="qual") +
			ggplot2::xlim(0, max(bootstrap$Q, bootstrap$Qdash)) +
			ggplot2::ylim(0, max(bootstrap$Q, bootstrap$Qdash)) +
			ggplot2::geom_abline(slope=1, colour="grey") +
			ggplot2::geom_abline(slope=1, intercept=-qchisq(Qthresh, 1, lower.tail=FALSE), linetype="dotted") +
			ggplot2::geom_hline(yintercept = qchisq(Qthresh, nsnp - 2, lower.tail=FALSE), linetype="dotted") +
			ggplot2::geom_vline(xintercept = qchisq(Qthresh, nsnp - 1, lower.tail=FALSE), linetype="dotted") +
			ggplot2::labs(x="Q", y="Q'")

		modsel$model_name <- "IVW"
		modsel$model_name[modsel$model %in% c("C", "D")] <- "Egger"

		p2 <- ggplot2::ggplot(modsel, ggplot2::aes_string(x="Estimate")) +
			ggplot2::geom_density(ggplot2::aes_string(fill="model_name"), alpha=0.4) +
			ggplot2::geom_vline(data=res, ggplot2::aes_string(xintercept="Estimate", colour="Method")) +
			ggplot2::scale_colour_brewer(type="qual") +
			ggplot2::scale_fill_brewer(type="qual") + 
			ggplot2::labs(fill="Bootstrap estimates", colour="")

		return(list(rucker=rucker, res=res, bootstrap_estimates=modsel, boostrap_q=bootstrap, q_plot=p1, e_plot=p2))
	}
}



#' MR Rucker with outliers automatically detected and removed
#'
#' Uses Cook's distance D > 4/nsnp to iteratively remove outliers.
#'
#' @md
#' @param dat Output from [`harmonise_data`].
#' @param parameters List of parameters. The default is `default_parameters()`.
#'
#' @return List
#' @export
mr_rucker_cooksdistance <- function(dat, parameters=default_parameters())
{

	if("mr_keep" %in% names(dat)) dat <- subset(dat, mr_keep)

	dat_orig <- dat
	rucker_orig <- mr_rucker(dat_orig, parameters)
	rucker <- rucker_orig
	cooks_threshold <- 4/nrow(dat)
	index <- rucker_orig$cooksdistance > cooks_threshold

	i <- 1
	l <- list()
	while(any(index) & sum(!index) > 3)
	{
		dat <- dat[!index, ]
		cooks_threshold <- 4/nrow(dat)
		rucker <- mr_rucker(dat, parameters)
		l[[i]] <- rucker
		index <- rucker$cooksdistance > cooks_threshold
		i <- i + 1
	}
	
	rucker$removed_snps <- dat_orig$SNP[! dat_orig$SNP %in% dat$SNP]
	rucker$selected$Method <- "Rucker (CD)"
	rucker$rucker$Method <- paste0(rucker$rucker$Method, " (CD)")
	return(rucker)
}
