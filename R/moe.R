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

system_metrics <- function(dat)
{
	library(car)

	# Number of SNPs
	# Sample size outcome
	# Sample size exposure
	metrics <- list()
	metrics$nsnp <- nrow(dat)
	metrics$nout <- mean(dat$samplesize.outcome, na.rm=TRUE)
	metrics$nexp <- mean(dat$samplesize.exposure, na.rm=TRUE)

	# F stats
	Fstat <- qf(dat$pval.exposure, 1, dat$samplesize.exposure, lower.tail=FALSE)
	Fstat[is.infinite(Fstat)] <- 300
	metrics$meanF <- mean(Fstat, na.rm=TRUE)
	metrics$varF <- var(Fstat, na.rm=TRUE)
	metrics$medianF <- median(Fstat, na.rm=TRUE)

	# IF more than 1 SNP

	if(nrow(dat) > 1)
	{
		# Egger-Isq
		metrics$egger_isq <- Isq(dat$beta.exposure, dat$se.exposure)
	}

	if(nrow(dat) > 2)
	{	
		# IF more than 2 SNP
		ruck <- mr_rucker(dat)[[1]]

		# Q_ivw
		# Q_egger
		# Q_diff
		metrics$Isq <- (ruck$Q$Q[1] - (ruck$Q$df[1]-1))/ruck$Q$Q[1]
		metrics$Isqe <- (ruck$Q$Q[2] - (ruck$Q$df[2]-1))/ruck$Q$Q[2]
		metrics$Qdiff <- ruck$Q$Q[3]

		# Intercept / se
		metrics$intercept <- abs(ruck$intercept$Estimate[1]) / ruck$intercept$SE[1]

		# Influential outliers
		dfbeta_thresh <- 2 * nrow(dat)^-0.5
		cooksthresh1 <- 4 / (nrow(dat) - 2)
		cooksthresh2 <- 4 / (nrow(dat) - 3)
		inf1 <- influence.measures(ruck$lmod_ivw)$infmat
		inf2 <- influence.measures(ruck$lmod_egger)$infmat
		metrics$dfb1_ivw <- sum(inf1[,1] > dfbeta_thresh) / nrow(dat)
		metrics$dfb2_ivw <- sum(inf1[,2] > dfbeta_thresh) / nrow(dat)
		metrics$dfb3_ivw <- sum(inf1[,3] > dfbeta_thresh) / nrow(dat)
		metrics$cooks_ivw <- sum(inf1[,4] > cooksthresh1) / nrow(dat)
		metrics$dfb1_egger <- sum(inf2[,1] > dfbeta_thresh) / nrow(dat)
		metrics$dfb2_egger <- sum(inf2[,2] > dfbeta_thresh) / nrow(dat)
		metrics$dfb3_egger <- sum(inf2[,3] > dfbeta_thresh) / nrow(dat)
		metrics$cooks_egger <- sum(inf2[,4] > cooksthresh2) / nrow(dat)

		# Homoscedasticity
		metrics$homosc_ivw <- car::ncvTest(ruck$lmod_ivw)$ChiSquare
		metrics$homosc_egg <- car::ncvTest(ruck$lmod_egger)$ChiSquare

		# Normality of residuals
		metrics$shap_ivw <- shapiro.test(residuals(ruck$lmod_ivw))$statistic
		metrics$shap_egger <- shapiro.test(residuals(ruck$lmod_egger))$statistic
		metrics$ks_ivw <- ks.test(residuals(ruck$lmod_ivw), "pnorm")$statistic
		metrics$ks_egger <- ks.test(residuals(ruck$lmod_egger), "pnorm")$statistic

	}
	return(metrics)
}

get_metrics <- function(dat)
{
	metrics <- system_metrics(dat)
	# Steiger
	steiger_keep <- dat$steiger_dir
	metrics$st_correct <- sum(steiger_keep) / nrow(dat)
	metrics$st_unknown <- sum(dat$steiger_pval < 0.05) / nrow(dat)
	metrics$st_incorrect <-  sum(!dat$steiger_dir & dat$steiger_pval < 0.05) / nrow(dat)

	dat2 <- dat[steiger_keep, ]
	if(nrow(dat2) > 0)
	{
		metrics2 <- system_metrics(dat2)
		names(metrics2) <- paste0(names(metrics2), "_after_steiger")
		metrics <- c(metrics, metrics2)
	}
	return(metrics)
}

get_rf_method <- function(m1, m2, ds, rf)
{
	nom <- names(rf)
	l <- list()
		res1 <- m1$out
		res1$Method <- paste0(res1$Method, " - tophits")
		res2 <- m2$out
		res2$Method <- paste0(res2$Method, " - steiger")
		res <- dplyr::bind_rows(res1, res2)
		if(nrow(m2$dat) > 5)
		{
			message("RF")
			met <- as.data.frame(get_metrics(ds))
			pr <- rep(0, length(rf))
			for(i in 1:length(pr))
			{
				pr[i] <- predict(rf[[i]], met)
			}
			pr <- data.frame(pr=pr, Method=nom)
			res <- dplyr::left_join(res, pr, by="Method")
			ress <- res[which.max(res$pr)[1], ]
			# n[[j]] <- data.frame(id.exposure=res$id.exposure[1], id.outcome=res$id.outcome[1], exposure=res$exposure[1], outcome=res$outcome[1], selected_method=ress$Method)
			ress$Method <- "RF"
			res <- rbind(res, ress)
		}

	return(res)
}

get_rsq <- function(dat)
{
	stopifnot(length(unique(dat$exposure)) == 1)
	stopifnot(length(unique(dat$outcome)) == 1)
	stopifnot(length(unique(dat$units.exposure)) == 1)
	stopifnot(length(unique(dat$units.outcome)) == 1)


	dat$pval.exposure[dat$pval.exposure < 1e-300] <- 1e-300
	if(dat$units.exposure[1] == "log odds")
	{
		ind1 <- !is.na(dat$beta.exposure) &
			!is.na(dat$eaf.exposure) &
			!is.na(dat$ncase.exposure) &
			!is.na(dat$ncontrol.exposure)
		dat$rsq.exposure <- NA
		if(sum(ind1) > 0)
		{
			dat$rsq.exposure[ind1] <- get_r_from_lor(
				dat$beta.exposure[ind1],
				dat$eaf.exposure[ind1],
				dat$ncase.exposure[ind1],
				dat$ncontrol.exposure[ind1],
				0.1
			)^2
		}
	} else {
		ind1 <- !is.na(dat$pval.exposure) & !is.na(dat$samplesize.exposure)
		dat$rsq.exposure <- NA
		if(sum(ind1) > 0)
		{		
			dat$rsq.exposure[ind1] <- get_r_from_pn(
				dat$pval.exposure[ind1],
				dat$samplesize.exposure[ind1]
			)
		}
	}


	dat$pval.outcome[dat$pval.outcome < 1e-300] <- 1e-300
	if(dat$units.outcome[1] == "log odds")
	{
		ind1 <- !is.na(dat$beta.outcome) &
			!is.na(dat$eaf.outcome) &
			!is.na(dat$ncase.outcome) &
			!is.na(dat$ncontrol.outcome)
		dat$rsq.outcome <- NA
		if(sum(ind1) > 0)
		{
			dat$rsq.outcome[ind1] <- get_r_from_lor(
				dat$beta.outcome[ind1],
				dat$eaf.outcome[ind1],
				dat$ncase.outcome[ind1],
				dat$ncontrol.outcome[ind1],
				0.1
			)^2
		}
	} else {
		ind1 <- !is.na(dat$pval.outcome) & !is.na(dat$samplesize.outcome)
		dat$rsq.outcome <- NA
		if(sum(ind1) > 0)
		{		
			dat$rsq.outcome[ind1] <- get_r_from_pn(
				dat$pval.outcome[ind1],
				dat$samplesize.outcome[ind1]
			)
		}
	}
	return(dat)
}



#' Mixture of experts
#'
#' @param dat Output from harmonise_data. Ensure that for continuous traits the p-value and sample size are not missing, and for binary traits the number of cases, number of controls and allele frequencies are provided. Must also ensure that the units are provided, with binary traits in units of 'log odds'.
#' @param rf The trained random forest for the methods. This is available to download at https://www.dropbox.com/s/k0grrhh0ak8er7q/rf.rdata?dl=0
#'
#' @export
#' @return List
mr_moe <- function(dat, rf)
{
	require(randomForest)
	dat <- suppressMessages(get_rsq(dat))
	st <- psych::r.test(
		n = dat$samplesize.exposure, 
		n2 = dat$samplesize.outcome, 
		r12 = sqrt(dat$rsq.exposure), 
		r34 = sqrt(dat$rsq.outcome)
	)

	dat$steiger_dir <- dat$rsq.exposure > dat$rsq.outcome
	dat$steiger_pval <- st$p

	dat <- subset(dat, mr_keep)
	d <- subset(dat, !duplicated(paste(id.exposure, " - ", id.outcome)), 
		select=c(exposure, outcome, id.exposure, id.outcome)
	)
	res <- list()

	for(j in 1:nrow(d))
	{
		res[[j]] <- list()		
		x <- subset(dat, exposure == d$exposure[j] & outcome == d$outcome[j])
		res[[j]]$exposure <- x$exposure[1]
		res[[j]]$outcome <- x$outcome[1]
		res[[j]]$id.exposure <- x$id.exposure[1]
		res[[j]]$id.outcome <- x$id.outcome[1]

		message(d$exposure[j], " - ", d$outcome[j])
		message("Basic methods")
		m1 <- suppressMessages(run_mr(subset(x, id.exposure != id.outcome))[[1]])
		message("Applying Steiger filtering")
		m2 <- suppressMessages(run_mr(subset(x, id.exposure != id.outcome & steiger_dir))[[1]])
		message("Mixture of experts")
		m3 <- get_rf_method(m1, m2, x, rf)
		res[[j]]$m1 <- m1
		res[[j]]$m2 <- m2
		res[[j]]$m3 <- m3
	}

	attributes(res) <- d

	return(res)

}