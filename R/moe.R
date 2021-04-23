#' @importFrom stats influence.measures ks.test median pnorm residuals shapiro.test var
system_metrics <- function(dat)
{
	# Number of SNPs
	# Sample size outcome
	# Sample size exposure
	metrics <- list()
	metrics$nsnp <- nrow(dat)
	metrics$nout <- mean(dat$samplesize.outcome, na.rm=TRUE)
	metrics$nexp <- mean(dat$samplesize.exposure, na.rm=TRUE)

	# F stats
	# Fstat <- qf(dat$pval.exposure, 1, dat$samplesize.exposure, lower.tail=FALSE)
	Fstat <- dat$beta.exposure^2 / dat$se.exposure^2
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
		sct <- mr_sign(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome)
		metrics$sct <- -log10(sct$pval) * sign(sct$b)
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
#' Based on the method described here \url{https://www.biorxiv.org/content/early/2017/08/23/173682}.
#' Once all MR methods have been applied to a summary set, you can then use the mixture of experts to predict the method most likely to be the most accurate.
#'
#' @md
#' @param res Output from [`mr_wrapper`]. 
#' @param rf The trained random forest for the methods. This is available to download at <https://www.dropbox.com/s/5la7y38od95swcf/rf.rdata?dl=0>.
#' 
#' @md
#' @details
#' The `mr_moe` function modifies the `estimates` item in the list of results from the [`mr_wrapper`] function. It does three things:
#' 1. Adds the MOE column, which is a predictor for each method for how well it performs in terms of high power and low type 1 error (scaled 0-1, where 1 is best performance). 
#' 2. It renames the methods to be the estimating method + the instrument selection method. There are 4 instrument selection methods: Tophits (i.e. no filtering), directional filtering (DF, an unthresholded version of Steiger filtering), heterogeneity filtering (HF, removing instruments that make substantial (p < 0.05) contributions to Cochran's Q statistic), and DF + HF which is where DF is applied and the HF applied on top of that. 
#' 3. It orders the table to be in order of best performing method.
#' 
#' Note that the mixture of experts has only been trained on datasets with at least 5 SNPs. If your dataset has fewer than 5 SNPs this function might return errors.
#' 
#' @export
#' @return List
#' @examples
#' \dontrun{
#' # Load libraries
#' library(dplyr)
#' library(randomForest)
#' library(car)
#' 
#' # Example of body mass index on coronary heart disease
#' # Extract and harmonise data
#' a <- extract_instruments(2)
#' b <- extract_outcome_data(a$SNP, 7)
#' dat <- harmonise_data(a,b)
#' 
#' # Apply all MR methods
#' r <- mr_wrapper(dat)
#' 
#' # Load the rf object containing the trained models
#' load("rf.rdata")
#' # Update the results with mixture of experts
#' r <- mr_moe(r, rf)
#' 
#' # Now you can view the estimates, and see that they have 
#' # been sorted in order from most likely to least likely to 
#' # be accurate, based on MOE prediction
#' r[[1]]$estimates
#'}
mr_moe <- function(res, rf)
{
	requireNamespace("dplyr", quietly = TRUE)
	requireNamespace("randomForest", quietly = TRUE)
	lapply(res, function(x)
	{
		o <- try(mr_moe_single(x, rf))
		if(class(o) == "try-error")
		{
			return(x)
		} else {
			return(o)
		}
	})
}

mr_moe_single <- function(res, rf)
{
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("randomForest", quietly = TRUE)
	metric <- res$info[1,] %>% dplyr::select(-c(id.exposure, id.outcome, steiger_filtered, outlier_filtered, nsnp_removed))

	methodlist <- names(rf)
	pred <- lapply(methodlist, function(m)
	{
		d <- dplyr::tibble(
			method = m,
			MOE = predict(rf[[m]], metric, type="prob")[,2]
		)
		return(d)
	}) %>% bind_rows %>% arrange(desc(MOE))
	if("MOE" %in% names(res$estimates))
	{
		message("Overwriting previous MOE estimate")
		res$estimates <- subset(res$estimates, select=-c(MOE, method2))
	}
	res$estimates$selection <- "DF + HF"
	res$estimates$selection[!res$estimates$outlier_filtered & res$estimates$steiger_filtered] <- "DF"
	res$estimates$selection[res$estimates$outlier_filtered & !res$estimates$steiger_filtered] <- "HF"
	res$estimates$selection[!res$estimates$outlier_filtered & !res$estimates$steiger_filtered] <- "Tophits"
	res$estimates$method2 <- paste(res$estimates$method, "-", res$estimates$selection)
	res$estimates <- dplyr::left_join(res$estimates, pred, by=c("method2"="method")) %>% dplyr::arrange(dplyr::desc(MOE))
	return(res)
}
