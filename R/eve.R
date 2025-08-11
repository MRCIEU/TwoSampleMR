mr_mean_ivw <- function(d) {
	d <- subset(d, mr_keep)
	stopifnot(nrow(d) >= 1)
	b_exp <- d$beta.exposure
	b_out <- d$beta.outcome
	se_exp <- d$se.exposure
	se_out <- d$se.outcome
	ratios <- b_out / b_exp
	stopifnot(length(unique(d$id.exposure)) == 1)
	stopifnot(length(unique(d$id.outcome)) == 1)
	id.exposure = d$id.exposure[1]
	id.outcome = d$id.outcome[1]

	if (nrow(d) == 1) {
		res <- mr_wald_ratio(b_exp, b_out, se_exp, se_out)
		out <- dplyr::tibble(
			id.exposure = id.exposure,
			id.outcome = id.outcome,
			method = "Wald ratio",
			nsnp = nrow(d),
			b=res$b,
			se=res$b,
			ci_low=b-se*1.96,
			ci_upp=b+se*1.96,
			pval=stats::pt(abs(b/se), 1, lower.tail=FALSE) * 2
		)
		return(list(estimates=out))
	}

	unw <- summary(stats::lm(b_out ~ -1 + b_exp))
	unw_out <- dplyr::tibble(
		id.exposure = id.exposure,
		id.outcome = id.outcome,
		method = "Simple mean",
		nsnp = nrow(d),
		b = unw$coefficients[1],
		se = unw$coefficients[2],
		ci_low=b-se*1.96,
		ci_upp=b+se*1.96,
		pval=stats::pt(abs(b/se), nrow(d)-1, lower.tail=FALSE) * 2
	)

	weights1 <- sqrt(b_exp^2 / se_out^2)
	y1 <- ratios * weights1

	ivw1 <- stats::lm(y1 ~ -1 + weights1)
	weights2 <- sqrt(((se_out^2 + stats::coefficients(ivw1)[1]^2 * se_exp^2) / b_exp^2)^-1)

	y2 <- ratios * weights2
	ivw2 <- summary(stats::lm(y2 ~ -1 + weights2))

	ivwoutliers <- dplyr::tibble(id.exposure = id.exposure, id.outcome = id.outcome, SNP=d$SNP, Qj=weights2^2 * (ratios - stats::coefficients(ivw2)[1])^2, Qpval=stats::pchisq(Qj,1,lower.tail=FALSE))

	Qivw2 <- sum(ivwoutliers$Qj)
	Qivw2pval <- stats::pchisq(Qivw2, nrow(d)-1, lower.tail=FALSE)

	# Collate
	re_out <- dplyr::tibble(
		id.exposure = id.exposure,
		id.outcome = id.outcome,
		method = c("RE IVW"),
		nsnp = nrow(d),
		b = c(ivw2$coefficients[1]),
		se = c(ivw2$coefficients[2]),
		ci_low = b - se * 1.96,
		ci_upp = b + se * 1.96,
		pval = stats::pt(abs(b / se), nrow(d)-1, lower.tail=FALSE) * 2
	)

	fe_out <- re_out
	fe_out$se <- fe_out$se / max(1, ivw2$sigma)
	fe_out$pval <- stats::pnorm(abs(fe_out$b / fe_out$se), lower.tail=FALSE) * 2
	fe_out$method <- c("FE IVW")

	out <- dplyr::bind_rows(unw_out, fe_out, re_out)

	# Pleiotropy
	heterogeneity <- dplyr::tibble(
		id.exposure = id.exposure,
		id.outcome = id.outcome,
		method = c("IVW"),
		Q = c(Qivw2),
		df = c(nrow(d)-1),
		pval = stats::pchisq(Q, df, lower.tail=FALSE)
	)

	ret <- list(estimates=out, heterogeneity = heterogeneity, outliers=ivwoutliers)
	return(ret)
}

mr_mean_egger <- function(d) {
	d <- subset(d, mr_keep)
	stopifnot(nrow(d) >= 3)
	b_exp <- d$beta.exposure
	b_out <- d$beta.outcome
	se_exp <- d$se.exposure
	se_out <- d$se.outcome
	stopifnot(length(unique(d$id.exposure)) == 1)
	stopifnot(length(unique(d$id.outcome)) == 1)
	id.exposure = d$id.exposure[1]
	id.outcome = d$id.outcome[1]

	ratios <- b_out / b_exp
	weights1 <- sqrt(b_exp^2 / se_out^2)
	y1 <- ratios * weights1

	# Egger
	egger1 <- stats::lm(y1 ~ weights1)
	weights2 <- sqrt(((se_out^2 + stats::coefficients(egger1)[2]^2 * se_exp^2) / b_exp^2)^-1)

	y2 <- ratios * weights2
	egger2 <- summary(stats::lm(y2 ~ weights2))
	eggeroutliers <- dplyr::tibble(
		SNP=d$SNP,
		Qj = weights2^2 * (ratios - stats::coefficients(egger2)[1,1] / weights2 - stats::coefficients(egger2)[2,1])^2,
		Qpval=stats::pchisq(Qj,1,lower.tail=FALSE)
	)

	Qegger2 <- sum(eggeroutliers$Qj)
	Qegger2pval <- stats::pchisq(Qegger2, nrow(d)-1, lower.tail=FALSE)


	# Collate
	re_out <- dplyr::tibble(
		id.exposure = id.exposure,
		id.outcome = id.outcome,
		method = c("RE Egger"),
		nsnp = nrow(d),
		b = c(egger2$coefficients[2,1]),
		se = c(egger2$coefficients[2,2]),
		ci_low = b - se * 1.96,
		ci_upp = b + se * 1.96,
		pval = stats::pt(abs(b / se), nrow(d)-2, lower.tail=FALSE) * 2
	)

	fe_out <- re_out
	fe_out$se <- fe_out$se / max(1, egger2$sigma)
	fe_out$pval <- stats::pnorm(abs(fe_out$b / fe_out$se), lower.tail=FALSE) * 2
	fe_out$method <- c("FE Egger")

	out <- dplyr::bind_rows(fe_out, re_out)

	# Pleiotropy
	heterogeneity <- dplyr::tibble(
		id.exposure = id.exposure,
		id.outcome = id.outcome,
		method = c("Egger"),
		Q = c(Qegger2),
		df = c(nrow(d)-2),
		pval = stats::pchisq(Q, df, lower.tail=FALSE)
	)

	directional_pleiotropy <- dplyr::tibble(
		id.exposure = id.exposure,
		id.outcome = id.outcome,
		method = c("FE Egger intercept", "RE Egger intercept"),
		nsnp = nrow(d),
		b = c(egger2$coefficients[1,1], egger2$coefficients[1,1]),
		se = c(egger2$coefficients[1,2] / max(1, egger2$sigma), egger2$coefficients[1,2]),
		pval = stats::pt(abs(b/se), nrow(d)-2, lower.tail=FALSE) * 2
	)

	ret <- list(estimates=out, heterogeneity = heterogeneity, directional_pleiotropy = directional_pleiotropy, outliers=eggeroutliers)
	return(ret)
}

mr_mean <- function(dat, parameters=default_parameters()) {
	m1 <- try(mr_mean_ivw(dat))
	m2 <- try(mr_mean_egger(dat))
	if (inherits(m1, "try-error")) {
		return(NULL)
	} else {
		if (inherits(m2, "try-error")) {
			return(m1)
		} else {
			out <- list(
				estimates = dplyr::bind_rows(m1$estimates, m2$estimates),
				heterogeneity = dplyr::bind_rows(m1$heterogeneity, m2$heterogeneity),
				directional_pleiotropy = m2$directional_pleiotropy,
				outliers = m1$outliers
			)
			temp <- dplyr::tibble(
				id.exposure = dat$id.exposure[1],
				id.outcome = dat$id.outcome[1],
				method = "Rucker",
				Q = out$heterogeneity$Q[1] - out$heterogeneity$Q[2],
				df = 1,
				pval = stats::pchisq(Q, df, lower.tail=FALSE)
			)
			out$heterogeneity <- dplyr::bind_rows(out$heterogeneity, temp)
			return(out)
		}
	}
}

mr_all <- function(dat, parameters=default_parameters()) {
	m1 <- mr_mean(dat)
	if (sum(dat$mr_keep) > 3) {
		m2 <- try(mr_median(dat, parameters=parameters))
		m3 <- try(mr_mode(dat, parameters=parameters)[1:3,])
		m1$estimates <- dplyr::bind_rows(m1$estimates, m2, m3)
	}
	m1$info <- c(list(
			id.exposure = dat$id.exposure[1], id.outcome = dat$id.outcome[1]),
			system_metrics(dat)
		) %>% dplyr::as_tibble()
	return(m1)
}

mr_wrapper_single <- function(dat, parameters=default_parameters()) {
	dat <- steiger_filtering(dat)
	m <- list()
	snps_retained <- dplyr::tibble(
		SNP = dat$SNP,
		outlier = FALSE, steiger = FALSE, both = FALSE
	)
	m[[1]] <- mr_all(dat, parameters=parameters)
	if (!is.null(m[[1]])) {
		if ("outliers" %in% names(m[[1]])) {
			temp <- subset(dat, ! SNP %in% subset(m[[1]]$outliers, Qpval < 0.05)$SNP)
			m[[2]] <- mr_all(temp, parameters=parameters)
			snps_retained$outlier[snps_retained$SNP %in% temp$SNP] <- TRUE
		} else {
			m[[2]] <- m[[1]]
			snps_retained$outlier <- TRUE
		}

		dat_st <- subset(dat, steiger_dir)
		snps_retained$steiger[snps_retained$SNP %in% dat_st$SNP] <- TRUE
		if (nrow(dat_st) == 0) {
			m[[3]] <- m[[4]] <- list(
				estimates=dplyr::tibble(method="Steiger null", nsnp = 0, b=0, se=NA, ci_low=NA, ci_upp=NA, pval=1)
			)
		} else {
			m[[3]] <- mr_all(dat_st, parameters=parameters)
			if ("outliers" %in% names(m[[3]])) {
				temp <- subset(dat_st, ! SNP %in% subset(m[[3]]$outliers, Qpval < 0.05)$SNP)
				m[[4]] <- mr_all(temp, parameters=parameters)
				snps_retained$both[snps_retained$SNP %in% temp$SNP] <- TRUE
			} else {
				m[[4]] <- m[[3]]
				snps_retained$both <- TRUE
			}
		}
	}

	if (!is.null(m[[1]])) {
		m[[1]] <- lapply(m[[1]], function(x) {
			x$steiger_filtered <- FALSE
			x$outlier_filtered <- FALSE
			x$id.exposure <- dat$id.exposure[1]
			x$id.outcome <- dat$id.outcome[1]
			return(x)
		})
	}
	if (!is.null(m[[2]])) {
		m[[2]] <- lapply(m[[2]], function(x) {
			x$steiger_filtered <- FALSE
			x$outlier_filtered <- TRUE
			x$id.exposure <- dat$id.exposure[1]
			x$id.outcome <- dat$id.outcome[1]
			return(x)
		})
	}
	if (!is.null(m[[3]])) {
		m[[3]] <- lapply(m[[3]], function(x) {
			x$steiger_filtered <- TRUE
			x$outlier_filtered <- FALSE
			x$id.exposure <- dat$id.exposure[1]
			x$id.outcome <- dat$id.outcome[1]
			return(x)
		})
	}
	if (!is.null(m[[4]])) {
		m[[4]] <- lapply(m[[4]], function(x) {
			x$steiger_filtered <- TRUE
			x$outlier_filtered <- TRUE
			x$id.exposure <- dat$id.exposure[1]
			x$id.outcome <- dat$id.outcome[1]
			return(x)
		})
	}

	nom <- lapply(m, names) %>% unlist %>% unique %>% as.list
	nom <- nom[nom != "outliers"]
	o <- lapply(nom, function(i) {
		lapply(m, function(y) y[[i]]) %>% dplyr::bind_rows()
	})
	names(o) <- nom
	o$info <- o$info %>% dplyr::mutate(nsnp_removed = dplyr::first(nsnp)-nsnp)
	o$snps_retained <- snps_retained

	return(o)
}

#' Perform full set of MR analyses
#'
#'
#' @param dat Output from [harmonise_data()].
#' @param parameters Parameters to pass to MR functions. Output from [default_parameters()] used as default.
#'
#' @export
#' @return list
mr_wrapper <- function(dat, parameters=default_parameters()) {
	plyr::dlply(dat, c("id.exposure", "id.outcome"), function(x) {
		message("Performing MR analysis of '", x$id.exposure[1], "' on '", x$id.outcome[1], "'")
		d <- subset(x, mr_keep)
		o <- mr_wrapper_single(d, parameters=parameters)
		o
	})
}
