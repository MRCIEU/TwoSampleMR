#' Perform Rucker, Median and Mode
#'
#' @param dat Output from harmonise_data
#' @param parameters default_parameters()
#' @param methods Named methods to use. Options are "rucker", "rucker jackknife", "rucker bootstrap", "rucker cooksdistance", "mode", "median"
#'
#' @export
#' @return Data frame
run_mr <- function(dat, parameters=default_parameters(), methods=c("rucker jackknife", "mode", "median"), plots=FALSE)
{
	dat <- subset(dat, mr_keep)
	d <- subset(dat, !duplicated(paste(id.exposure, " - ", id.outcome)), select=c(exposure, outcome, id.exposure, id.outcome))
	res <- list()
	for(j in 1:nrow(d))
	{
		res[[j]] <- list()
		x <- subset(dat, exposure == d$exposure[j] & outcome == d$outcome[j])
		message(x$exposure[1], " - ", x$outcome[1])
		x <- dplyr::select(x, SNP, a1=effect_allele.exposure, a2=other_allele.exposure, beta.exposure, se.exposure, pval.exposure, beta.outcome, se.outcome, pval.outcome, samplesize.exposure, samplesize.outcome, mr_keep)
		res[[j]]$exposure <- x$exposure[1]
		res[[j]]$outcome <- x$outcome[1]
		res[[j]]$id.exposure <- x$id.exposure[1]
		res[[j]]$id.outcome <- x$id.outcome[1]

		if(nrow(x) == 1)
		{
			a <- mr_wald_ratio(x$beta.exposure, x$beta.outcome, x$se.exposure, x$se.outcome)
			out <- tibble(
				# exposure = x$exposure[1],
				# outcome = x$outcome[1],
				# id.exposure = x$id.exposure[1],
				# id.outcome = x$id.outcome[1],
				Method = "Wald ratio",
				Estimate = a$b,
				SE = a$se,
				CI_low = a$b - 1.96 * a$se,
				CI_upp = a$b + 1.96 * a$se,
				P = a$pval,
				nsnp = 1
			)
			res[[j]]$out <- out
		} else if(nrow(x) <= 3) {
			a <- mr_ivw(x$beta.exposure, x$beta.outcome, x$se.exposure, x$se.outcome)
			out <- tibble(
				# exposure = x$exposure[1],
				# outcome = x$outcome[1],
				# id.exposure = x$id.exposure[1],
				# id.outcome = x$id.outcome[1],
				Method = "IVW fixed effects",
				Estimate = a$b,
				SE = a$se,
				CI_low = a$b - 1.96 * a$se,
				CI_upp = a$b + 1.96 * a$se,
				P = a$pval,
				nsnp = nrow(x)
			)
			res[[j]]$out <- out
		} else {
			l <- list()
			r <- list()
			p <- list()
			i <- 1
			if("rucker" %in% methods & !"rucker jackknife" %in% methods & FALSE)
			{
				message("r")
				temp <- try(mr_rucker_internal(x, parameters)$rucker)
				if(class(temp) != "try-error")
				{
					l[[i]] <- temp
					i <- i + 1
				}
			}
			if("rucker jackknife" %in% methods | TRUE)
			{
				message("rj")
				temp <- try(mr_rucker_jackknife_internal(x, parameters))
				if(class(temp) != "try-error")
				{

					l[[i]] <- temp$res
					p$rucker1 <- temp$q_plot
					p$rucker2 <- temp$e_plot
					r$intercept <- temp$rucker$intercept
					r$Q <- temp$Q
					i <- i + 1
				}
			}
			if("median" %in% methods | TRUE)
			{
				message("med")
				temp <- try(mr_median(x, parameters))
				if(class(temp) != "try-error")
				{
					l[[i]] <- rename_result_cols(temp)
					i <- i + 1
				}
			}
			if("mode" %in% methods | TRUE)
			{
				message("mod")
				temp <- try(mr_mode(x, parameters))
				if(class(temp) != "try-error")
				{
					l[[i]] <- rename_result_cols(temp)
					i <- i + 1
				}
			}
			if(length(l) > 0)
			{
				out <- suppressWarnings(dplyr::bind_rows(l))
				# out$exposure <- x$exposure[1]
				# out$outcome <- x$outcome[1]
				# out$id.exposure <- x$id.exposure[1]
				# out$id.outcome <- x$id.outcome[1]

				ind <- sign(x$beta.exposure) == -1
				x$beta.exposure <- abs(x$beta.exposure)
				x$beta.outcome[ind] <- x$beta.outcome[ind] * -1

				temp <- merge(out, dplyr::select(r$intercept, Method=Method, intercept=Estimate), by="Method")
				temp$intercept[is.na(temp$intercept)] <- 0

				temp <- subset(temp, Method %in% c("IVW random effects", "Egger random effects", "Weighted mode", "Weighted median"))


				res[[j]]$out <- out
				res[[j]]$r <- r
				if(plots)
				{
					p$sc <- ggplot2::ggplot(x, ggplot2::aes(x=beta.exposure, y=beta.outcome)) +
						ggplot2::geom_errorbar(ggplot2::aes(ymin=beta.outcome-se.outcome * 1.96, ymax=beta.outcome+se.outcome * 1.96), width=0, colour="grey") +
						ggplot2::geom_errorbarh(ggplot2::aes(xmin=beta.exposure-se.exposure * 1.96, xmax=beta.exposure+se.exposure * 1.96), height=0, colour="grey") +
						ggplot2::geom_point() +
						ggplot2::geom_abline(data=temp, ggplot2::aes(slope=Estimate, intercept=intercept, colour=Method)) +
						ggplot2::scale_colour_brewer(type="qual") +
						ggplot2::labs(x=res[[j]]$exposure, y=res[[j]]$outcome)
					res[[j]]$p <- p
				}
			}
		}
		res[[j]]$isq <- Isq(x$beta.exposure, x$se.exposure)
		res[[j]]$dat <- as.data.frame(x)
	}

	# temp <- res
	attributes(res) <- d

	# ggplot(dat, aes(x=))

	return(res)
}



rename_result_cols <- function(x)
{
	names(x)[names(x) == "method"] <- "Method"
	names(x)[names(x) == "b"] <- "Estimate"
	names(x)[names(x) == "se"] <- "SE"
	names(x)[names(x) == "ci_low"] <- "CI_low"
	names(x)[names(x) == "ci_upp"] <- "CI_upp"
	names(x)[names(x) == "pval"] <- "P"
	return(x)
}
