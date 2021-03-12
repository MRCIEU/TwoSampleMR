#' Perform 2 sample MR on each SNP individually
#'
#' @md
#' @param dat Output from [`harmonise_data`].
#' @param parameters List of parameters. The default is `default_parameters()`.
#' @param single_method Function to use for MR analysis. The default is `"mr_wald_ratio"`.
#' @param all_method Functions to use for MR analysis. The default is `c("mr_ivw", "mr_egger_regression")`.
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


#' Forest plot
#'
#' @param singlesnp_results from [`mr_singlesnp`].
#' @param exponentiate Plot on exponential scale. The default is `FALSE`.
#'
#' @export
#' @return List of plots
mr_forest_plot <- function(singlesnp_results, exponentiate=FALSE)
{
	requireNamespace("ggplot2", quietly=TRUE)
	requireNamespace("plyr", quietly=TRUE)
	res <- plyr::dlply(singlesnp_results, c("id.exposure", "id.outcome"), function(d)
	{
		d <- plyr::mutate(d)
		if(sum(!grepl("All", d$SNP)) < 2) {
			return(
				blank_plot("Insufficient number of SNPs")
			)
		}
		levels(d$SNP)[levels(d$SNP) == "All - Inverse variance weighted"] <- "All - IVW"
		levels(d$SNP)[levels(d$SNP) == "All - MR Egger"] <- "All - Egger"
		am <- grep("All", d$SNP, value=TRUE)
		d$up <- d$b + 1.96 * d$se
		d$lo <- d$b - 1.96 * d$se
		d$tot <- 0.01
		d$tot[d$SNP %in% am] <- 1
		d$SNP <- as.character(d$SNP)
		nom <- d$SNP[! d$SNP %in% am]
		nom <- nom[order(d$b)]
		d <- rbind(d, d[nrow(d),])
		d$SNP[nrow(d)-1] <- ""
		d$b[nrow(d)-1] <- NA
		d$up[nrow(d)-1] <- NA
		d$lo[nrow(d)-1] <- NA
		d$SNP <- ordered(d$SNP, levels=c(am, "", nom))

		xint <- 0
		if(exponentiate)
		{
			d$b <- exp(d$b)
			d$up <- exp(d$up)
			d$lo <- exp(d$lo)
			xint <- 1
		}

		ggplot2::ggplot(d, ggplot2::aes(y=SNP, x=b)) +
		ggplot2::geom_vline(xintercept=xint, linetype="dotted") +
		# ggplot2::geom_errorbarh(ggplot2::aes(xmin=pmax(lo, min(d$b, na.rm=T)), xmax=pmin(up, max(d$b, na.rm=T)), size=as.factor(tot), colour=as.factor(tot)), height=0) +
		ggplot2::geom_errorbarh(ggplot2::aes(xmin=lo, xmax=up, size=as.factor(tot), colour=as.factor(tot)), height=0) +
		ggplot2::geom_point(ggplot2::aes(colour=as.factor(tot))) +
		ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% "")), colour="grey") +
		ggplot2::scale_colour_manual(values=c("black", "red")) +
		ggplot2::scale_size_manual(values=c(0.3, 1)) +
		# xlim(c(min(c(0, d$b), na.rm=T), max(c(0, d$b), na.rm=T))) +
		ggplot2::theme(
			legend.position="none", 
			axis.text.y=ggplot2::element_text(size=8), 
			axis.ticks.y=ggplot2::element_line(size=0),
			axis.title.x=ggplot2::element_text(size=8)) +
		ggplot2::labs(y="", x=paste0("MR effect size for\n'", d$exposure[1], "' on '", d$outcome[1], "'"))
	})
	res
}


#' Density plot
#'
#' @param singlesnp_results from [`mr_singlesnp`].
#' @param mr_results Results from [`mr`].
#' @param exponentiate Plot on exponentiated scale. The default is `FALSE`.
#' @param bandwidth Density bandwidth parameter.
#'
#' @export
#' @return List of plots
mr_density_plot <- function(singlesnp_results, mr_results, exponentiate=FALSE, bandwidth="nrd0")
{
	requireNamespace("ggplot2", quietly=TRUE)
	requireNamespace("plyr", quietly=TRUE)
	res <- plyr::dlply(singlesnp_results, c("id.exposure", "id.outcome"), function(d)
	{
		d <- plyr::mutate(d)
		if(sum(!grepl("All", d$SNP)) < 2) {
			return(
				blank_plot("Insufficient number of SNPs")
			)
		}
		d$SNP <- as.character(d$SNP)

		d2 <- subset(d, !grepl("All - ", SNP))
		d1 <- subset(mr_results, id.exposure == d2$id.exposure[1] & id.outcome == d2$id.outcome[1])

		xint <- 0
		if(exponentiate)
		{
			d$b <- exp(d$b)
			d$up <- exp(d$up)
			d$lo <- exp(d$lo)
			xint <- 1
		}

		ggplot2::ggplot(d2, ggplot2::aes(x=b)) +
		ggplot2::geom_vline(xintercept=xint, linetype="dotted") +
			ggplot2::geom_density(ggplot2::aes(weight=1/se), bw=bandwidth) +
			ggplot2::geom_point(y=0, colour="red", ggplot2::aes(size=1/se)) +
			ggplot2::geom_vline(data=mr_results, ggplot2::aes(xintercept=b, colour=method)) +
			ggplot2::scale_colour_brewer(type="qual") +
			ggplot2::labs(x="Per SNP MR estimate")
	})
	res
}

#' Funnel plot
#'
#' Create funnel plot from single SNP analyses.
#'
#' @param singlesnp_results from [`mr_singlesnp`].
#'
#' @export
#' @return List of plots
mr_funnel_plot <- function(singlesnp_results)
{
	requireNamespace("ggplot2", quietly=TRUE)
	requireNamespace("plyr", quietly=TRUE)
	res <- plyr::dlply(singlesnp_results, c("id.exposure", "id.outcome"), function(d)
	{
		d <- plyr::mutate(d)
		if(sum(!grepl("All", d$SNP)) < 2) {
			return(
				blank_plot("Insufficient number of SNPs")
			)
		}
		am <- grep("All", d$SNP, value=TRUE)
		d$SNP <- gsub("All - ", "", d$SNP)
		am <- gsub("All - ", "", am)
		ggplot2::ggplot(subset(d, ! SNP %in% am), ggplot2::aes(y = 1/se, x=b)) +
		ggplot2::geom_point() +
		ggplot2::geom_vline(data=subset(d, SNP %in% am), ggplot2::aes(xintercept=b, colour = SNP)) +
		# ggplot2::scale_colour_brewer(type="qual") +
		ggplot2::scale_colour_manual(values = c("#a6cee3", 
                  "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", 
                  "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", 
                  "#6a3d9a", "#ffff99", "#b15928")) +
		ggplot2::labs(y=expression(1/SE[IV]), x=expression(beta[IV]), colour="MR Method") +
		ggplot2::theme(legend.position="top", legend.direction="vertical")
	})
	res
}

