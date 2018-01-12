

#' Create scatter plot with lines showing the causal estimate for different MR tests
#'
#' Requires dev version of ggplot2
#' 
#' @param mr_results Output from \code{mr}
#' @param dat Output from \code{harmonise_exposure_outcome}
#' @export
#' @return List of plots
mr_scatter_plot <- function(mr_results, dat)
{
	# dat <- subset(dat, paste(id.outcome, id.exposure) %in% paste(mr_results$id.outcome, mr_results$id.exposure))
	requireNamespace("ggplot2", quietly=TRUE)
	requireNamespace("plyr", quietly=TRUE)
	mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"), function(d)
	{
		d <- plyr::mutate(d)
		if(nrow(d) < 2 | sum(d$mr_keep) == 0)
		{
			return(blank_plot("Insufficient number of SNPs"))
		}
		d <- subset(d, mr_keep)
		index <- d$beta.exposure < 0
		d$beta.exposure[index] <- d$beta.exposure[index] * -1
		d$beta.outcome[index] <- d$beta.outcome[index] * -1
		mrres <- subset(mr_results, id.exposure == d$id.exposure[1] & id.outcome == d$id.outcome[1])
		mrres$a <- 0
		if("MR Egger" %in% mrres$method)
		{
			temp <- mr_egger_regression(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
			mrres$a[mrres$method == "MR Egger"] <- temp$b_i
		}

		if("MR Egger (bootstrap)" %in% mrres$method)
		{
			temp <- mr_egger_regression_bootstrap(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
			mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
		}

		ggplot2::ggplot(data=d, ggplot2::aes(x=beta.exposure, y=beta.outcome)) +
			ggplot2::geom_errorbar(ggplot2::aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="grey", width=0) +
			ggplot2::geom_errorbarh(ggplot2::aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="grey", height=0) +
			ggplot2::geom_point(ggplot2::aes(text=paste("SNP:", SNP))) +
			ggplot2::geom_abline(data=mrres, ggplot2::aes(intercept=a, slope=b, colour=method), show.legend=TRUE) +
			ggplot2::scale_colour_manual(values=c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")) +
			ggplot2::labs(colour="MR Test", x=paste("SNP effect on", d$exposure[1]), y=paste("SNP effect on", d$outcome[1])) +
			ggplot2::theme(legend.position="top", legend.direction="vertical") +
			ggplot2::guides(colour=ggplot2::guide_legend(ncol=2))
	})
	mrres
}


#' Plot results from leaveoneout analysis
#'
#' @param leaveoneout_results Output from \code{mr_leaveoneout}
#'
#' @export
#' @return List of plots
mr_leaveoneout_plot <- function(leaveoneout_results)
{
	requireNamespace("ggplot2", quietly=TRUE)
	requireNamespace("plyr", quietly=TRUE)
	res <- plyr::dlply(leaveoneout_results, c("id.exposure", "id.outcome"), function(d)
	{
		d <- plyr::mutate(d)
		# Need to have at least 3 SNPs because IVW etc methods can't be performed with fewer than 2 SNPs
		if(sum(!grepl("All", d$SNP)) < 3) {
			return(
				blank_plot("Insufficient number of SNPs")
			)
		}
		d$up <- d$b + 1.96 * d$se
		d$lo <- d$b - 1.96 * d$se
		d$tot <- 1
		d$tot[d$SNP != "All"] <- 0.01
		d$SNP <- as.character(d$SNP)
		nom <- d$SNP[d$SNP != "All"]
		nom <- nom[order(d$b)]
		d <- rbind(d, d[nrow(d),])
		d$SNP[nrow(d)-1] <- ""
		d$b[nrow(d)-1] <- NA
		d$up[nrow(d)-1] <- NA
		d$lo[nrow(d)-1] <- NA
		d$SNP <- ordered(d$SNP, levels=c("All", "", nom))

		ggplot2::ggplot(d, ggplot2::aes(y=SNP, x=b)) +
		ggplot2::geom_vline(xintercept=0, linetype="dotted") +
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
		ggplot2::labs(y="", x=paste0("MR leave-one-out sensitivity analysis for\n'", d$exposure[1], "' on '", d$outcome[1], "'"))
	})
	res
}


blank_plot <- function(message)
{
	requireNamespace("ggplot2", quietly=TRUE)
	ggplot2::ggplot(data.frame(a=0,b=0,n=message)) + ggplot2::geom_text(ggplot2::aes(x=a,y=b,label=n)) + ggplot2::labs(x=NULL,y=NULL) + ggplot2::theme(axis.text=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank())
}


#' Forest plot
#'
#' @param singlesnp_results from \code{mr_singlesnp}
#' @param exponentiate Plot on exponential scale. Default=FALSE
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
#' @param singlesnp_results from \code{mr_singlesnp}
#' @param mr_results Results from \code{mr}
#' @param exponentiate Plot on exponential scale. Default=FALSE
#' @param bandwidth Density bandwidth parameter
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
			ggplot2::geom_point(y=0, colour="red", ggplot2::aes(size=1/se, text=paste("SNP:", SNP))) +
			ggplot2::geom_vline(data=mr_results, aes(xintercept=b, colour=method)) +
			ggplot2::scale_colour_brewer(type="qual") +
			ggplot2::labs(x="Per SNP MR estiamte")
	})
	res
}

#' Funnel plot
#'
#' Create funnel plot from single SNP analyses
#'
#' @param singlesnp_results from \code{mr_singlesnp}
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

